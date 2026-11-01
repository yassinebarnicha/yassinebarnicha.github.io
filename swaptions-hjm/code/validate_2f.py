"""
Local validation of the claims made in the theory text, the page and the report (run it before publishing anything).

Every check recomputes a number quoted in the documentation and compares it with the quoted value and tolerance.
Exit code 0 only if all checks pass.   python validate_2f.py
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from hjm1f.calibration_1f import CalibrationInstrument, CoterminalBootstrap1F
from hjm1f.curve import ZeroCurve
from hjm1f.hjm_1f import HJM1FModel, PiecewiseConstantSigma
from hjm1f.instruments import BermudanCoterminalSwaption, CoterminalSwaption
from hjm1f.market_quotes import forward_swap_rate
from hjm1f.mc_1f import BermudanLSMC1F, MonteCarloSpec
from hjm1f.pde_1f import BermudanPDE1F, PDEGrid
from hjm2f.calibration_2f import CoterminalBootstrap2F
from hjm2f.calibration_global import fit_one_factor, fit_two_factor, from_params
from hjm2f.caps import atm_cap_strike, bachelier_cap_price, model_cap_price, model_cap_vol
from hjm2f.market_data import MarketData, summarise
from hjm2f.hjm_2f import HJM2FModel, normal_grid_2d, psd_sqrt
from hjm2f.mc_2f import BermudanLSMC2F, MonteCarloSpec2F
from hjm2f.quad_2f import BermudanQuad2F

TENORS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
RESULTS = []


def check(name: str, value: float, expected: float, tol: float, relative: bool = False) -> None:
    err = abs(value - expected) / (abs(expected) if relative else 1.0)
    ok = err <= tol
    RESULTS.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name:<70} value {value:>12.6g}  claimed {expected:>10.6g}  err {err:.2e} (tol {tol:.0e}{' rel' if relative else ''})")


def at_most(name: str, value: float, bound: float) -> None:
    ok = value <= bound
    RESULTS.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name:<70} value {value:>12.6g}  bound   {bound:>10.6g}")


def main() -> int:
    t_start = time.time()
    curve = ZeroCurve([0, 1, 2, 3, 5, 7, 10], [1.0, 0.978, 0.954, 0.930, 0.885, 0.842, 0.780])
    ex, pay, acc = [1.0, 2.0, 3.0], [2.0, 3.0, 4.0, 5.0, 6.0, 7.0], [1.0] * 6
    sw = lambda e, payer=True: CoterminalSwaption(e, 7.0, 0.03, pay, acc, payer)
    inst = [CalibrationInstrument(sw(e), v) for e, v in zip(ex, [0.008, 0.0085, 0.009])]
    berm = BermudanCoterminalSwaption(ex, 7.0, 0.03, pay, acc, True)

    print("== 1. No arbitrage (theory 2F, section 2) ==")
    m = HJM2FModel(curve, (0.03, 0.6), PiecewiseConstantSigma(ex, [.009, .010, .011]), PiecewiseConstantSigma(ex, [.006, .005, .007]), -0.6)
    worst = max(abs(m.martingale_defect(t, b, T)) for t, b, T in [(1., 3., 7.), (2., 5., 7.), (3., 4., 7.), (1., 7., 7.)])
    at_most("martingale defect, generic 2F parameters", worst, 1e-12)
    m1 = HJM1FModel(curve, 0.03, PiecewiseConstantSigma(ex, [.009, .010, .011]))
    t, b, T = 2.0, 5.0, 7.0
    mean = m1.mean_under_forward(t, T)
    var = m1.alpha(t)
    d = m1.D(t, b) - m1.D(t, T)
    c_b = 0.5 * (m1.Bbar(b) + m1.Bbar(t)) * var - m1.gamma(t)
    c_T = 0.5 * (m1.Bbar(T) + m1.Bbar(t)) * var - m1.gamma(t)
    lhs = math.exp(-m1.D(t, b) * c_b + m1.D(t, T) * c_T - d * mean + 0.5 * d * d * var) * m1.curve.discount(b) / m1.curve.discount(T)
    check("1F martingale: E^T[P(t,b)/P(t,T)] = P(0,b)/P(0,T)", lhs, m1.curve.discount(b) / m1.curve.discount(T), 1e-12, relative=True)

    print("\n== 2. One-factor code after correction (theory 1F, section 8) ==")
    calib = CoterminalBootstrap1F(curve, 0.03, ex).calibrate(inst)
    model = HJM1FModel(curve, 0.03, calib.sigma_curve)
    worst_ratio = 1.0
    for payer in (True, False):
        for e in ex:
            jam = model.european_swaption_price_jamshidian(sw(e, payer))
            pde = BermudanPDE1F(model, BermudanCoterminalSwaption([e], 7.0, 0.03, pay, acc, payer),
                                PDEGrid.around_factor_std(model, e, 5.0, nx=601, nt=800)).price()
            worst_ratio = min(worst_ratio, pde / jam)
    check("PDE / Jamshidian, payers and receivers, worst ratio (claimed 0.9994)", worst_ratio, 0.9994, 2e-4)
    pde_b = BermudanPDE1F(model, berm, PDEGrid.around_factor_std(model, 3.0, 5.0, nx=801, nt=1500)).price()
    lsmc_b = BermudanLSMC1F(model, berm, MonteCarloSpec(n_paths=200000, seed=1)).price()
    check("1F Bermudan (PDE, implicit) = 0.016741", pde_b, 0.016741, 5e-6)
    check("LSMC within 0.6 % of the PDE (claimed 0.9944)", lsmc_b / pde_b, 0.9944, 3e-3)
    pde_cn = BermudanPDE1F(model, berm, PDEGrid.around_factor_std(model, 3.0, 5.0, nx=801, nt=1500), theta=0.5, rannacher=2).price()
    check("Crank-Nicolson + Rannacher within 0.05 % of the (first-order) implicit price", pde_cn / pde_b, 1.0, 5e-4)

    print("\n== 3. Two factors against one factor (theory 2F, sections 3 and 7) ==")
    s1 = list(calib.sigma_curve.values)
    m0 = HJM2FModel(curve, (0.03, 0.6), PiecewiseConstantSigma(ex, s1), PiecewiseConstantSigma(ex, [0, 0, 0]), -0.6)
    rel = max(abs(m0.european_swaption_price(sw(e), n_grid=161) / model.european_swaption_price_jamshidian(sw(e)) - 1) for e in ex)
    at_most("sigma2 = 0: European 2-D integral vs Jamshidian, relative (claimed 4e-4)", rel, 5e-4)
    a, bb = [0.005, 0.006, 0.007], [0.0036, 0.0036, 0.0036]
    mm = HJM2FModel(curve, (0.03, 0.03), PiecewiseConstantSigma(ex, a), PiecewiseConstantSigma(ex, bb), 1.0)
    m1t = HJM1FModel(curve, 0.03, PiecewiseConstantSigma(ex, [x + y for x, y in zip(a, bb)]))
    rel = max(abs(mm.european_swaption_price(sw(e), n_grid=161) / m1t.european_swaption_price_jamshidian(sw(e)) - 1) for e in ex)
    at_most("rho = 1, lambda1 = lambda2, sigma1+sigma2 = sigma: relative error", rel, 5e-4)
    # grid convergence of the Bermudan quadrature in the one-factor limit (claimed 1 % at 41x61, 0.2 % at 101x141)
    r0 = CoterminalBootstrap2F(curve, (0.03, 0.5), -0.5, 1e-3, ex).calibrate(inst)
    mq = HJM2FModel(curve, (0.03, 0.5), r0.sigma1, r0.sigma2, -0.5)
    check("quadrature 41x61 vs 1F PDE, relative error (claimed 1.2 %)", BermudanQuad2F(mq, berm, n_grid=41, n_quad=61).price() / pde_b - 1, 0.012, 1.5e-3)
    check("quadrature 101x141 vs 1F PDE, relative error (claimed 0.2 %)", BermudanQuad2F(mq, berm, n_grid=101, n_quad=141).price() / pde_b - 1, 0.002, 1.5e-3)
    # quadrature against Monte Carlo on genuine 2F cases (claimed: within 0.7 %)
    worst = 0.0
    for lam2, rho, ratio in [(0.5, -0.5, 0.6), (0.5, -0.9, 0.6), (0.5, 0.5, 0.6), (1.0, -0.7, 1.0)]:
        r = CoterminalBootstrap2F(curve, (0.03, lam2), rho, ratio, ex).calibrate(inst)
        mo = HJM2FModel(curve, (0.03, lam2), r.sigma1, r.sigma2, rho)
        q = BermudanQuad2F(mo, berm, n_grid=81, n_quad=101).price()
        l = BermudanLSMC2F(mo, berm, MonteCarloSpec2F(200000, 7)).price()
        worst = max(worst, abs(l / q - 1))
    at_most("quadrature vs LSMC, worst relative gap over 4 scenarios (claimed 0.4 %)", worst, 0.005)

    print("\n== 4. Caps and smile (theory 2F, sections 4 and 8) ==")
    md = MarketData()
    check("10-year par rate = 0.532 %", md.bermudan(10, [1])[0] * 100, 0.532, 5e-4)
    check("ATM cap strike, 10 years = 0.5257 % (course workbook), within 1 bp", atm_cap_strike(md.curve, 10.0) * 100, 0.5257, 1e-2)
    mq = HJM2FModel(curve, (0.03, 0.5), PiecewiseConstantSigma(ex, [.009, .010, .011]), PiecewiseConstantSigma(ex, [.004, .005, .006]), -0.6)
    k = 1.0 / 1.0145
    put_cf = mq.zcb_put(2.0, 2.5, k)
    # independent check: the same payoff max(k - P(2, 2.5), 0) integrated on a 2-D normal grid under the 2.5-forward measure
    z, w = normal_grid_2d(401)
    x = mq.mean_under_forward(2.0, 2.0) + z @ psd_sqrt(mq.cov(2.0)).T
    pb = np.array([mq.discount_bond(2.0, 2.5, xi) for xi in x])
    put_grid = curve.discount(2.0) * float(np.sum(w * np.maximum(k - pb, 0.0)))
    check("bond put: closed form vs 2-D integral of the same payoff", put_cf, put_grid, 1e-4, relative=True)
    kk = atm_cap_strike(curve, 5.0)
    price = model_cap_price(mq, 5.0, kk)
    at_most("cap price -> implied vol -> cap price round trip", abs(bachelier_cap_price(curve, 5.0, kk, model_cap_vol(mq, 5.0)) - price), 1e-10)
    check("smile ratio at the money within 1.5 % of 1 (kernel mixes smiles whose ATM point differs from the matrix)", md.smile.ratio(1.0, 9.0, 0.0), 1.0, 1.5e-2)
    ends, ins = md.bucket_plan(10, [1, 2, 3, 4, 5], use_smile=True)
    check("bucket plan: 11 volatility segments up to 5 years", float(len(ends)), 11.0, 0.0)
    atm_q = md.coterminal_instruments(10, [1, 2, 3, 4, 5], use_smile=False)
    sm_q = md.coterminal_instruments(10, [1, 2, 3, 4, 5], use_smile=True)
    strike = md.bermudan(10, [1])[0]
    gaps = [(strike - forward_swap_rate(md.curve, float(e), 10.0, [float(i) for i in range(e + 1, 11)], [1.0] * (10 - e))) * 1e4 for e in range(1, 6)]
    check("co-terminal strike minus forward, 1y: -7.6 bp", gaps[0], -7.6, 0.1)
    check("co-terminal strike minus forward, 5y: -52.6 bp", gaps[4], -52.6, 0.1)
    check("co-terminal quote 1y is left unchanged by the smile: 60.2 bp", sm_q[0].normal_vol * 1e4, 60.2, 0.1)
    check("co-terminal quote 2y moved by the smile: 61.7 -> 61.0 bp", sm_q[1].normal_vol * 1e4, 61.0, 0.1)
    check("co-terminal quote 5y moved by the smile: 68.6 -> 67.3 bp", sm_q[4].normal_vol * 1e4, 67.3, 0.1)
    at_most("quotes only move by up to 1.4 bp", max(abs(a.normal_vol - b.normal_vol) for a, b in zip(atm_q, sm_q)) * 1e4, 1.4)

    print("\n== 5. The whole market, 5-year window (page tab 'Fit to the market'; theory 2F, section 4) ==")
    exf = [1.0, 2.0, 3.0, 4.0, 5.0]
    bm = BermudanCoterminalSwaption(exf, 10.0, strike, [float(i) for i in range(1, 11)], [1.0] * 10, True)
    both = md.targets(5.0, TENORS, with_caps=True)
    check("143 ATM swaptions fitted", float(len(both.swaptions)), 143.0, 0.0)
    check("6 ATM caps fitted", float(len(both.caps)), 6.0, 0.0)
    lam1, rmse1, _ = fit_one_factor(md.curve, ends, ins, both)
    check("1F: lambda = -0.039 (swaptions and caps)", lam1, -0.0389, 2e-3)
    c1 = CoterminalBootstrap1F(md.curve, lam1, ends).calibrate(ins)
    zero = PiecewiseConstantSigma(ends, [0.0] * len(ends))
    s1 = summarise(md.report_card(HJM2FModel(md.curve, (lam1, 0.5), c1.sigma_curve, zero, 0.0), 5.0, TENORS, True))
    check("1F: fitted swaptions 12.4 bp", s1["swaptions_fitted"], 12.39, 0.1)
    check("1F: fitted caps 21.7 bp", s1["caps_fitted"], 21.73, 0.2)
    check("1F: swaptions beyond the window 21.4 bp", s1["swaptions_beyond"], 21.36, 0.2)
    check("1F: caps beyond the window 7.6 bp", s1["caps_beyond"], 7.56, 0.2)
    m1 = HJM1FModel(md.curve, lam1, c1.sigma_curve)
    b1 = BermudanPDE1F(m1, bm, PDEGrid.around_factor_std(m1, 5.0, 5.0, nx=801, nt=1500), theta=0.5, rannacher=2).price()
    eur1 = max(p for p, i in zip(c1.model_prices, ins) if i.swaption.expiry in exf)
    check("1F Bermudan = 0.050057", b1, 0.050057, 3e-5)
    check("1F exercise premium = 63.6 bp", (b1 - eur1) * 1e4, 63.6, 0.3)

    lam, rho, ratio = (-0.003414, 0.533868), -0.95, 1.483256
    f = fit_two_factor(md.curve, ends, ins, both, from_params(lam, rho, ratio), n_grid=61, max_fev=60)
    at_most("2F: a local search around the claimed optimum finds nothing below 3.85 bp", 3.85 - f.rmse_bp, 0.0)
    at_most("2F: ... and does not drift above 4.1 bp", f.rmse_bp, 4.1)
    at_most("2F: rho stays at the bound -0.95", abs(f.rho) - 0.9500001, 0.0)
    c2 = CoterminalBootstrap2F(md.curve, lam, rho, ratio, ends).calibrate(ins, n_grid=161)
    m2 = HJM2FModel(md.curve, lam, c2.sigma1, c2.sigma2, rho)
    s2 = summarise(md.report_card(m2, 5.0, TENORS, True))
    check("2F: fitted swaptions 3.9 bp", s2["swaptions_fitted"], 3.92, 0.1)
    check("2F: fitted caps 6.0 bp", s2["caps_fitted"], 5.96, 0.2)
    check("2F: swaptions beyond the window 11.8 bp", s2["swaptions_beyond"], 11.83, 0.2)
    check("2F: caps beyond the window 9.6 bp", s2["caps_beyond"], 9.59, 0.2)
    b2 = BermudanQuad2F(m2, bm, n_grid=81, n_quad=101).price()
    eur2 = max(p for p, i in zip(c2.model_prices, ins) if i.swaption.expiry in exf)
    l2 = BermudanLSMC2F(m2, bm, MonteCarloSpec2F(200000, 11)).price()
    check("2F Bermudan = 0.048161", b2, 0.048161, 3e-5)
    check("2F exercise premium = 44.6 bp", (b2 - eur2) * 1e4, 44.6, 0.3)
    check("2F Bermudan vs 1F = -3.8 %", (b2 / b1 - 1) * 100, -3.8, 0.1)
    check("2F LSMC within 0.6 % of the quadrature (below)", l2 / b2, 0.995, 6e-3)
    check("RMSE ratio on the fitted swaptions, 1F / 2F = 3.2", s1["swaptions_fitted"] / s2["swaptions_fitted"], 3.2, 0.1)

    print("\n== 6. Same Europeans, different Bermudan: rho sweep on the real case ==")
    prem, bermu = {}, {}
    for rr in (-0.95, 0.9):
        c = CoterminalBootstrap2F(md.curve, lam, rr, ratio, ends).calibrate(ins, n_grid=121)
        mo = HJM2FModel(md.curve, lam, c.sigma1, c.sigma2, rr)
        bermu[rr] = BermudanQuad2F(mo, bm, n_grid=61, n_quad=81).price()
        prem[rr] = (bermu[rr] - max(p for p, i in zip(c.model_prices, ins) if i.swaption.expiry in exf)) * 1e4
    check("premium at rho = -0.95 is about 46 bp", prem[-0.95], 46.05, 0.6)
    check("premium at rho = +0.9 is about 102 bp", prem[0.9], 101.9, 0.8)
    check("Bermudan moves by 12 % between rho = -0.95 and +0.9", (bermu[0.9] / bermu[-0.95] - 1) * 100, 11.6, 0.8)

    print(f"\n{sum(RESULTS)}/{len(RESULTS)} checks passed in {time.time() - t_start:.0f} s")
    return 0 if all(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
