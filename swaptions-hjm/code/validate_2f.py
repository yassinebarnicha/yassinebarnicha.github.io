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
from hjm2f.calibration_2f import CoterminalBootstrap2F, atm_swaption
from hjm2f.calibration_global import fit_one_factor, fit_two_factor, from_params
from hjm2f.hjm_2f import HJM2FModel
from hjm2f.mc_2f import BermudanLSMC2F, MonteCarloSpec2F
from hjm2f.quad_2f import BermudanQuad2F

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

    print("\n== 4. Real market data (theory 2F, section 4; report, table 1) ==")
    d = json.loads((Path(__file__).resolve().parent / "data" / "eur_2016-03-30.json").read_text(encoding="utf-8"))
    cv = ZeroCurve([c["t"] for c in d["curve"]], [c["df"] for c in d["curve"]])
    mx = d["swaptionNormalVolATM"]
    vol = lambda e, t: mx["vol"][mx["expiries"].index(float(e))][mx["tenors"].index(int(t))]
    P10, A10 = [float(i) for i in range(1, 11)], [1.0] * 10
    K = forward_swap_rate(cv, 0.0, 10.0, P10, A10)
    check("10-year par rate = 0.532 %", K * 100, 0.532, 5e-4)
    EX = [1, 2, 3, 4, 5]
    exf = [float(e) for e in EX]
    ins = [CalibrationInstrument(CoterminalSwaption(e, 10.0, K, P10, A10, True), vol(e, 10 - e)) for e in EX]
    bm = BermudanCoterminalSwaption(exf, 10.0, K, P10, A10, True)
    targets = [(atm_swaption(cv, float(e), float(t)), vol(e, t)) for e in EX for t in range(1, 11)]
    lam1f, rmse1, err1 = fit_one_factor(cv, exf, ins, targets)
    check("1F: best lambda = 0 (lower bound)", lam1f, 0.0, 1e-3)
    check("1F: RMSE on the 50 swaptions = 11.2 bp", rmse1, 11.2, 0.05)
    check("1F: largest gap = 35 bp (1y x 1y)", max(abs(x) for x in err1), 35.0, 0.5)
    fit = fit_two_factor(cv, exf, ins, targets, from_params((0.01, 1.0), -0.9, 1.5), n_grid=61, max_fev=250)
    check("2F: RMSE = 2.6 bp (claimed 2.58)", fit.rmse_bp, 2.58, 0.03)
    check("2F: lambda_2 = 0.92 (claimed 0.916)", fit.lam[1], 0.916, 0.01)
    check("2F: rho = -0.98", fit.rho, -0.98, 0.005)
    check("2F: sigma2/sigma1 = 2.65", fit.ratio, 2.65, 0.02)
    check("2F: largest gap = 6 bp (claimed 'about 6')", max(abs(x) for x in fit.errors_bp), 6.2, 0.3)
    # other starting points (claimed: same fit with different parameters, or stuck at a bound)
    alt = fit_two_factor(cv, exf, ins, targets, [-3.0, 1.5, -0.9, 0.5], n_grid=61, max_fev=300)
    check("2F: second start reaches the same fit, 2.62 bp", alt.rmse_bp, 2.62, 0.05)
    at_most("2F: ... with lambda_2 <= 1.05", alt.lam[1], 1.05)
    at_most("2F: ... and rho <= -0.98 (close to -1)", alt.rho, -0.98)
    check("2F: ... and sigma2/sigma1 = 3.0", alt.ratio, 2.99, 0.1)
    stuck = fit_two_factor(cv, exf, ins, targets, [0.0, -1.0, 0.5, -0.5], n_grid=61, max_fev=300)
    check("2F: a bad start stalls at 7.2 bp against lambda_2 = 3, rho = -1", stuck.rmse_bp, 7.19, 0.2)
    c1 = CoterminalBootstrap1F(cv, lam1f, exf).calibrate(ins)
    m1x = HJM1FModel(cv, lam1f, c1.sigma_curve)
    b1 = BermudanPDE1F(m1x, bm, PDEGrid.around_factor_std(m1x, 5.0, 5.0, nx=801, nt=1500), theta=0.5, rannacher=2).price()
    c2 = CoterminalBootstrap2F(cv, fit.lam, fit.rho, fit.ratio, exf).calibrate(ins)
    m2x = HJM2FModel(cv, fit.lam, c2.sigma1, c2.sigma2, fit.rho)
    b2 = BermudanQuad2F(m2x, bm, n_grid=81, n_quad=101).price()
    eur = max(c2.model_prices)
    check("Bermudan 1F = 0.0518 (claimed 0.051765)", b1, 0.051765, 2e-5)
    check("Bermudan 2F = 0.0490 (claimed 0.048973)", b2, 0.048973, 3e-5)
    check("exercise premium 1F = 75 bp", (b1 - eur) * 1e4, 75.2, 0.3)
    check("exercise premium 2F = 47 bp", (b2 - eur) * 1e4, 47.2, 0.4)
    check("Bermudan 2F vs 1F = -5.4 %", (b2 / b1 - 1) * 100, -5.4, 0.1)

    print("\n== 5. Exercise premium against rho, illustrative case (page and report) ==")
    prem = {}
    for rho in (-0.95, 0.9):
        r = CoterminalBootstrap2F(curve, (0.03, 0.5), rho, 0.6, ex).calibrate(inst)
        mo = HJM2FModel(curve, (0.03, 0.5), r.sigma1, r.sigma2, rho)
        prem[rho] = (BermudanQuad2F(mo, berm, n_grid=81, n_quad=101).price() - max(r.model_prices)) * 1e4
    check("premium at rho = -0.95 is about 19 bp", prem[-0.95], 19.3, 0.5)
    check("premium at rho = +0.9 is about 29 bp", prem[0.9], 28.7, 0.5)
    p1 = BermudanPDE1F(model, berm, PDEGrid.around_factor_std(model, 3.0, 5.0, nx=801, nt=1500), theta=0.5, rannacher=2).price()
    check("1F premium = 24 bp", (p1 - max(calib.model_prices)) * 1e4, 24.2, 0.3)

    print(f"\n{sum(RESULTS)}/{len(RESULTS)} checks passed in {time.time() - t_start:.0f} s")
    return 0 if all(RESULTS) else 1


if __name__ == "__main__":
    sys.exit(main())
