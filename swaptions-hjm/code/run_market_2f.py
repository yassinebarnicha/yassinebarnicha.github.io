"""
1F versus 2F on the whole real market of the course data (EUR, 30 March 2016): the stripped curve, the 21 x 14 matrix of ATM
swaption normal vols, the 16 ATM cap vols and the strike smiles (data/eur_2016-03-30.json).

Product: Bermudan payer swaption into a 10-year swap (annual fixed leg), exercisable every year from 1 year to LAST years, struck at
the 10-year par rate. The model has one volatility segment at EVERY expiry of the matrix up to LAST (1M, 2M, 3M, 6M, 9M, 1Y, 18M, 2Y,
...), each fitted exactly: on the co-terminal swaption quoted at the Bermudan's strike along the smile, or on an ATM swaption for the
other expiries. The structural parameters (mean reversions, correlation, volatility ratio) are then fitted on EVERY ATM swaption and
ATM cap that expires within the window, and the models are judged on the WHOLE market, including what lies beyond the window.

    python run_market_2f.py [LAST]       LAST = last exercise date in years, default 5 (about 10 minutes; 9 takes longer)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from hjm1f.calibration_1f import CoterminalBootstrap1F
from hjm1f.hjm_1f import HJM1FModel, PiecewiseConstantSigma
from hjm1f.instruments import BermudanCoterminalSwaption
from hjm1f.pde_1f import BermudanPDE1F, PDEGrid
from hjm2f.calibration_2f import CoterminalBootstrap2F
from hjm2f.calibration_global import fit_one_factor, fit_two_factor, from_params
from hjm2f.hjm_2f import HJM2FModel
from hjm2f.market_data import MarketData, summarise
from hjm2f.mc_2f import BermudanLSMC2F, MonteCarloSpec2F
from hjm2f.quad_2f import BermudanQuad2F

TERMINAL = 10
TENORS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]


def show(name, params, fit_rmse, s, berm, eur):
    print(f"{name:<19} {params:<46} fit {fit_rmse:5.2f} | fitted: swpt {s['swaptions_fitted']:5.2f} caps {s['caps_fitted']:5.2f}"
          f" | beyond: swpt {s['swaptions_beyond']:6.2f} caps {s['caps_beyond']:6.2f} | Bermudan {berm:.6f}, premium {(berm - eur) * 1e4:5.1f} bp", flush=True)


def main() -> None:
    last = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    exercise = list(range(1, last + 1))
    md = MarketData()
    strike, pay, acc = md.bermudan(TERMINAL, exercise)
    exf = [float(e) for e in exercise]
    window = float(last)
    ends, ins = md.bucket_plan(TERMINAL, exercise, use_smile=True)
    atm_quotes = md.coterminal_instruments(TERMINAL, exercise, use_smile=False)
    smile_quotes = md.coterminal_instruments(TERMINAL, exercise, use_smile=True)
    berm = BermudanCoterminalSwaption(exf, float(TERMINAL), strike, pay, acc, True)
    n_total = sum(1 for e in md.expiries for t in md.tenors if md.has(e, t))
    both = md.targets(window, TENORS, with_caps=True)
    sw_only = md.targets(window, TENORS, with_caps=False)
    print(f"Window {last} years: {len(ends)} volatility segments; fitted on {len(both.swaptions)} ATM swaptions and {len(both.caps)} ATM caps; "
          f"judged on all {n_total} swaptions and {len(md.cap_maturities)} caps.")
    print(f"Bermudan strike {strike * 1e4:.1f} bp; co-terminal quotes, ATM -> along the smile (bp): "
          + ", ".join(f"{a.normal_vol * 1e4:.1f}->{b.normal_vol * 1e4:.1f}" for a, b in zip(atm_quotes, smile_quotes)) + "\n", flush=True)

    out = {"last_exercise": last, "n_fitted_swaptions": len(both.swaptions), "n_fitted_caps": len(both.caps), "n_swaptions": n_total}
    zero = PiecewiseConstantSigma(ends, [0.0] * len(ends))
    for label, targets, with_caps in (("swaptions only", sw_only, False), ("swaptions + caps", both, True)):
        t0 = time.time()
        lam1, r1, _ = fit_one_factor(md.curve, ends, ins, targets)
        c1 = CoterminalBootstrap1F(md.curve, lam1, ends).calibrate(ins)
        m1 = HJM1FModel(md.curve, lam1, c1.sigma_curve)
        s1 = summarise(md.report_card(HJM2FModel(md.curve, (lam1, 0.5), c1.sigma_curve, zero, 0.0), window, TENORS, with_caps))
        b1 = BermudanPDE1F(m1, berm, PDEGrid.around_factor_std(m1, window, 5.0, nx=801, nt=1500), theta=0.5, rannacher=2).price()
        eur1 = max(p for p, i in zip(c1.model_prices, ins) if i.swaption.expiry in exf)
        show(f"1F, {label}", f"lambda={lam1:+.4f}", r1, s1, b1, eur1)
        out[f"1F {label}"] = {"lam": lam1, "fit_rmse": r1, "bermudan": b1, "summary": s1}

        best = None
        for x0 in (from_params((0.0, 1.0), -0.8, 1.5), from_params((-0.1, 0.5), -0.5, 1.0)):
            f = fit_two_factor(md.curve, ends, ins, targets, x0, n_grid=61, max_fev=300)
            print(f"   start -> lam=({f.lam[0]:+.3f},{f.lam[1]:+.3f}) rho={f.rho:+.3f} ratio={f.ratio:.2f} rmse {f.rmse_bp:.2f} [{time.time() - t0:.0f}s]", flush=True)
            if best is None or f.rmse_bp < best.rmse_bp:
                best = f
        f = best
        c2 = CoterminalBootstrap2F(md.curve, f.lam, f.rho, f.ratio, ends).calibrate(ins, n_grid=161)
        m2 = HJM2FModel(md.curve, f.lam, c2.sigma1, c2.sigma2, f.rho)
        s2 = summarise(md.report_card(m2, window, TENORS, with_caps))
        b2 = BermudanQuad2F(m2, berm, n_grid=81, n_quad=101).price()
        l2 = BermudanLSMC2F(m2, berm, MonteCarloSpec2F(200000, 11)).price()
        eur2 = max(p for p, i in zip(c2.model_prices, ins) if i.swaption.expiry in exf)
        show(f"2F, {label}", f"lam=({f.lam[0]:+.3f},{f.lam[1]:+.3f}) rho={f.rho:+.3f} ratio={f.ratio:.2f}", f.rmse_bp, s2, b2, eur2)
        print(f"{'':<19} LSMC {l2:.6f}   sigma1 (bp) {[round(v * 1e4, 1) for v in c2.sigma1.values]}   [{time.time() - t0:.0f}s]\n", flush=True)
        out[f"2F {label}"] = {"lam": list(f.lam), "rho": f.rho, "ratio": f.ratio, "fit_rmse": f.rmse_bp, "bermudan": b2, "lsmc": l2,
                              "sigma1": list(c2.sigma1.values), "summary": s2}
    Path(f"market_results_{last}y.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
