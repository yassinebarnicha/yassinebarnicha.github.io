"""
Third model: the fitted 2F with a stochastic-volatility multiplier, on the whole real market (EUR, 30 March 2016).

    1F        flat normal vol in the strike, one factor
    2F        flat normal vol in the strike, two factors
    2F-SV     two factors + log-normal stochastic volatility (kappa, xi, rho_v) fitted on the swaption smiles

The structural parameters of the 2F (lambda_1, lambda_2, rho, sigma_2/sigma_1) are those of run_market_2f.py (5-year window,
swaptions + caps). The stochastic volatility is fitted on the smiles that expire within the window, judged on the others,
then the 11 volatility segments are re-bootstrapped under the SV dynamics on the same co-terminal quotes (along the smile).

    python run_sv.py          (about 10 minutes)
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from hjm1f.instruments import BermudanCoterminalSwaption
from hjm1f.hjm_1f import PiecewiseConstantSigma
from hjm2f.calibration_2f import CoterminalBootstrap2F, atm_swaption, implied_normal_vol
from hjm2f.caps import model_cap_vol
from hjm2f.hjm_2f import HJM2FModel
from hjm2f.market_data import MarketData, rmse
from hjm2f.quad_2f import BermudanQuad2F
from hjm3sv.model import SVParams, europeans
from hjm3sv.pricing import (ZERO, SmileEngine, bermudan_sv, caps_vols, fit_sv, flat_smile_errors, smile_points)

LAM, RHO, RATIO = (-0.003414, 0.533868), -0.95, 1.483256
WINDOW = 5.0
GAIN = 1.6   # over-relaxation of the proportional update (the segments are coupled: a later quote also depends on the earlier levels)
TENORS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
IN_SAMPLE = [(0.25, 2), (0.25, 10), (1.0, 1), (1.0, 5), (1.0, 10), (2.0, 2), (2.0, 10), (5.0, 5), (5.0, 10)]
OUT_SAMPLE = [(10.0, 10), (10.0, 20), (15.0, 15)]


def log(msg: str, t0: float) -> None:
    print(f"[{time.time() - t0:5.0f}s] {msg}", flush=True)


def model_vols(base: HJM2FModel, sv: SVParams, swaptions, t_star: float, n_paths: int, seed: int, n_grid: int = 121) -> np.ndarray:
    exact = np.array([base.european_swaption_price(s, n_grid=n_grid) for s in swaptions])
    price = exact + (europeans(base, sv, swaptions, t_star, n_paths, seed) - europeans(base, ZERO, swaptions, t_star, n_paths, seed))
    out = []
    for s, p in zip(swaptions, price):
        try:
            out.append(implied_normal_vol(base.curve, s, p) if p > 0 else float("nan"))
        except Exception:
            out.append(float("nan"))
    return np.array(out)


def rebootstrap(md: MarketData, sv: SVParams, ends, ins, sigma0, n_paths: int, t0: float, iterations: int = 14):
    """Level of sigma_1 in each segment such that the SV model reproduces every quote (proportional fixed-point iteration)."""
    s = list(sigma0)
    sws = [i.swaption for i in ins]
    target = np.array([i.normal_vol for i in ins])
    for it in range(iterations):
        base = HJM2FModel(md.curve, LAM, PiecewiseConstantSigma(ends, s), PiecewiseConstantSigma(ends, [RATIO * v for v in s]), RHO)
        v = model_vols(base, sv, sws, 10.0, n_paths, 17)
        err = (v - target) * 1e4
        log(f"  re-bootstrap iteration {it}: worst |error| {np.max(np.abs(err)):.3f} bp", t0)
        if np.max(np.abs(err)) < 0.05:
            break
        s = [x * (tg / m) ** GAIN for x, tg, m in zip(s, target, v)]
    return base, err


def main() -> None:
    t0 = time.time()
    md = MarketData()
    ends, ins = md.bucket_plan(10, [1, 2, 3, 4, 5], use_smile=True)
    exf = [1.0, 2.0, 3.0, 4.0, 5.0]
    strike, pay, acc = md.bermudan(10, [1, 2, 3, 4, 5])
    berm = BermudanCoterminalSwaption(exf, 10.0, strike, pay, acc, True)
    res2 = CoterminalBootstrap2F(md.curve, LAM, RHO, RATIO, ends).calibrate(ins, n_grid=161)
    base2 = HJM2FModel(md.curve, LAM, res2.sigma1, res2.sigma2, RHO)
    out = {"lam": LAM, "rho": RHO, "ratio": RATIO}

    # ---------------------------------------------------------------- 1. the smile
    pts_in, pts_out = smile_points(md, IN_SAMPLE), smile_points(md, OUT_SAMPLE)
    flat_in, flat_out = flat_smile_errors(pts_in), flat_smile_errors(pts_out)
    r_flat_in, r_flat_out = rmse([e for _, e in flat_in]), rmse([e for _, e in flat_out])
    log(f"Gaussian models (flat in the strike): smile RMSE {r_flat_in:.2f} bp in sample ({len(flat_in)} points), {r_flat_out:.2f} bp out of sample ({len(flat_out)})", t0)
    eng = SmileEngine(base2, pts_in, 10.0, n_paths=40000)
    best = None
    for start in (SVParams(0.3, 0.8, 0.5), SVParams(1.0, 1.5, 0.2)):
        sv, f, n = fit_sv(eng, start, max_fev=100)
        log(f"SV fit from {start}: kappa {sv.kappa:.3f} xi {sv.xi:.3f} rho_v {sv.rho_v:.3f}  smile RMSE {f:.3f} bp ({n} evaluations)", t0)
        if best is None or f < best[1]:
            best = (sv, f)
    sv = best[0]
    out["sv"] = {"kappa": sv.kappa, "xi": sv.xi, "rho_v": sv.rho_v, "smile_rmse_in": best[1], "flat_in": r_flat_in, "flat_out": r_flat_out}

    # ---------------------------------------------------------------- 2. segments under the SV dynamics
    base3, err = rebootstrap(md, sv, ends, ins, res2.sigma1.values, 60000, t0)
    out["sigma1_2F"], out["sigma1_SV"] = list(res2.sigma1.values), list(base3.sigma1.values)
    log("sigma1 (bp), 2F : " + " ".join(f"{v * 1e4:.1f}" for v in res2.sigma1.values), t0)
    log("sigma1 (bp), SV : " + " ".join(f"{v * 1e4:.1f}" for v in base3.sigma1.values), t0)

    # smile after the re-bootstrap, in and out of sample (shape only: the level is the segments' job)
    eng3_in = SmileEngine(base3, pts_in, 10.0, n_paths=60000)
    e_in = [e for _, e in eng3_in.errors_bp(sv)]
    eng3_out = SmileEngine(base3, pts_out, 30.0, n_paths=60000)
    e_out = [e for _, e in eng3_out.errors_bp(sv)]
    out["smile"] = {"flat_in": r_flat_in, "flat_out": r_flat_out, "sv_in": rmse(e_in), "sv_out": rmse(e_out), "n_in": len(e_in), "n_out": len(e_out)}
    log(f"smile RMSE: in sample flat {r_flat_in:.2f} -> SV {rmse(e_in):.2f} bp ; out of sample flat {r_flat_out:.2f} -> SV {rmse(e_out):.2f} bp", t0)

    # ---------------------------------------------------------------- 3. report card: the 240 ATM swaptions and 16 caps
    rows = [(e, t) for e in md.expiries for t in md.tenors if md.has(e, t)]
    sws = [atm_swaption(md.curve, e, float(t)) for e, t in rows]
    v2 = np.array([implied_normal_vol(md.curve, s, base2.european_swaption_price(s, n_grid=121)) for s in sws])
    v3 = model_vols(base3, sv, sws, 30.0, 40000, 31)
    mk = np.array([md.swaption_vol(e, t) for e, t in rows])
    group = np.array(["fitted" if (e <= WINDOW + 1e-9 and t in TENORS) else ("window" if e <= WINDOW + 1e-9 else "beyond") for e, t in rows])
    card = {}
    for name, v in (("2F", v2), ("2F-SV", v3)):
        err = (v - mk) * 1e4
        card[name] = {g: rmse([x for x, gg in zip(err, group) if gg == g and not np.isnan(x)]) for g in ("fitted", "beyond")}
        card[name]["n_nan"] = int(np.isnan(err).sum())
    cm = md.cap_maturities
    cap2 = np.array([model_cap_vol(base2, m) for m in cm])
    cap3 = np.array(caps_vols(base3, sv, cm, 40000, 41, 30.0))
    capm = np.array(md.cap_vols)
    cg = np.array(["fitted" if m - 0.5 <= WINDOW + 1e-9 else "beyond" for m in cm])
    for name, v in (("2F", cap2), ("2F-SV", cap3)):
        err = (v - capm) * 1e4
        card[name]["caps_fitted"] = rmse([x for x, g in zip(err, cg) if g == "fitted"])
        card[name]["caps_beyond"] = rmse([x for x, g in zip(err, cg) if g == "beyond"])
    out["card"] = card
    for name in ("2F", "2F-SV"):
        c = card[name]
        log(f"{name:6s} swaptions fitted {c['fitted']:.2f} beyond {c['beyond']:.2f} | caps fitted {c['caps_fitted']:.2f} beyond {c['caps_beyond']:.2f}", t0)

    # ---------------------------------------------------------------- 4. the Bermudan
    quad2 = BermudanQuad2F(base2, berm, n_grid=81, n_quad=101).price()
    quad3 = BermudanQuad2F(base3, berm, n_grid=81, n_quad=101).price()
    b3, l_sv, l_0 = bermudan_sv(base3, sv, berm, quad3, n_paths=150000)
    eur2 = max(res2.model_prices)
    out["bermudan"] = {"2F": quad2, "2F_premium_bp": (quad2 - eur2) * 1e4, "SV": b3, "SV_premium_bp": (b3 - eur2) * 1e4, "SV_lsmc": l_sv, "zero_lsmc": l_0, "base3_gaussian_quad": quad3,
                       "max_european": eur2}
    log(f"Bermudan: 2F {quad2:.6f} (premium {(quad2 - eur2) * 1e4:.1f} bp) | 2F-SV {b3:.6f} (premium {(b3 - eur2) * 1e4:.1f} bp); LSMC_SV {l_sv:.6f}, LSMC_0 {l_0:.6f}, quad on SV segments {quad3:.6f}", t0)
    Path("sv_results.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    log("done", t0)


if __name__ == "__main__":
    main()
