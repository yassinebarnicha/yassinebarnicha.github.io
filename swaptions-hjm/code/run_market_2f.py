"""
1F versus 2F on real market data: the EUR curve and the ATM normal-volatility swaption matrix of 30 March 2016
(data/eur_2016-03-30.json, from the ENSAE course material).

Product: Bermudan payer swaption into a 10-year swap (annual fixed leg), exercisable at 1, 2, 3, 4 and 5 years,
struck at the 10-year par rate. The Bermudan's co-terminal Europeans (1y x 9y, ..., 5y x 5y) are matched exactly;
the other 45 ATM swaptions (expiries 1-5 years, tenors 1-10 years) are what the extra factor has to explain.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from hjm1f.calibration_1f import CalibrationInstrument, CoterminalBootstrap1F
from hjm1f.curve import ZeroCurve
from hjm1f.hjm_1f import HJM1FModel
from hjm1f.instruments import BermudanCoterminalSwaption, CoterminalSwaption
from hjm1f.market_quotes import forward_swap_rate
from hjm1f.pde_1f import BermudanPDE1F, PDEGrid
from hjm2f.calibration_2f import atm_swaption
from hjm2f.calibration_global import fit_one_factor, fit_two_factor, from_params
from hjm2f.hjm_2f import HJM2FModel
from hjm2f.mc_2f import BermudanLSMC2F, MonteCarloSpec2F
from hjm2f.quad_2f import BermudanQuad2F
from hjm2f.calibration_2f import CoterminalBootstrap2F

DATA = Path(__file__).resolve().parent / "data" / "eur_2016-03-30.json"
EXPIRIES = [1, 2, 3, 4, 5]
TENORS = list(range(1, 11))


def main() -> None:
    d = json.loads(DATA.read_text(encoding="utf-8"))
    curve = ZeroCurve([c["t"] for c in d["curve"]], [c["df"] for c in d["curve"]])
    m = d["swaptionNormalVolATM"]
    vol = lambda e, t: m["vol"][m["expiries"].index(float(e))][m["tenors"].index(int(t))]

    pay = [float(i) for i in range(1, 11)]
    acc = [1.0] * 10
    strike = forward_swap_rate(curve, 0.0, 10.0, pay, acc)
    ex = [float(e) for e in EXPIRIES]
    instruments = [CalibrationInstrument(CoterminalSwaption(e, 10.0, strike, pay, acc, True), vol(e, 10 - e)) for e in EXPIRIES]
    berm = BermudanCoterminalSwaption(ex, 10.0, strike, pay, acc, True)
    targets = [(atm_swaption(curve, float(e), float(t)), vol(e, t)) for e in EXPIRIES for t in TENORS]

    print(f"EUR {d['asof']}: 10y par rate {strike * 100:.4f} %, co-terminal vols (bp) "
          f"{[round(i.normal_vol * 1e4, 1) for i in instruments]}\n")

    t0 = time.time()
    lam_1f, rmse_1f, err_1f = fit_one_factor(curve, ex, instruments, targets)
    print(f"1F: best mean reversion lambda = {lam_1f:.4f}   RMSE on the 50 swaptions = {rmse_1f:.2f} bp   [{time.time() - t0:.0f}s]")

    t0 = time.time()
    x0 = from_params((0.01, 1.0), -0.9, 1.5)
    fit = fit_two_factor(curve, ex, instruments, targets, x0, n_grid=61, max_fev=250)
    print(f"2F: lambda = ({fit.lam[0]:.4f}, {fit.lam[1]:.4f})  rho = {fit.rho:.4f}  sigma2/sigma1 = {fit.ratio:.3f}   "
          f"RMSE = {fit.rmse_bp:.2f} bp   [{fit.evaluations} evaluations, {time.time() - t0:.0f}s]\n")

    print("Model - market, in bp of normal vol (rows: expiry; columns: tenor)")
    print(f"{'':>8}" + "".join(f"{t:>7}y" for t in TENORS))
    for name, errs in (("1F", err_1f), ("2F", fit.errors_bp)):
        print(name)
        for i, e in enumerate(EXPIRIES):
            print(f"{e:>6}y " + "".join(f"{errs[i * len(TENORS) + j]:>8.1f}" for j in range(len(TENORS))))

    # Bermudan under each calibration
    calib1 = CoterminalBootstrap1F(curve, lam_1f, ex).calibrate(instruments)
    m1 = HJM1FModel(curve, lam_1f, calib1.sigma_curve)
    pde = BermudanPDE1F(m1, berm, PDEGrid.around_factor_std(m1, ex[-1], 5.0, nx=801, nt=1500), theta=0.5, rannacher=2).price()
    res2 = CoterminalBootstrap2F(curve, fit.lam, fit.rho, fit.ratio, ex).calibrate(instruments)
    m2 = HJM2FModel(curve, fit.lam, res2.sigma1, res2.sigma2, fit.rho)
    quad = BermudanQuad2F(m2, berm, n_grid=81, n_quad=101).price()
    lsmc = BermudanLSMC2F(m2, berm, MonteCarloSpec2F(200000, 11)).price()
    eur = max(i for i in res2.model_prices)
    print(f"\nBermudan 10y, exercise 1..5y, strike {strike * 100:.3f} %:")
    print(f"  1F (PDE, Crank-Nicolson)     {pde:.6f}   sigma(bp) {[round(v * 1e4, 1) for v in calib1.sigma_curve.values]}")
    print(f"  2F (2-D quadrature)          {quad:.6f}   sigma1(bp) {[round(v * 1e4, 1) for v in res2.sigma1.values]}")
    print(f"  2F (LSMC, 200k paths)        {lsmc:.6f}")
    print(f"  most expensive European      {eur:.6f}")
    print(f"  switch premium: 1F {(pde - eur) * 1e4:.1f} bp, 2F {(quad - eur) * 1e4:.1f} bp of notional")


if __name__ == "__main__":
    main()
