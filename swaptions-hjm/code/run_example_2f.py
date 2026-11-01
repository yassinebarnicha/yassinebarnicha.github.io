from __future__ import annotations

import time

from hjm1f.calibration_1f import CalibrationInstrument
from hjm1f.curve import ZeroCurve
from hjm1f.instruments import BermudanCoterminalSwaption, CoterminalSwaption
from hjm2f.calibration_2f import CoterminalBootstrap2F, atm_swaption, implied_normal_vol
from hjm2f.hjm_2f import HJM2FModel
from hjm2f.mc_2f import BermudanLSMC2F, MonteCarloSpec2F
from hjm2f.quad_2f import BermudanQuad2F

LAM1 = 0.03
# (label, lambda_2, rho, sigma_2/sigma_1). ratio ~ 0 is the one-factor limit.
SCENARIOS = [
    ("1F limit", 0.5, -0.5, 1e-3),
    ("2F rho=-0.9", 0.5, -0.9, 0.6),
    ("2F rho=-0.5", 0.5, -0.5, 0.6),
    ("2F rho= 0.0", 0.5, 0.0, 0.6),
    ("2F rho=+0.5", 0.5, 0.5, 0.6),
    ("2F lam2=1, rho=-0.7", 1.0, -0.7, 1.0),
]


def main() -> None:
    curve = ZeroCurve([0, 1, 2, 3, 5, 7, 10], [1.0, 0.978, 0.954, 0.930, 0.885, 0.842, 0.780])
    exercise_times = [1.0, 2.0, 3.0]
    payments = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    accruals = [1.0] * len(payments)
    strike = 0.03
    normal_vols = [0.0080, 0.0085, 0.0090]

    instruments = [
        CalibrationInstrument(CoterminalSwaption(e, 7.0, strike, payments, accruals, True), v)
        for e, v in zip(exercise_times, normal_vols)
    ]
    bermudan = BermudanCoterminalSwaption(exercise_times, 7.0, strike, payments, accruals, True)

    print("Every scenario is calibrated EXACTLY to the same 3 coterminal swaption quotes.\n")
    header = f"{'scenario':<22}{'sigma1 (1y,2y,3y)':<28}{'fit err':>9}{'LSMC':>12}{'quad 2D':>12}{'max Europ.':>12}{'switch (bp)':>12}"
    print(header)
    models = {}
    for label, lam2, rho, ratio in SCENARIOS:
        t0 = time.time()
        res = CoterminalBootstrap2F(curve, (LAM1, lam2), rho, ratio, exercise_times).calibrate(instruments)
        model = HJM2FModel(curve, (LAM1, lam2), res.sigma1, res.sigma2, rho)
        models[label] = model
        fit_err = max(abs(a - b) for a, b in zip(res.model_prices, res.market_prices))
        lsmc = BermudanLSMC2F(model, bermudan, MonteCarloSpec2F(200000, 7)).price()
        quad = BermudanQuad2F(model, bermudan, n_grid=81, n_quad=101).price()
        max_eur = max(res.model_prices)
        sig = ",".join(f"{v * 1e4:.1f}" for v in res.sigma1.values) + " bp"
        print(f"{label:<22}{sig:<28}{fit_err:>9.0e}{lsmc:>12.6f}{quad:>12.6f}{max_eur:>12.6f}{(quad - max_eur) * 1e4:>12.2f}   [{time.time() - t0:.0f}s]")

    print("\nImplied normal vols (bp) of ATM European swaptions NOT used in calibration:")
    grid = [(1, 1), (1, 4), (1, 7), (3, 1), (3, 4), (3, 7)]
    print(f"{'scenario':<22}" + "".join(f"{f'{e}y x {t}y':>10}" for e, t in grid))
    for label, model in models.items():
        row = []
        for e, t in grid:
            sw = atm_swaption(curve, float(e), float(t))
            row.append(implied_normal_vol(curve, sw, model.european_swaption_price(sw)) * 1e4)
        print(f"{label:<22}" + "".join(f"{v:>10.1f}" for v in row))


if __name__ == "__main__":
    main()
