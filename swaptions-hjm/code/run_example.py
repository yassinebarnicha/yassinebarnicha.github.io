from __future__ import annotations

# Import the calibration bootstrap.
from hjm1f.calibration_1f import CalibrationInstrument, CoterminalBootstrap1F

# Import the zero curve.
from hjm1f.curve import ZeroCurve

# Import the HJM model.
from hjm1f.hjm_1f import HJM1FModel

# Import the instruments.
from hjm1f.instruments import BermudanCoterminalSwaption, CoterminalSwaption

# Import the simple Monte-Carlo control pricer.
from hjm1f.mc_1f import BermudanLSMC1F, MonteCarloSpec

# Import the PDE engine and the PDE grid helper.
from hjm1f.pde_1f import BermudanPDE1F, PDEGrid


def main() -> None:
    # ----------------------------------------------------------------------
    # 1) INPUT ZERO CURVE
    # ----------------------------------------------------------------------
    curve = ZeroCurve(
        times=[0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0],
        discounts=[1.0, 0.978, 0.954, 0.930, 0.885, 0.842, 0.780],
    )

    # ----------------------------------------------------------------------
    # 2) BERMUDAN / COTERMINAL SETUP
    # ----------------------------------------------------------------------
    # Exercise dates of the Bermudan.
    exercise_times = [1.0, 2.0, 3.0]

    # Final maturity common to all coterminal underlyings.
    terminal = 7.0

    # Fixed payment dates of the full terminal swap.
    fixed_payment_times = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

    # Yearly accruals for simplicity.
    fixed_accruals = [1.0] * len(fixed_payment_times)

    # Fixed strike of the swap / swaption.
    strike = 0.03

    # Market quotes assumed in Bachelier normal vols.
    normal_vols = [0.0080, 0.0085, 0.0090]

    # ----------------------------------------------------------------------
    # 3) BUILD THE COTERMINAL CALIBRATION INSTRUMENTS
    # ----------------------------------------------------------------------
    instruments = []
    for ex, nv in zip(exercise_times, normal_vols):
        swpt = CoterminalSwaption(
            expiry=ex,
            terminal=terminal,
            strike=strike,
            fixed_payment_times=fixed_payment_times,
            fixed_accruals=fixed_accruals,
            payer=True,
        )
        instruments.append(CalibrationInstrument(swpt, nv))

    # ----------------------------------------------------------------------
    # 4) BOOTSTRAP PIECEWISE-CONSTANT SIGMA(t)
    # ----------------------------------------------------------------------
    bootstrap = CoterminalBootstrap1F(curve, lam=0.03, bucket_ends=exercise_times, initial_sigma=0.01)
    calib = bootstrap.calibrate(instruments)

    print("Calibrated sigmas:", calib.sigma_curve.values)
    print("Market prices    :", [round(x, 8) for x in calib.market_prices])
    print("Model prices     :", [round(x, 8) for x in calib.model_prices])

    # ----------------------------------------------------------------------
    # 5) BUILD THE FINAL CALIBRATED MODEL
    # ----------------------------------------------------------------------
    model = HJM1FModel(curve, lam=0.03, sigma_curve=calib.sigma_curve)

    # ----------------------------------------------------------------------
    # 6) BUILD THE BERMUDAN INSTRUMENT
    # ----------------------------------------------------------------------
    berm = BermudanCoterminalSwaption(
        exercise_times=exercise_times,
        terminal=terminal,
        strike=strike,
        fixed_payment_times=fixed_payment_times,
        fixed_accruals=fixed_accruals,
        payer=True,
    )

    # ----------------------------------------------------------------------
    # 7) PDE GRID
    # ----------------------------------------------------------------------
    # IMPORTANT:
    # The natural scale is sqrt(alpha(T_max)), because alpha(T_max) is the
    # variance of X_{T_max}. So a standard practical choice is +/- 5 standard
    # deviations, not +/- 5 times alpha.
    grid = PDEGrid.around_factor_std(model, t_max=exercise_times[-1], n_std=5.0, nx=301, nt=500)

    print("PDE x-range      :", (round(grid.x_min, 6), round(grid.x_max, 6)))
    print("sqrt(alpha(Tmax)):", round(model.alpha(exercise_times[-1]) ** 0.5, 6))

    # ----------------------------------------------------------------------
    # 8) PRICE THE BERMUDAN WITH THE PDE ENGINE
    # ----------------------------------------------------------------------
    pde = BermudanPDE1F(model, berm, grid)
    pde_price = pde.price()
    print("Bermudan PDE price :", round(pde_price, 8))

    # ----------------------------------------------------------------------
    # 9) OPTIONAL MONTE-CARLO CROSS-CHECK
    # ----------------------------------------------------------------------
    lsmc = BermudanLSMC1F(model, berm, MonteCarloSpec(n_paths=10000, seed=1234))
    mc_price = lsmc.price()
    print("Bermudan LSMC price:", round(mc_price, 8))


if __name__ == "__main__":
    main()
