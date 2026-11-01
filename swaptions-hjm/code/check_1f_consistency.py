"""
No-arbitrage regression checks of the 1F code.

1. A single-date European swaption must cost the same through the PDE (discounting with
   short_rate, exercising with swap_value) as through Jamshidian, for payers and receivers.
2. The Bermudan PDE and the Longstaff-Schwartz engine must agree.
3. With sigma_2 = 0 the two-factor model must reproduce the one-factor prices.
"""
from __future__ import annotations

from hjm1f.calibration_1f import CalibrationInstrument, CoterminalBootstrap1F
from hjm1f.curve import ZeroCurve
from hjm1f.hjm_1f import HJM1FModel, PiecewiseConstantSigma
from hjm1f.instruments import BermudanCoterminalSwaption, CoterminalSwaption
from hjm1f.mc_1f import BermudanLSMC1F, MonteCarloSpec
from hjm1f.pde_1f import BermudanPDE1F, PDEGrid
from hjm2f.hjm_2f import HJM2FModel


def main() -> None:
    curve = ZeroCurve([0, 1, 2, 3, 5, 7, 10], [1.0, 0.978, 0.954, 0.930, 0.885, 0.842, 0.780])
    ex, pay, acc = [1.0, 2.0, 3.0], [2.0, 3.0, 4.0, 5.0, 6.0, 7.0], [1.0] * 6
    inst = [CalibrationInstrument(CoterminalSwaption(e, 7.0, 0.03, pay, acc, True), v) for e, v in zip(ex, [0.008, 0.0085, 0.009])]
    calib = CoterminalBootstrap1F(curve, 0.03, ex).calibrate(inst)
    model = HJM1FModel(curve, 0.03, calib.sigma_curve)

    print("1) European swaption, PDE vs Jamshidian")
    for payer in (True, False):
        for e in ex:
            jam = model.european_swaption_price_jamshidian(CoterminalSwaption(e, 7.0, 0.03, pay, acc, payer))
            grid = PDEGrid.around_factor_std(model, e, 5.0, nx=601, nt=800)
            pde = BermudanPDE1F(model, BermudanCoterminalSwaption([e], 7.0, 0.03, pay, acc, payer), grid).price()
            print(f"   {'payer' if payer else 'receiver'} {e:.0f}y  Jamshidian {jam:.8f}  PDE {pde:.8f}  ratio {pde / jam:.4f}")

    print("2) Bermudan payer, PDE vs LSMC")
    berm = BermudanCoterminalSwaption(ex, 7.0, 0.03, pay, acc, True)
    pde = BermudanPDE1F(model, berm, PDEGrid.around_factor_std(model, 3.0, 5.0, nx=801, nt=1500)).price()
    lsmc = BermudanLSMC1F(model, berm, MonteCarloSpec(n_paths=200000, seed=1)).price()
    print(f"   PDE {pde:.8f}  LSMC {lsmc:.8f}  ratio {lsmc / pde:.4f}")

    print("3) 2F model with sigma_2 = 0 vs 1F Jamshidian")
    zero = PiecewiseConstantSigma(ex, [0.0, 0.0, 0.0])
    model2 = HJM2FModel(curve, (0.03, 0.5), calib.sigma_curve, zero, -0.5)
    for e in ex:
        sw = CoterminalSwaption(e, 7.0, 0.03, pay, acc, True)
        print(f"   {e:.0f}y  1F {model.european_swaption_price_jamshidian(sw):.8f}  2F {model2.european_swaption_price(sw):.8f}")


if __name__ == "__main__":
    main()
