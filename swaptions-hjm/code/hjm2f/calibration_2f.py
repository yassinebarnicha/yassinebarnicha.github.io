from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

from hjm1f.calibration_1f import CalibrationInstrument
from hjm1f.curve import ZeroCurve
from hjm1f.hjm_1f import PiecewiseConstantSigma
from hjm1f.instruments import CoterminalSwaption
from hjm1f.market_quotes import bachelier_swaption_price, forward_swap_rate
from hjm1f.solvers import brent

from .hjm_2f import HJM2FModel


@dataclass
class CalibrationResult2F:
    sigma1: PiecewiseConstantSigma
    sigma2: PiecewiseConstantSigma
    model_prices: List[float]
    market_prices: List[float]


class CoterminalBootstrap2F:
    """
    Sequential bootstrap on coterminal swaptions, as in 1F, with sigma_2(t) = ratio * sigma_1(t).

    The structural parameters (lambda_1, lambda_2, rho, ratio) are inputs: a strip of coterminal
    swaptions only pins down the *level* of volatility at each exercise date, not the correlation
    structure. Identifying them needs swaptions on other tenors (see `run_example_2f.py`).
    """

    def __init__(self, curve: ZeroCurve, lam: Tuple[float, float], rho: float, ratio: float,
                 bucket_ends: Sequence[float], initial_sigma: float = 0.01):
        self.curve = curve
        self.lam = lam
        self.rho = rho
        self.ratio = ratio
        self.bucket_ends = list(bucket_ends)
        self.initial_sigma = initial_sigma

    def _model(self, sigma_vals: Sequence[float]) -> HJM2FModel:
        s1 = PiecewiseConstantSigma(self.bucket_ends, list(sigma_vals))
        s2 = PiecewiseConstantSigma(self.bucket_ends, [self.ratio * v for v in sigma_vals])
        return HJM2FModel(self.curve, self.lam, s1, s2, self.rho)

    def calibrate(self, instruments: Sequence[CalibrationInstrument], lower: float = 1e-6,
                  upper: float = 0.20, n_grid: int = 241) -> CalibrationResult2F:
        sigma_vals = [self.initial_sigma] * len(self.bucket_ends)
        model_prices: List[float] = []
        market_prices: List[float] = []
        for i, inst in enumerate(instruments):
            market_price = bachelier_swaption_price(self.curve, inst.swaption, inst.normal_vol)

            def objective(sig_i: float) -> float:
                trial = list(sigma_vals)
                trial[i] = sig_i
                return self._model(trial).european_swaption_price(inst.swaption, n_grid=n_grid) - market_price

            sigma_vals[i] = brent(objective, lower, upper, tol=1e-10, max_iter=200)
            model = self._model(sigma_vals)
            model_prices.append(model.european_swaption_price(inst.swaption, n_grid=n_grid))
            market_prices.append(market_price)
        model = self._model(sigma_vals)
        return CalibrationResult2F(model.sigma1, model.sigma2, model_prices, market_prices)


def implied_normal_vol(curve: ZeroCurve, swaption: CoterminalSwaption, price: float) -> float:
    """Bachelier normal vol that reproduces `price`."""
    return brent(lambda v: bachelier_swaption_price(curve, swaption, v) - price, 1e-6, 0.05, tol=1e-12, max_iter=200)


def atm_swaption(curve: ZeroCurve, expiry: float, tenor: float) -> CoterminalSwaption:
    """Annual-pay payer swaption struck at the forward swap rate (expiry x tenor)."""
    pay = [expiry + i for i in range(1, int(round(tenor)) + 1)]
    acc = [1.0] * len(pay)
    strike = forward_swap_rate(curve, expiry, pay[-1], pay, acc)
    return CoterminalSwaption(expiry, pay[-1], strike, pay, acc, True)
