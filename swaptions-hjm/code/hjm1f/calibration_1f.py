from __future__ import annotations

# Standard-library imports.
from dataclasses import dataclass
from typing import List, Sequence

# Local imports.
from .curve import ZeroCurve
from .hjm_1f import HJM1FModel, PiecewiseConstantSigma
from .instruments import CoterminalSwaption
from .market_quotes import bachelier_swaption_price
from .solvers import brent


@dataclass
class CalibrationInstrument:
    """One market calibration quote: a coterminal swaption + its normal vol."""

    swaption: CoterminalSwaption
    normal_vol: float


@dataclass
class CalibrationResult:
    """Returned object after the bootstrap."""

    sigma_curve: PiecewiseConstantSigma
    model_prices: List[float]
    market_prices: List[float]


class CoterminalBootstrap1F:
    """
    Sequential bootstrap on coterminal swaptions.

    For a Bermudan with exercise dates T1,...,TN we typically choose one sigma
    bucket per exercise interval and bootstrap in increasing expiry order.
    """

    def __init__(self, curve: ZeroCurve, lam: float, bucket_ends: Sequence[float], initial_sigma: float = 0.01):
        # Store the curve.
        self.curve = curve

        # Store lambda / mean reversion.
        self.lam = lam

        # Store the sigma bucket ends.
        self.bucket_ends = list(bucket_ends)

        # Store the initial guess for each bucket.
        self.initial_sigma = initial_sigma

    def calibrate(self, instruments: Sequence[CalibrationInstrument], lower: float = 1e-6, upper: float = 0.20) -> CalibrationResult:
        # Start from a flat initial sigma guess.
        sigma_vals = [self.initial_sigma] * len(self.bucket_ends)

        # We store model and market prices for diagnostics.
        model_prices = []
        market_prices = []

        # Bootstrap sequentially in expiry order.
        for i, inst in enumerate(instruments):
            # Convert the market normal vol quote into a market price.
            market_price = bachelier_swaption_price(self.curve, inst.swaption, inst.normal_vol)

            # Define the one-dimensional root function for the current bucket.
            def objective(sig_i: float) -> float:
                # Freeze the already calibrated buckets.
                sigma_vals_local = list(sigma_vals)

                # Replace only the current bucket value.
                sigma_vals_local[i] = sig_i

                # Rebuild the sigma curve.
                sigma_curve = PiecewiseConstantSigma(self.bucket_ends, sigma_vals_local)

                # Build the corresponding 1F HJM model.
                model = HJM1FModel(self.curve, self.lam, sigma_curve)

                # Model price minus market price.
                return model.european_swaption_price_jamshidian(inst.swaption) - market_price

            # Solve for the current sigma bucket.
            sigma_vals[i] = brent(objective, lower, upper, tol=1e-10, max_iter=200)

            # Rebuild the calibrated curve up to the current step.
            sigma_curve = PiecewiseConstantSigma(self.bucket_ends, sigma_vals)
            model = HJM1FModel(self.curve, self.lam, sigma_curve)

            # Store diagnostics.
            model_prices.append(model.european_swaption_price_jamshidian(inst.swaption))
            market_prices.append(market_price)

        # Return the final calibrated object.
        return CalibrationResult(PiecewiseConstantSigma(self.bucket_ends, sigma_vals), model_prices, market_prices)
