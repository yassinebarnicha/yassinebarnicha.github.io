from __future__ import annotations

# Standard-library imports.
import math
from dataclasses import dataclass

# Third-party imports.
import numpy as np

# Local imports.
from .hjm_1f import HJM1FModel
from .instruments import BermudanCoterminalSwaption


@dataclass
class MonteCarloSpec:
    """Configuration of the Monte-Carlo control pricer."""

    n_paths: int = 20000
    seed: int = 1234


class BermudanLSMC1F:
    """
    Simple Longstaff-Schwartz control pricer.

    This is not the main engine of the project. It is included mainly as a
    rough cross-check against the PDE price.
    """

    def __init__(self, model: HJM1FModel, bermudan: BermudanCoterminalSwaption, spec: MonteCarloSpec):
        # Store the model.
        self.model = model

        # Store the Bermudan spec.
        self.bermudan = bermudan

        # Store the Monte-Carlo settings.
        self.spec = spec

    def price(self) -> float:
        # Bermudan exercise dates.
        ex_times = list(self.bermudan.exercise_times)

        # Random generator.
        rng = np.random.default_rng(self.spec.seed)

        # Number of exercise dates.
        n = len(ex_times)

        # Simulated factor paths X_{T_j}.
        x = np.zeros((self.spec.n_paths, n))

        # alpha(T_{j-1}) to build independent Gaussian increments.
        alpha_prev = 0.0

        # Simulate X at the exercise dates only.
        for j, t in enumerate(ex_times):
            # Variance Var[X_t].
            alpha_t = self.model.alpha(t)

            # Incremental variance because X_t is a Gaussian martingale.
            var = max(alpha_t - alpha_prev, 0.0)

            # Increment the paths.
            x[:, j] = (x[:, j - 1] if j > 0 else 0.0) + math.sqrt(var) * rng.standard_normal(self.spec.n_paths)

            # Update previous variance.
            alpha_prev = alpha_t

        # Realized cashflow on each path.
        cashflows = np.zeros(self.spec.n_paths)

        # Exercise time chosen on each path.
        exercise_time = np.full(self.spec.n_paths, ex_times[-1])

        # Last underlying swaption.
        last_swaption = self.bermudan.underlying_at(ex_times[-1])

        # Final-time exercise payoff.
        for i in range(self.spec.n_paths):
            cashflows[i] = max(self.model.swap_value(last_swaption, float(x[i, -1])), 0.0)

        # Backward Longstaff-Schwartz regression.
        for j in reversed(range(n - 1)):
            # Current exercise date.
            t = ex_times[j]

            # Underlying coterminal swaption at this date.
            und = self.bermudan.underlying_at(t)

            # Immediate exercise values on all paths.
            ex_payoff = np.array([max(self.model.swap_value(und, float(xx)), 0.0) for xx in x[:, j]])

            # Very simple discounting approximation from the current date to the future stopping time.
            df = np.array([
                math.exp(-self.model.short_rate(t, float(x[i, j])) * (exercise_time[i] - t)) for i in range(self.spec.n_paths)
            ])

            # Discounted continuation values.
            y = cashflows * df

            # In-the-money paths only.
            itm = ex_payoff > 1e-14

            if np.any(itm):
                # Basis regressor = [1, X, X^2].
                X = x[itm, j]
                A = np.column_stack([np.ones_like(X), X, X * X])

                # Least-squares fit of the continuation value.
                beta, *_ = np.linalg.lstsq(A, y[itm], rcond=None)

                # Predicted continuation values.
                continuation = beta[0] + beta[1] * X + beta[2] * X * X

                # Exercise decision.
                exercise = ex_payoff[itm] > continuation

                # Indices of paths where exercise is optimal.
                idx = np.where(itm)[0][exercise]

                # Update the pathwise cashflows and stopping times.
                cashflows[idx] = ex_payoff[idx]
                exercise_time[idx] = t

        # Discount all selected cashflows back to today.
        discounts = np.array([math.exp(-self.model.short_rate(0.0, 0.0) * exercise_time[i]) for i in range(self.spec.n_paths)])

        # Monte-Carlo estimate.
        return float(np.mean(cashflows * discounts))
