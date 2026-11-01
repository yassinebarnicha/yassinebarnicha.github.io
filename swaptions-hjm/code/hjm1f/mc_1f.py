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
    # degree of the regression polynomial (2: basis [1, X, X^2])
    degree: int = 2
    # regress on one set of paths, price on another: removes the in-sample (high) bias
    independent: bool = False


class BermudanLSMC1F:
    """
    Longstaff-Schwartz control pricer.

    This is not the main engine of the project: it is a cross-check of the PDE price.

    It runs under the terminal forward measure with numeraire P(t, T_N). Values deflated by this
    numeraire are martingales, so no path-dependent discount factor is needed, and because X is
    a Gaussian process with a deterministic drift it can be sampled exactly at the exercise dates.
    """

    def __init__(self, model: HJM1FModel, bermudan: BermudanCoterminalSwaption, spec: MonteCarloSpec):
        # Store the model.
        self.model = model

        # Store the Bermudan spec.
        self.bermudan = bermudan

        # Store the Monte-Carlo settings.
        self.spec = spec

    def _deflated_payoff(self, t: float, x: np.ndarray) -> np.ndarray:
        und = self.bermudan.underlying_at(t)
        terminal = self.bermudan.terminal
        out = np.empty(len(x))
        for i, xi in enumerate(x):
            out[i] = max(self.model.swap_value(und, float(xi)), 0.0) / self.model.discount_bond(t, terminal, float(xi))
        return out

    def _simulate(self, rng: np.random.Generator) -> np.ndarray:
        """X at the exercise dates under the T_N-forward measure (exact Gaussian increments)."""
        ex_times = list(self.bermudan.exercise_times)
        terminal = self.bermudan.terminal
        n_paths = self.spec.n_paths
        x = np.zeros((n_paths, len(ex_times)))
        prev_mean, prev_var, cur = 0.0, 0.0, np.zeros(n_paths)
        for j, t in enumerate(ex_times):
            mean = self.model.mean_under_forward(t, terminal)
            var = self.model.alpha(t)
            cur = cur + (mean - prev_mean) + math.sqrt(max(var - prev_var, 0.0)) * rng.standard_normal(n_paths)
            x[:, j] = cur
            prev_mean, prev_var = mean, var
        return x

    def price(self) -> float:
        ex_times = list(self.bermudan.exercise_times)
        n = len(ex_times)
        terminal = self.bermudan.terminal
        rng = np.random.default_rng(self.spec.seed)
        degree = max(1, int(self.spec.degree))

        x = self._simulate(rng)
        payoff = [self._deflated_payoff(t, x[:, j]) for j, t in enumerate(ex_times)]
        scale = [math.sqrt(max(self.model.alpha(t), 1e-300)) for t in ex_times]

        def basis(xs: np.ndarray, j: int) -> np.ndarray:
            z = xs / scale[j]
            return np.column_stack([z ** k for k in range(degree + 1)])

        # Deflated value of the strategy, starting from the last exercise date.
        value = payoff[-1].copy()
        betas = [None] * n

        # Backward Longstaff-Schwartz regression.
        for j in reversed(range(n - 1)):
            # In-the-money paths only.
            itm = payoff[j] > 1e-14
            if not np.any(itm):
                continue
            a = basis(x[itm, j], j)
            beta, *_ = np.linalg.lstsq(a, value[itm], rcond=None)
            betas[j] = beta
            exercise = payoff[j][itm] > a @ beta
            idx = np.where(itm)[0][exercise]
            value[idx] = payoff[j][idx]

        if self.spec.independent:
            # new paths: the exercise rule learnt above is applied forwards
            x2 = self._simulate(rng)
            payoff2 = [self._deflated_payoff(t, x2[:, j]) for j, t in enumerate(ex_times)]
            value = payoff2[-1].copy()
            done = np.zeros(self.spec.n_paths, dtype=bool)
            for j in range(n - 1):
                if betas[j] is None:
                    continue
                itm = (payoff2[j] > 1e-14) & ~done
                if not np.any(itm):
                    continue
                cont = basis(x2[itm, j], j) @ betas[j]
                go = payoff2[j][itm] > cont
                idx = np.where(itm)[0][go]
                value[idx] = payoff2[j][idx]
                done[idx] = True

        # Monte-Carlo estimate: V(0) = P(0,T_N) * E^{T_N}[deflated cashflow].
        return float(self.model.curve.discount(terminal) * np.mean(value))
