from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from hjm1f.instruments import BermudanCoterminalSwaption

from .hjm_2f import HJM2FModel, psd_sqrt


@dataclass
class MonteCarloSpec2F:
    n_paths: int = 100000
    seed: int = 1234


class BermudanLSMC2F:
    """
    Longstaff-Schwartz in two factors, run under the terminal forward measure with numeraire P(t, T_N).

    Deflated values V/P(t,T_N) are martingales, so no path-dependent discount factor is needed and
    X can be sampled exactly at the exercise dates (the drift of X is deterministic).
    """

    def __init__(self, model: HJM2FModel, bermudan: BermudanCoterminalSwaption, spec: MonteCarloSpec2F):
        self.model = model
        self.bermudan = bermudan
        self.spec = spec

    def _deflated_payoff(self, t: float, x: np.ndarray) -> np.ndarray:
        und = self.bermudan.underlying_at(t)
        c, d1, d2 = self.model.bond_coeffs(t, self.bermudan.terminal)
        pn = np.exp(c - d1 * x[:, 0] - d2 * x[:, 1])
        return np.maximum(self.model.swap_value_vec(und, x[:, 0], x[:, 1]), 0.0) / pn

    def price(self) -> float:
        ex = list(self.bermudan.exercise_times)
        n_ex = len(ex)
        n = self.spec.n_paths
        tn = self.bermudan.terminal
        rng = np.random.default_rng(self.spec.seed)

        x = np.zeros((n, n_ex, 2))
        prev_m, prev_a, cur = np.zeros(2), np.zeros((2, 2)), np.zeros((n, 2))
        for j, t in enumerate(ex):
            m, a = self.model.mean_under_forward(t, tn), self.model.cov(t)
            chol = psd_sqrt(a - prev_a)
            half = rng.standard_normal((n // 2, 2))
            z = np.vstack([half, -half])
            cur = cur + (m - prev_m) + z @ chol.T
            x[:, j, :] = cur
            prev_m, prev_a = m, a

        payoff = [self._deflated_payoff(t, x[:, j, :]) for j, t in enumerate(ex)]
        value = payoff[-1].copy()
        for j in reversed(range(n_ex - 1)):
            itm = payoff[j] > 1e-14
            if not np.any(itm):
                continue
            xs = x[itm, j, :]
            s = xs.std(axis=0) + 1e-12
            u, v = xs[:, 0] / s[0], xs[:, 1] / s[1]
            basis = np.column_stack([np.ones_like(u), u, v, u * u, u * v, v * v, u ** 3, u * u * v, u * v * v, v ** 3])
            beta, *_ = np.linalg.lstsq(basis, value[itm], rcond=None)
            exercise = payoff[j][itm] > basis @ beta
            idx = np.where(itm)[0][exercise]
            value[idx] = payoff[j][idx]
        return self.model.curve.discount(tn) * float(np.mean(value))
