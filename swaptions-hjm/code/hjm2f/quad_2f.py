from __future__ import annotations

import numpy as np

from hjm1f.instruments import BermudanCoterminalSwaption

from .hjm_2f import HJM2FModel, normal_grid_2d, psd_sqrt


def _bilinear(ax1: np.ndarray, ax2: np.ndarray, f: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
    n1, n2 = len(ax1), len(ax2)
    u = np.clip((p1 - ax1[0]) / (ax1[1] - ax1[0]), 0.0, n1 - 1 - 1e-9)
    v = np.clip((p2 - ax2[0]) / (ax2[1] - ax2[0]), 0.0, n2 - 1 - 1e-9)
    i, k = u.astype(int), v.astype(int)
    fu, fv = u - i, v - k
    return (f[i, k] * (1 - fu) * (1 - fv) + f[i + 1, k] * fu * (1 - fv)
            + f[i, k + 1] * (1 - fu) * fv + f[i + 1, k + 1] * fu * fv)


class BermudanQuad2F:
    """
    Backward induction on a 2-D grid of (X^1, X^2) at the exercise dates -- the 2F counterpart of the 1F PDE.

    Because X is Gaussian with deterministic drift, the transition between two exercise dates is
    known exactly, so there is no time discretisation: only the grid (and quadrature) error remains.
    It is run under the T_N-forward measure on values deflated by P(t, T_N).
    The cost grows like (grid points)^2 per factor, which is why a 3rd factor would call for LSMC.
    """

    def __init__(self, model: HJM2FModel, bermudan: BermudanCoterminalSwaption,
                 n_grid: int = 61, n_std: float = 4.5, n_quad: int = 81):
        self.model = model
        self.bermudan = bermudan
        self.n_grid = n_grid
        self.n_std = n_std
        self.n_quad = n_quad

    def _deflated_payoff(self, t: float, x1: np.ndarray, x2: np.ndarray) -> np.ndarray:
        und = self.bermudan.underlying_at(t)
        c, d1, d2 = self.model.bond_coeffs(t, self.bermudan.terminal)
        pn = np.exp(c - d1 * x1 - d2 * x2)
        return np.maximum(self.model.swap_value_vec(und, x1, x2), 0.0) / pn

    def price(self) -> float:
        ex = list(self.bermudan.exercise_times)
        tn = self.bermudan.terminal
        means = [self.model.mean_under_forward(t, tn) for t in ex]
        covs = [self.model.cov(t) for t in ex]
        z, w = normal_grid_2d(self.n_quad)

        axes = []
        for m, a in zip(means, covs):
            sd = np.sqrt(np.maximum(np.diag(a), 1e-18))
            axes.append((np.linspace(m[0] - self.n_std * sd[0], m[0] + self.n_std * sd[0], self.n_grid),
                         np.linspace(m[1] - self.n_std * sd[1], m[1] + self.n_std * sd[1], self.n_grid)))

        g1, g2 = np.meshgrid(*axes[-1], indexing="ij")
        u = self._deflated_payoff(ex[-1], g1, g2)

        for j in reversed(range(len(ex) - 1)):
            a1, a2 = axes[j]
            b1, b2 = axes[j + 1]
            chol = psd_sqrt(covs[j + 1] - covs[j])
            shift = means[j + 1] - means[j]
            inc = shift + z @ chol.T
            g1, g2 = np.meshgrid(a1, a2, indexing="ij")
            cont = np.zeros_like(g1)
            for r in range(self.n_grid):
                p1 = g1[r][:, None] + inc[None, :, 0]
                p2 = g2[r][:, None] + inc[None, :, 1]
                cont[r] = _bilinear(b1, b2, u, p1, p2) @ w
            u = np.maximum(self._deflated_payoff(ex[j], g1, g2), cont)

        chol0 = psd_sqrt(covs[0])
        pts = means[0] + z @ chol0.T
        v0 = _bilinear(axes[0][0], axes[0][1], u, pts[:, 0], pts[:, 1]) @ w
        return self.model.curve.discount(tn) * float(v0)
