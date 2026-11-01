from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np

from hjm1f.curve import ZeroCurve
from hjm1f.hjm_1f import PiecewiseConstantSigma
from hjm1f.instruments import CoterminalSwaption

_GL_X, _GL_W = np.polynomial.legendre.leggauss(24)


def _B(lam: float, t):
    """B_lambda(t) = (1 - exp(-lambda t)) / lambda, with the lambda -> 0 limit t."""
    if abs(lam) < 1e-12:
        return t
    return -np.expm1(-lam * t) / lam


def psd_sqrt(m: np.ndarray) -> np.ndarray:
    """L with L @ L.T = m for a symmetric positive semi-definite m (rank-deficient allowed)."""
    w, v = np.linalg.eigh(m)
    return v * np.sqrt(np.clip(w, 0.0, None))


def normal_grid_2d(n: int, half_width: float = 6.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Uniform n x n grid on [-half_width, half_width]^2 with normalised bivariate-normal weights.

    Payoffs like max(swap, 0) have a kink: on a kink a fine uniform rule (error ~ h^2) is more
    accurate than Gauss-Hermite, whose nodes are far apart near the centre.
    """
    z = np.linspace(-half_width, half_width, n)
    w = np.exp(-0.5 * z * z)
    w = w / w.sum()
    z1, z2 = np.meshgrid(z, z, indexing="ij")
    return np.column_stack([z1.ravel(), z2.ravel()]), np.outer(w, w).ravel()


@dataclass
class HJM2FModel:
    """
    Two-factor Gaussian HJM with exponentially decaying volatilities

        sigma_k(t, T) = sigma_k(t) exp(-lambda_k (T - t)),   k = 1, 2,   corr(dW1, dW2) = rho.

    State variables X^k_t = int_0^t sigma_k(u) exp(lambda_k u) dW^k_u are driftless Gaussian, and

        P(t,T) = P(0,T)/P(0,t) * exp( - sum_k D_k(t,T) X^k_t + det(t,T) ).

    With sigma_2 = 0 this is exactly the one-factor model, with the convexity term of the
    no-arbitrage drift (checked by `martingale_defect`).
    """

    curve: ZeroCurve
    lam: Tuple[float, float]
    sigma1: PiecewiseConstantSigma
    sigma2: PiecewiseConstantSigma
    rho: float
    _cache: Dict = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self) -> None:
        if list(self.sigma1.bucket_ends) != list(self.sigma2.bucket_ends):
            raise ValueError("both factors must share the same volatility buckets")
        if not -1.0 <= self.rho <= 1.0:
            raise ValueError("rho must lie in [-1, 1]")

    def _rho(self, k: int, l: int) -> float:
        return 1.0 if k == l else self.rho

    def _segments(self, t: float):
        ends = list(self.sigma1.bucket_ends)
        left = 0.0
        for i, end in enumerate(ends):
            right = t if (i == len(ends) - 1 or end >= t) else end
            if right > left:
                yield left, right, (self.sigma1.values[i], self.sigma2.values[i])
            if right >= t:
                return
            left = right

    def _integrate(self, t: float, fn) -> float:
        total = 0.0
        for a, b, s in self._segments(t):
            half = 0.5 * (b - a)
            u = 0.5 * (a + b) + half * _GL_X
            total += half * float(np.sum(_GL_W * fn(u, s)))
        return total

    def _AG(self, t: float) -> Tuple[np.ndarray, np.ndarray]:
        """A_kl(t) = Cov[X^k_t, X^l_t] and G_kl(t) = int rho_kl s_k s_l e^{(lk+ll)u} B_l(u) du."""
        key = ("AG", t)
        if key not in self._cache:
            a = np.zeros((2, 2))
            g = np.zeros((2, 2))
            for k in range(2):
                for l in range(2):
                    rl, lk, ll = self._rho(k, l), self.lam[k], self.lam[l]
                    a[k, l] = rl * self._integrate(t, lambda u, s: s[k] * s[l] * np.exp((lk + ll) * u))
                    g[k, l] = rl * self._integrate(t, lambda u, s: s[k] * s[l] * np.exp((lk + ll) * u) * _B(ll, u))
            self._cache[key] = (a, g)
        return self._cache[key]

    def cov(self, t: float) -> np.ndarray:
        return self._AG(t)[0]

    def bond_coeffs(self, t: float, T: float) -> Tuple[float, float, float]:
        """(const, D1, D2) with P(t,T;x) = exp(const - D1 x1 - D2 x2)."""
        key = ("bc", t, T)
        if key not in self._cache:
            a, g = self._AG(t)
            bt = [float(_B(self.lam[k], t)) for k in range(2)]
            bT = [float(_B(self.lam[k], T)) for k in range(2)]
            d = [bT[k] - bt[k] for k in range(2)]
            det = 0.0
            for k in range(2):
                for l in range(2):
                    det += (bT[k] * bT[l] - bt[k] * bt[l]) * a[k, l] - d[k] * g[k, l] - d[l] * g[l, k]
            det *= -0.5
            const = math.log(self.curve.discount(T) / self.curve.discount(t)) + det
            self._cache[key] = (const, d[0], d[1])
        return self._cache[key]

    def discount_bond(self, t: float, T: float, x) -> float:
        if T <= t:
            return 1.0
        c, d1, d2 = self.bond_coeffs(t, T)
        return math.exp(c - d1 * x[0] - d2 * x[1])

    def short_rate(self, t: float, x) -> float:
        a, g = self._AG(t)
        b = [math.exp(-self.lam[k] * t) for k in range(2)]
        bb = [float(_B(self.lam[k], t)) for k in range(2)]
        conv = 0.0
        for k in range(2):
            for l in range(2):
                conv += (b[k] * bb[l] + bb[k] * b[l]) * a[k, l] - b[k] * g[k, l] - b[l] * g[l, k]
        return self.curve.inst_fwd(t) + b[0] * x[0] + b[1] * x[1] + 0.5 * conv

    def mean_under_forward(self, t: float, T: float) -> np.ndarray:
        """E[X_t] under the T-forward measure (the drift is deterministic in a Gaussian model)."""
        key = ("mean", t, T)
        if key not in self._cache:
            m = np.zeros(2)
            for k in range(2):
                for l in range(2):
                    lk, ll = self.lam[k], self.lam[l]
                    m[k] -= self._rho(k, l) * self._integrate(
                        t, lambda u, s: s[k] * s[l] * np.exp(lk * u) * _B(ll, T - u)
                    )
            self._cache[key] = m
        return self._cache[key]

    def swap_value_vec(self, swaption: CoterminalSwaption, x1, x2):
        """Value at expiry of the underlying coterminal swap, vectorised over states."""
        times = swaption.remaining_payment_times()
        accruals = swaption.remaining_accruals()
        c, d1, d2 = self.bond_coeffs(swaption.expiry, swaption.terminal)
        pn = np.exp(c - d1 * x1 - d2 * x2)
        ann = 0.0
        for t_i, a_i in zip(times, accruals):
            ci, e1, e2 = self.bond_coeffs(swaption.expiry, t_i)
            ann = ann + a_i * np.exp(ci - e1 * x1 - e2 * x2)
        payer = 1.0 - pn - swaption.strike * ann
        return payer if swaption.payer else -payer

    def european_swaption_price(self, swaption: CoterminalSwaption, n_grid: int = 241) -> float:
        """
        Price = P(0,T) E^T[(swap value)^+] by 2-D quadrature of the Gaussian X_T.

        This replaces Jamshidian: with two factors the swap value is no longer monotone in a
        single scalar, so there is no single critical x* and no decomposition into bond options.
        """
        T = swaption.expiry
        if not swaption.remaining_payment_times():
            return 0.0
        m = self.mean_under_forward(T, T)
        chol = psd_sqrt(self.cov(T))
        z, w = normal_grid_2d(n_grid)
        x = m + z @ chol.T
        payoff = np.maximum(self.swap_value_vec(swaption, x[:, 0], x[:, 1]), 0.0)
        return self.curve.discount(T) * float(np.sum(w * payoff))

    def martingale_defect(self, t: float, b: float, T: float) -> float:
        """E^T[P(t,b)/P(t,T)] - P(0,b)/P(0,T): must be 0 if the bond formula and drift are consistent."""
        cb, d1b, d2b = self.bond_coeffs(t, b)
        cT, d1T, d2T = self.bond_coeffs(t, T)
        dd = np.array([d1b - d1T, d2b - d2T])
        m = self.mean_under_forward(t, T)
        a = self.cov(t)
        expo = (cb - cT) - dd @ m + 0.5 * dd @ a @ dd
        return math.exp(expo) - self.curve.discount(b) / self.curve.discount(T)
