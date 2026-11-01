"""
Two-factor Gaussian HJM with a stochastic volatility multiplier (third model, next to hjm1f and hjm2f).

Both forward-rate volatilities are scaled by sqrt(V_t),

    sigma_k(t, T) = sigma_k(t) m_t exp(-lambda_k (T - t)),      m_0 = 1,

with the multiplier written V_t = m_t^2 and m_t = exp(y_t - v(t)), a log-normal (SABR-like) Ornstein-Uhlenbeck volatility:
dy = -kappa y dt + xi dZ,  v(t) = xi^2 (1 - e^{-2 kappa t}) / (2 kappa)  (so that E[V_t] = 1 under the spot measure), and
dZ = rho_v dW^1 + sqrt(1 - rho_v^2) dB (B independent of W^1, W^2). The no-arbitrage drift of the forward rates holds path by
path, so the bond formula of the 2F model survives with RANDOM covariances A_kl(t), G_kl(t) (integrals of V against the same
deterministic weights):

    P(t,T) = P(0,T)/P(0,t) exp( - sum_k D_k(t,T) X^k_t + kappa(t,T; A, G) ).

Everything is simulated under the T*-forward measure, where X has the drift  -V_t sigma_k e^{lambda_k t} sum_l rho_kl sigma_l B_l(T*-t)
and y acquires the extra drift  -xi rho_v m_t (sigma_1 B_1(T*-t) + rho sigma_2 B_2(T*-t)).  Prices are deflated by P(t,T*).

With xi = 0 the simulation reproduces the exact Gaussian transition of the 2F model (the increments of A, G and of the mean are taken from the
2F model), so a run with xi = 0 on the same random numbers is a control variate for the exact 2F integral.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np

from hjm2f.hjm_2f import HJM2FModel, _B


@dataclass(frozen=True)
class SVParams:
    kappa: float = 1.0
    xi: float = 0.0
    rho_v: float = 0.0


@dataclass
class State:
    t: float
    x: np.ndarray       # (n, 2)
    v: np.ndarray       # (n,)
    a: np.ndarray       # (n, 3): A11, A12, A22
    g: np.ndarray       # (n, 4): G11, G12, G21, G22


def bond(base: HJM2FModel, t: float, T: float, st: State) -> np.ndarray:
    """P(t,T) on every path."""
    lam = base.lam
    bt = [float(_B(lam[k], t)) for k in range(2)]
    bT = [float(_B(lam[k], T)) for k in range(2)]
    d = [bT[k] - bt[k] for k in range(2)]
    a11, a12, a22 = st.a[:, 0], st.a[:, 1], st.a[:, 2]
    g11, g12, g21, g22 = st.g[:, 0], st.g[:, 1], st.g[:, 2], st.g[:, 3]
    det = ((bT[0] ** 2 - bt[0] ** 2) * a11 + 2.0 * (bT[0] * bT[1] - bt[0] * bt[1]) * a12 + (bT[1] ** 2 - bt[1] ** 2) * a22
           - 2.0 * (d[0] * g11 + d[0] * g12 + d[1] * g21 + d[1] * g22))
    det = -0.5 * det
    return base.curve.discount(T) / base.curve.discount(t) * np.exp(-d[0] * st.x[:, 0] - d[1] * st.x[:, 1] + det)


def swap_value(base: HJM2FModel, sw, st: State, p_terminal: np.ndarray) -> np.ndarray:
    ann = 0.0
    for t_i, a_i in zip(sw.remaining_payment_times(), sw.remaining_accruals()):
        ann = ann + a_i * bond(base, st.t, t_i, st)
    payer = 1.0 - p_terminal - sw.strike * ann
    return payer if sw.payer else -payer


def _segment_values(base: HJM2FModel, t: float) -> Tuple[float, float]:
    ends = list(base.sigma1.bucket_ends)
    for i, e in enumerate(ends):
        if t < e or i == len(ends) - 1:
            return base.sigma1.values[i], base.sigma2.values[i]
    raise AssertionError


def make_grid(obs_times: Sequence[float], fine_until: float = 5.0, dt_fine: float = 1.0 / 24.0, dt_coarse: float = 1.0 / 12.0) -> List[float]:
    last = max(obs_times)
    pts = set(round(float(t), 10) for t in obs_times)
    t = 0.0
    while t < last - 1e-12:
        t += dt_fine if t < fine_until - 1e-12 else dt_coarse
        pts.add(round(min(t, last), 10))
    pts.add(0.0)
    return sorted(pts)


def simulate(base: HJM2FModel, sv: SVParams, t_star: float, obs_times: Sequence[float], n_paths: int, seed: int,
             observer: Callable[[State], None]) -> None:
    """Simulate (X, V, A, G) under the T*-forward measure; call `observer` at each observation time."""
    grid = make_grid(obs_times)
    obs = set(round(float(t), 10) for t in obs_times)
    rng = np.random.default_rng(seed)
    n = n_paths - n_paths % 2
    x = np.zeros((n, 2))
    y = np.zeros(n)
    v = np.ones(n)
    a = np.zeros((n, 3))
    g = np.zeros((n, 4))
    rho = base.rho
    s_v = math.sqrt(max(1.0 - sv.rho_v ** 2, 0.0))
    a_prev, g_prev, m_prev = np.zeros((2, 2)), np.zeros((2, 2)), np.zeros(2)
    for i in range(1, len(grid)):
        t0, t1 = grid[i - 1], grid[i]
        dt = t1 - t0
        a1, g1 = base._AG(t1)
        m1 = base.mean_under_forward(t1, t_star)
        da, dg, dm = a1 - a_prev, g1 - g_prev, m1 - m_prev
        a_prev, g_prev, m_prev = a1, g1, m1
        l11 = math.sqrt(max(da[0, 0], 0.0))
        l21 = da[0, 1] / l11 if l11 > 1e-300 else 0.0
        l22 = math.sqrt(max(da[1, 1] - l21 * l21, 0.0))
        half = rng.standard_normal((n // 2, 3))
        u = np.vstack([half, -half])
        sqv = np.sqrt(v)
        x_new = x + dm[None, :] * v[:, None] + sqv[:, None] * np.column_stack([l11 * u[:, 0], l21 * u[:, 0] + l22 * u[:, 1]])
        y_new, v_new = y, v
        if sv.xi != 0.0:
            tm = 0.5 * (t0 + t1)
            s1, s2 = _segment_values(base, tm)
            c = s1 * float(_B(base.lam[0], t_star - tm)) + rho * s2 * float(_B(base.lam[1], t_star - tm))
            zv = sv.rho_v * u[:, 0] + s_v * u[:, 2]
            decay = math.exp(-sv.kappa * dt)
            sd = sv.xi * math.sqrt((1.0 - math.exp(-2.0 * sv.kappa * dt)) / (2.0 * sv.kappa))
            # clipped at +-4 so that the (rare) explosion of the forward-measure drift for rho_v < 0 cannot produce overflows
            y_new = np.clip(y * decay + sd * zv - sv.xi * sv.rho_v * sqv * c * dt, -4.0, 4.0)
            vt = sv.xi ** 2 * (1.0 - math.exp(-2.0 * sv.kappa * t1)) / (2.0 * sv.kappa)
            v_new = np.exp(2.0 * (y_new - vt))
        vbar = 0.5 * (v + v_new)
        a = a + vbar[:, None] * np.array([da[0, 0], da[0, 1], da[1, 1]])[None, :]
        g = g + vbar[:, None] * np.array([dg[0, 0], dg[0, 1], dg[1, 0], dg[1, 1]])[None, :]
        x, y, v = x_new, y_new, v_new
        if t1 in obs:
            observer(State(t1, x, v, a, g))


def europeans(base: HJM2FModel, sv: SVParams, swaptions: Sequence, t_star: float, n_paths: int, seed: int) -> np.ndarray:
    """Monte Carlo prices of European swaptions (deflated payoff, T* measure), shape (len,). One simulation for all."""
    out = np.zeros(len(swaptions))
    by_time: Dict[float, List[int]] = {}
    for i, s in enumerate(swaptions):
        by_time.setdefault(round(float(s.expiry), 10), []).append(i)

    def observer(st: State) -> None:
        idx = by_time.get(round(st.t, 10))
        if not idx:
            return
        pn = bond(base, st.t, t_star, st)
        for i in idx:
            out[i] = _price_one(base, swaptions[i], st, pn, t_star)

    simulate(base, sv, t_star, list(by_time.keys()), n_paths, seed, observer)
    return out


def _price_one(base: HJM2FModel, sw, st: State, pn: np.ndarray, t_star: float) -> float:
    p_term = bond(base, st.t, sw.terminal, st)
    return base.curve.discount(t_star) * float(np.mean(np.maximum(swap_value(base, sw, st, p_term), 0.0) / pn))


def cap_puts(base: HJM2FModel, sv: SVParams, fixings: Sequence[float], delta: float, strike_bond: float, t_star: float, n_paths: int, seed: int) -> np.ndarray:
    """MC value of the puts max(k - P(t, t+delta), 0) paid at t+delta, discounted to 0, for each fixing t."""
    out = np.zeros(len(fixings))
    where = {round(float(t), 10): i for i, t in enumerate(fixings)}

    def observer(st: State) -> None:
        i = where.get(round(st.t, 10))
        if i is None:
            return
        pn = bond(base, st.t, t_star, st)
        p = bond(base, st.t, st.t + delta, st)
        # put option on P(t, t+delta) paid at t: price = P(0,t) E^t[...]; here deflated by P(t, T*): payoff * P(t,t) = payoff at t
        out[i] = base.curve.discount(t_star) * float(np.mean(np.maximum(strike_bond - p, 0.0) / pn))

    simulate(base, sv, t_star, list(where.keys()), n_paths, seed, observer)
    return out
