from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List

import numpy as np

from .hjm_1f import HJM1FModel
from .instruments import BermudanCoterminalSwaption


def solve_tridiagonal(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Thomas algorithm for a tridiagonal linear system."""
    n = len(d)
    cp = np.zeros(n - 1)
    dp = np.zeros(n)
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n - 1):
        den = b[i] - a[i - 1] * cp[i - 1]
        cp[i] = c[i] / den
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / den
    dp[n - 1] = (d[n - 1] - a[n - 2] * dp[n - 2]) / (b[n - 1] - a[n - 2] * cp[n - 2])
    x = np.zeros(n)
    x[n - 1] = dp[n - 1]
    for i in reversed(range(n - 1)):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


@dataclass
class PDEGrid:
    x_min: float
    x_max: float
    nx: int
    nt: int
    t_max: float

    @staticmethod
    def around_factor_std(model: HJM1FModel, t_max: float, n_std: float = 5.0, nx: int = 301, nt: int = 500) -> "PDEGrid":
        var = max(model.alpha(t_max), 1e-16)
        std = math.sqrt(var)
        return PDEGrid(x_min=-n_std * std, x_max=n_std * std, nx=nx, nt=nt, t_max=t_max)


class BermudanPDE1F:
    """
    1D backward PDE pricer in the driftless HJM state variable X_t.

    Key fixes relative to the first version:
    1) exercise dates are inserted *exactly* in the time grid;
    2) boundary conditions are no longer built from a wrong underlying-at-time-t proxy;
    3) early exercise is applied exactly at the Bermudan dates.

    Time scheme: theta = 1 is implicit Euler (first order); theta = 0.5 is Crank-Nicolson (second order),
    for which `rannacher` implicit steps after the final date and after every exercise date damp the
    oscillations caused by the kink of the payoff.
    """

    def __init__(self, model: HJM1FModel, bermudan: BermudanCoterminalSwaption, grid: PDEGrid,
                 theta: float = 1.0, rannacher: int = 0):
        self.model = model
        self.bermudan = bermudan
        self.grid = grid
        self.theta = theta
        self.rannacher = rannacher

    def _exercise_values(self, t: float, xs: np.ndarray) -> np.ndarray:
        und = self.bermudan.underlying_at(t)
        return np.array([max(self.model.swap_value(und, float(x)), 0.0) for x in xs])

    def _future_exercise_times(self, t: float) -> List[float]:
        eps = 1e-12
        return [u for u in self.bermudan.exercise_times if u >= t - eps]

    def _boundary_values(self, t: float, v_next: np.ndarray, xs: np.ndarray) -> tuple[float, float]:
        """
        Conservative boundary treatment.

        We keep the continuation value inherited from the next time slice and also
        enforce that the boundary value cannot fall below the best intrinsic value
        over the remaining exercise dates.
        """
        left = float(v_next[0])
        right = float(v_next[-1])

        future_ex = self._future_exercise_times(t)
        if future_ex:
            left_intr = 0.0
            right_intr = 0.0
            x_left = float(xs[0])
            x_right = float(xs[-1])
            for ex in future_ex:
                und = self.bermudan.underlying_at(ex)
                left_intr = max(left_intr, max(self.model.swap_value(und, x_left), 0.0))
                right_intr = max(right_intr, max(self.model.swap_value(und, x_right), 0.0))
            left = max(left, left_intr)
            right = max(right, right_intr)
        return left, right

    def _build_time_grid(self) -> np.ndarray:
        """
        Build a grid that contains all Bermudan exercise dates exactly.

        Between two consecutive key dates, we insert a number of uniform substeps
        proportional to the interval length so that the overall resolution remains
        close to the user-requested nt.
        """
        key_times = sorted({0.0, self.grid.t_max, *self.bermudan.exercise_times})
        if len(key_times) < 2:
            return np.array([0.0, self.grid.t_max], dtype=float)

        total_horizon = key_times[-1] - key_times[0]
        base_dt = total_horizon / max(self.grid.nt, 1)
        times: List[float] = [key_times[0]]
        for t0, t1 in zip(key_times[:-1], key_times[1:]):
            interval = t1 - t0
            n_sub = max(1, int(round(interval / base_dt)))
            local = np.linspace(t0, t1, n_sub + 1)
            times.extend(float(x) for x in local[1:])
        return np.array(times, dtype=float)

    def price(self) -> float:
        xs = np.linspace(self.grid.x_min, self.grid.x_max, self.grid.nx)
        dx = xs[1] - xs[0]
        ts = self._build_time_grid()
        exercise_set = {round(t, 12) for t in self.bermudan.exercise_times}
        model = self.model

        last_ex = self.bermudan.exercise_times[-1]
        v = self._exercise_values(last_ex, xs)
        smoothing = self.rannacher

        for n in range(len(ts) - 2, -1, -1):
            t = float(ts[n])
            t1 = float(ts[n + 1])
            dt = t1 - t
            theta = 1.0 if smoothing > 0 else self.theta
            if smoothing > 0:
                smoothing -= 1

            # the short rate is affine in x: r(t, x) = A(t) + b(t) x
            a_t = model.short_rate(t, 0.0)
            b_t = model.b(t)
            m = self.grid.nx - 2
            v_in = v[1:-1]

            if self.theta == 1.0:
                eta = model.eta(t)
                diff = 0.5 * eta * eta / (dx * dx)
                diag = 1.0 + 2.0 * dt * diff + dt * (a_t + b_t * xs[1:-1])
                rhs = v_in.copy()
            else:
                # every coefficient at its own end of the step; sigma read in the middle of the step
                a_1 = model.short_rate(t1, 0.0)
                b_1 = model.b(t1)
                s = model.sigma_curve.sigma(0.5 * (t + t1))
                e0 = s * math.exp(model.lam * t)
                e1 = s * math.exp(model.lam * t1)
                diff = 0.5 * e0 * e0 / (dx * dx)
                diff1 = 0.5 * e1 * e1 / (dx * dx)
                diag = 1.0 + theta * dt * (2.0 * diff + a_t + b_t * xs[1:-1])
                rhs = v_in - (1.0 - theta) * dt * ((2.0 * diff1 + a_1 + b_1 * xs[1:-1]) * v_in - diff1 * (v[:-2] + v[2:]))
            lower = np.full(m - 1, -theta * dt * diff)
            upper = np.full(m - 1, -theta * dt * diff)

            left_bc, right_bc = self._boundary_values(t, v, xs)
            rhs[0] += theta * dt * diff * left_bc
            rhs[-1] += theta * dt * diff * right_bc

            inner = solve_tridiagonal(lower, diag, upper, rhs)
            v_new = np.zeros_like(v)
            v_new[0] = left_bc
            v_new[-1] = right_bc
            v_new[1:-1] = inner

            if round(t, 12) in exercise_set:
                v_new = np.maximum(v_new, self._exercise_values(t, xs))
                smoothing = self.rannacher

            v = v_new

        return float(np.interp(0.0, xs, v))
