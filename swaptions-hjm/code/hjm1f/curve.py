from __future__ import annotations

# Standard-library imports.
import bisect
import math
from dataclasses import dataclass
from typing import Sequence


# Small helper: linear interpolation.
# We use it both for the discount curve itself and for log-discounts.
def _linear_interp(x: float, xs: Sequence[float], ys: Sequence[float]) -> float:
    # Flat extrapolation to the left.
    if x <= xs[0]:
        return ys[0]

    # Flat extrapolation to the right.
    if x >= xs[-1]:
        return ys[-1]

    # Find the interval [xs[i], xs[i+1]] containing x.
    i = bisect.bisect_right(xs, x) - 1

    # Endpoints of the interval.
    x0, x1 = xs[i], xs[i + 1]
    y0, y1 = ys[i], ys[i + 1]

    # Linear interpolation weight.
    w = (x - x0) / (x1 - x0)

    # Return the interpolated value.
    return y0 * (1.0 - w) + y1 * w


@dataclass
class ZeroCurve:
    """
    Zero-coupon curve supplied directly as discount factors P(0,T).

    We intentionally keep this class minimal: the project assumes the curve is
    already available as an input, so we do not build it from deposits/swaps.
    """

    # Curve pillar times in years.
    times: Sequence[float]

    # Corresponding discount factors P(0,T).
    discounts: Sequence[float]

    def __post_init__(self) -> None:
        # Sanity check: same number of times and discounts.
        if len(self.times) != len(self.discounts):
            raise ValueError("times and discounts must have same length")

        # Sanity check: times must be strictly increasing.
        if any(t2 <= t1 for t1, t2 in zip(self.times[:-1], self.times[1:])):
            raise ValueError("curve times must be strictly increasing")

        # If time 0 is missing, prepend it explicitly.
        if self.times[0] > 0.0:
            self.times = [0.0, *self.times]
            self.discounts = [1.0, *self.discounts]

        # At time 0 the discount factor must be exactly 1.
        if self.discounts[0] != 1.0:
            raise ValueError("discount at time 0 must be 1.0")

    def discount(self, t: float) -> float:
        """
        Return P(0,t).

        We interpolate in log-discounts rather than in discounts directly.
        This is a standard pragmatic choice because it preserves positivity.
        """
        # By definition P(0,0)=1.
        if t <= 0.0:
            return 1.0

        # Convert the stored discounts into log-discounts.
        logs = [math.log(max(d, 1e-16)) for d in self.discounts]

        # Interpolate the log-discount and exponentiate back.
        return math.exp(_linear_interp(t, self.times, logs))

    def inst_fwd(self, t: float, h: float = 1e-5) -> float:
        """
        Approximate the instantaneous forward rate f(0,t).

        Numerically:
            f(0,t) = - d/dT log P(0,T) |_{T=t}

        We use a symmetric finite-difference for simplicity.
        """
        # Left point for the finite difference. We do not go below 0.
        t0 = max(0.0, t - h)

        # Right point for the finite difference.
        t1 = t + h

        # Discount factors at the two nearby maturities.
        p0 = self.discount(t0)
        p1 = self.discount(t1)

        # Central finite-difference on log-discounts.
        return -(math.log(p1) - math.log(p0)) / (t1 - t0)
