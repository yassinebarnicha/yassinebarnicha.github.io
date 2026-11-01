from __future__ import annotations

import math
from typing import List, Tuple

from hjm1f.curve import ZeroCurve
from hjm1f.solvers import brent

from .hjm_2f import HJM2FModel


def _phi(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def cap_fixings(maturity: float, delta: float = 0.5) -> List[float]:
    """Fixing dates of a cap on a 6-month index: delta, 2 delta, ..., maturity - delta (the first period is already fixed)."""
    n = int(round(maturity / delta))
    return [delta * i for i in range(1, n)]


def atm_cap_strike(curve: ZeroCurve, maturity: float, delta: float = 0.5) -> float:
    """At-the-money strike of a cap: the swap rate of its schedule."""
    n = int(round(maturity / delta))
    annuity = sum(delta * curve.discount(delta * i) for i in range(1, n + 1))
    return (1.0 - curve.discount(maturity)) / annuity


def bachelier_cap_price(curve: ZeroCurve, maturity: float, strike: float, vol: float, delta: float = 0.5) -> float:
    """Cap price (per unit notional) at a flat normal volatility: sum of Bachelier caplets on the forward index."""
    total = 0.0
    for t in cap_fixings(maturity, delta):
        p0, p1 = curve.discount(t), curve.discount(t + delta)
        fwd = (p0 / p1 - 1.0) / delta
        s = vol * math.sqrt(t)
        d = (fwd - strike) / s
        total += delta * p1 * ((fwd - strike) * _phi(d) + s * _pdf(d))
    return total


def model_cap_price(model: HJM2FModel, maturity: float, strike: float, delta: float = 0.5) -> float:
    """
    Cap in the model: each caplet delta (L - K)^+ paid at T + delta is (1 + delta K) times a put, struck at
    1/(1 + delta K), on the zero-coupon bond P(T, T + delta). Closed form with one or two factors.
    """
    k = 1.0 / (1.0 + delta * strike)
    return sum((1.0 + delta * strike) * model.zcb_put(t, t + delta, k) for t in cap_fixings(maturity, delta))


def implied_cap_vol(curve: ZeroCurve, maturity: float, strike: float, price: float, delta: float = 0.5) -> float:
    """Flat normal volatility of the cap that reproduces `price`."""
    return brent(lambda v: bachelier_cap_price(curve, maturity, strike, v, delta) - price, 1e-6, 0.05, tol=1e-12, max_iter=200)


def model_cap_vol(model: HJM2FModel, maturity: float, delta: float = 0.5) -> float:
    """Model ATM cap vol (flat normal vol of the cap struck at the swap rate of its schedule)."""
    k = atm_cap_strike(model.curve, maturity, delta)
    price = model_cap_price(model, maturity, k, delta)
    return implied_cap_vol(model.curve, maturity, k, price, delta) if price > 0 else 0.0
