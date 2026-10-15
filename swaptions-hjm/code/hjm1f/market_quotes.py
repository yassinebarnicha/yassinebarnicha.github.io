from __future__ import annotations

# Standard-library imports.
import math

# Local imports.
from .curve import ZeroCurve
from .instruments import CoterminalSwaption


# Standard normal density.
def _pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


# Standard normal cumulative distribution function.
def _phi(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


# Swap annuity A(0).
def annuity(curve: ZeroCurve, payment_times: list[float], accruals: list[float]) -> float:
    return sum(a * curve.discount(t) for t, a in zip(payment_times, accruals))


# Forward swap rate S(0) under the usual single-curve simplification used in the notes.
def forward_swap_rate(
    curve: ZeroCurve,
    expiry: float,
    terminal: float,
    payment_times: list[float],
    accruals: list[float],
) -> float:
    # Terminal discount factor P(0,T_N).
    p_terminal = curve.discount(terminal)

    # Annuity of the fixed leg.
    ann = annuity(curve, payment_times, accruals)

    # If the annuity is numerically zero, the forward is undefined.
    if ann <= 0.0:
        raise ValueError("annuity must be positive")

    # Standard single-curve formula.
    return (curve.discount(expiry) - p_terminal) / ann


# Bachelier price from a normal vol quote.
def bachelier_swaption_price(curve: ZeroCurve, swaption: CoterminalSwaption, normal_vol: float) -> float:
    # Remaining fixed-leg schedule of the underlying coterminal swap.
    times = swaption.remaining_payment_times()
    accruals = swaption.remaining_accruals()

    # Annuity A(0).
    ann = annuity(curve, times, accruals)

    # Forward swap rate S(0).
    fwd = forward_swap_rate(curve, swaption.expiry, swaption.terminal, times, accruals)

    # Strike K.
    k = swaption.strike

    # Time to expiry.
    tau = max(swaption.expiry, 1e-12)

    # Total normal standard deviation sigma_N * sqrt(tau).
    s = max(normal_vol, 0.0) * math.sqrt(tau)

    # Degenerate case: zero volatility -> intrinsic value.
    if s < 1e-16:
        intrinsic = max(fwd - k, 0.0) if swaption.payer else max(k - fwd, 0.0)
        return ann * intrinsic

    # Bachelier d parameter.
    d = (fwd - k) / s

    # Payer formula.
    if swaption.payer:
        return ann * ((fwd - k) * _phi(d) + s * _pdf(d))

    # Receiver formula.
    return ann * ((k - fwd) * _phi(-d) + s * _pdf(d))
