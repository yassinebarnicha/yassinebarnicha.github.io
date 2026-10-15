from __future__ import annotations

# Standard-library imports.
from dataclasses import dataclass
from typing import List


@dataclass
class CoterminalSwaption:
    """
    European coterminal swaption.

    expiry: exercise date of the swaption
    terminal: final maturity of the underlying swap
    strike: fixed rate K
    fixed_payment_times: fixed leg payment dates of the full terminal swap
    fixed_accruals: corresponding accrual factors
    payer: True for payer, False for receiver
    """

    expiry: float
    terminal: float
    strike: float
    fixed_payment_times: List[float]
    fixed_accruals: List[float]
    payer: bool = True

    def remaining_payment_times(self) -> List[float]:
        """Keep only coupons strictly after the expiry date."""
        return [t for t in self.fixed_payment_times if t > self.expiry]

    def remaining_accruals(self) -> List[float]:
        """Keep the accruals matching the remaining payment dates."""
        return [a for t, a in zip(self.fixed_payment_times, self.fixed_accruals) if t > self.expiry]


@dataclass
class BermudanCoterminalSwaption:
    """
    Bermudan coterminal swaption.

    At each exercise date T_k, the holder may enter a swap ending at the same
    final maturity 'terminal'. This is exactly the coterminal Bermudan setup.
    """

    exercise_times: List[float]
    terminal: float
    strike: float
    fixed_payment_times: List[float]
    fixed_accruals: List[float]
    payer: bool = True

    def underlying_at(self, exercise_time: float) -> CoterminalSwaption:
        """
        Build the European coterminal swaption corresponding to one Bermudan
        exercise date.
        """
        return CoterminalSwaption(
            expiry=exercise_time,
            terminal=self.terminal,
            strike=self.strike,
            fixed_payment_times=self.fixed_payment_times,
            fixed_accruals=self.fixed_accruals,
            payer=self.payer,
        )
