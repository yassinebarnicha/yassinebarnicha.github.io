from __future__ import annotations

# Standard-library imports.
import bisect
import math
from dataclasses import dataclass
from typing import Sequence

# Local imports.
from .curve import ZeroCurve
from .instruments import CoterminalSwaption
from .solvers import brent


@dataclass
class PiecewiseConstantSigma:
    """
    Piecewise-constant instantaneous volatility sigma(t).

    bucket_ends = [T1, T2, ..., TN]
    values      = [s1, s2, ..., sN]

    This means:
      sigma(t)=s1 on (0,T1]
      sigma(t)=s2 on (T1,T2]
      ...
    """

    bucket_ends: Sequence[float]
    values: Sequence[float]

    def __post_init__(self) -> None:
        # One volatility value per bucket end.
        if len(self.bucket_ends) != len(self.values):
            raise ValueError("bucket_ends and values must have same length")

        # Bucket ends must be strictly increasing in time.
        if any(t2 <= t1 for t1, t2 in zip(self.bucket_ends[:-1], self.bucket_ends[1:])):
            raise ValueError("bucket_ends must be strictly increasing")

    def sigma(self, t: float) -> float:
        """Return sigma(t)."""
        # For t <= 0 we just return the first bucket value.
        if t <= 0.0:
            return self.values[0]

        # Locate the bucket containing t.
        idx = bisect.bisect_left(self.bucket_ends, t)

        # Cap the index at the last available bucket.
        idx = min(idx, len(self.values) - 1)

        # Return the bucket value.
        return self.values[idx]

    def with_replaced_value(self, idx: int, value: float) -> "PiecewiseConstantSigma":
        """Convenience method used in bootstrap."""
        vals = list(self.values)
        vals[idx] = value
        return PiecewiseConstantSigma(self.bucket_ends, vals)


@dataclass
class HJM1FModel:
    """
    One-factor HJM model in the exact state variable used in the notes:

        X_t = int_0^t sigma(u) exp(lambda u) dW_u

    so that:

        dX_t = sigma(t) exp(lambda t) dW_t = eta(t) dW_t

    Hence X_t has *no drift*, and the PDE contains no first derivative term.
    """

    # Input zero curve P(0,T).
    curve: ZeroCurve

    # Mean-reversion parameter lambda appearing in exp(-lambda(T-t)).
    lam: float

    # Piecewise-constant sigma(t).
    sigma_curve: PiecewiseConstantSigma

    def Bbar(self, t: float) -> float:
        """
        B_lambda(t) from the notes:

            B_lambda(t) = int_0^t exp(-lambda u) du

        closed form = (1 - exp(-lambda t)) / lambda.
        """
        if abs(self.lam) < 1e-14:
            return t
        return (1.0 - math.exp(-self.lam * t)) / self.lam

    def b(self, t: float) -> float:
        """Short notation b(t)=exp(-lambda t)."""
        return math.exp(-self.lam * t)

    def D(self, t: float, T: float) -> float:
        """
        D(t,T)=B_lambda(T)-B_lambda(t).

        This is exactly the coefficient multiplying the factor X_t in the bond
        price formula of the notes.
        """
        return self.Bbar(T) - self.Bbar(t)

    def alpha(self, t: float) -> float:
        """
        Closed-form piecewise-constant evaluation of

            alpha(t) = int_0^t sigma(u)^2 exp(2 lambda u) du.

        This is Var[X_t].
        """
        if t <= 0.0:
            return 0.0

        total = 0.0
        left = 0.0

        # Loop over the piecewise-constant buckets.
        for end, sig in zip(self.sigma_curve.bucket_ends, self.sigma_curve.values):
            # Truncate the bucket to the target horizon t.
            right = min(end, t)

            # Only integrate if the truncated bucket is non-empty.
            if right > left:
                if abs(self.lam) < 1e-14:
                    # Limit lambda -> 0.
                    total += sig * sig * (right - left)
                else:
                    # Exact primitive of exp(2 lambda u).
                    total += sig * sig * (math.exp(2.0 * self.lam * right) - math.exp(2.0 * self.lam * left)) / (2.0 * self.lam)

            # Stop as soon as we have reached t.
            if end >= t:
                break

            # Move to the next bucket.
            left = end

        return total

    def gamma(self, t: float) -> float:
        """
        Closed-form piecewise-constant evaluation of

            gamma(t) = int_0^t sigma(u)^2 exp(2 lambda u) B_lambda(u) du.

        This is another deterministic quantity appearing in the bond formula.
        """
        if t <= 0.0:
            return 0.0

        total = 0.0
        left = 0.0

        # Loop bucket by bucket.
        for end, sig in zip(self.sigma_curve.bucket_ends, self.sigma_curve.values):
            # Truncate the bucket.
            right = min(end, t)

            # Integrate only if non-empty.
            if right > left:
                if abs(self.lam) < 1e-14:
                    # Limit lambda -> 0: B_lambda(u) -> u.
                    total += 0.5 * sig * sig * (right * right - left * left)
                else:
                    # Primitive of exp(2 lambda u) * B_lambda(u).
                    total += (sig * sig / self.lam) * (
                        (math.exp(2.0 * self.lam * right) - math.exp(2.0 * self.lam * left)) / (2.0 * self.lam)
                        - (math.exp(self.lam * right) - math.exp(self.lam * left)) / self.lam
                    )

            # Stop once the horizon t has been covered.
            if end >= t:
                break

            # Move to the next bucket.
            left = end

        return total

    def eta(self, t: float) -> float:
        """
        Instantaneous volatility of the state variable X_t:

            dX_t = eta(t) dW_t
            eta(t) = sigma(t) exp(lambda t)
        """
        return self.sigma_curve.sigma(t) * math.exp(self.lam * t)

    def discount_bond(self, t: float, T: float, x_t: float) -> float:
        """
        Bond price formula from the notes:

            P(t,T) = P(0,T)/P(0,t) * exp( D(t,T) * ( c(t,T) - X_t ) )

        where:
            D(t,T) = B_lambda(T) - B_lambda(t)
            c(t,T) = 0.5 * (B_lambda(T)+B_lambda(t)) * alpha(t) - gamma(t)
        """
        # By definition P(t,T)=1 for T<=t.
        if T <= t:
            return 1.0

        # Coefficient multiplying the state variable X_t.
        d = self.D(t, T)

        # Deterministic term in the exponent.
        c = 0.5 * (self.Bbar(T) + self.Bbar(t)) * self.alpha(t) - self.gamma(t)

        # Full exponent.
        expo = d * (c - x_t)

        # Final bond price.
        return (self.curve.discount(T) / self.curve.discount(t)) * math.exp(expo)

    def short_rate(self, t: float, x_t: float) -> float:
        """
        Recover r_t consistently from the bond formula.

        From the notes, taking T down to t yields:

            r_t = f(0,t) - b(t) * ( B_lambda(t) * alpha(t) - gamma(t) - X_t )
        """
        # Initial instantaneous forward curve f(0,t).
        f0 = self.curve.inst_fwd(t)

        # Short-rate formula consistent with the HJM state variable X_t.
        return f0 - self.b(t) * (self.Bbar(t) * self.alpha(t) - self.gamma(t) - x_t)

    def swap_value(self, swaption: CoterminalSwaption, x_expiry: float) -> float:
        """
        Value at the swaption expiry of the underlying coterminal swap as a
        function of the state X_expiry.
        """
        # Remaining fixed-leg schedule after the expiry date.
        times = swaption.remaining_payment_times()
        accruals = swaption.remaining_accruals()

        # Empty schedule -> zero value.
        if not times:
            return 0.0

        # Discount factor to the terminal date of the swap.
        pn = self.discount_bond(swaption.expiry, swaption.terminal, x_expiry)

        # Fixed-leg annuity valued at the expiry date.
        ann = sum(a * self.discount_bond(swaption.expiry, t, x_expiry) for t, a in zip(times, accruals))

        # Payer swap value under the single-curve simplification.
        payer_val = 1.0 - pn - swaption.strike * ann

        # Receiver is the opposite sign.
        return payer_val if swaption.payer else -payer_val

    def zcb_option_vol(self, option_expiry: float, bond_maturity: float) -> float:
        """
        Total volatility entering the bond-option formula.

        Because the bond is exponential-affine in X_t and Var[X_t]=alpha(t),
        the total bond-option standard deviation is:

            |D(expiry,maturity)| * sqrt(alpha(expiry))
        """
        d = self.D(option_expiry, bond_maturity)
        return abs(d) * math.sqrt(max(self.alpha(option_expiry), 0.0))

    def zcb_put(self, option_expiry: float, bond_maturity: float, strike: float) -> float:
        """
        Put on a zero-coupon bond.

        This is the analytical building block used by the Jamshidian
        decomposition for the European swaption.
        """
        # Immediate-expiry case.
        if option_expiry <= 0.0:
            return max(strike * self.curve.discount(option_expiry) - self.curve.discount(bond_maturity), 0.0)

        # Market bond prices today.
        p0t = self.curve.discount(option_expiry)
        p0T = self.curve.discount(bond_maturity)

        # Total lognormal-style bond volatility in this model.
        v = self.zcb_option_vol(option_expiry, bond_maturity)

        # Near-deterministic case.
        if v < 1e-14:
            return max(strike * p0t - p0T, 0.0)

        # Standard h variable for the closed-form bond-option formula.
        h = math.log(p0T / (strike * p0t)) / v + 0.5 * v

        # Gaussian probabilities.
        n1 = 0.5 * (1.0 + math.erf(-h / math.sqrt(2.0)))
        n2 = 0.5 * (1.0 + math.erf((v - h) / math.sqrt(2.0)))

        # Bond put price.
        return strike * p0t * n2 - p0T * n1

    def european_swaption_price_jamshidian(self, swaption: CoterminalSwaption) -> float:
        """
        Semi-analytical European swaption price via Jamshidian.

        'Semi-analytical' because the decomposition itself is analytical, but we
        first need to solve numerically for the unique root x* such that the
        underlying swap value at expiry is zero.
        """
        # Remaining schedule of the underlying coterminal swap.
        times = swaption.remaining_payment_times()
        accruals = swaption.remaining_accruals()

        # Empty schedule -> zero value.
        if not times:
            return 0.0

        # Root function for Jamshidian.
        def objective(x: float) -> float:
            return self.swap_value(swaption, x)

        # Solve for x* such that swap_value(expiry, x*) = 0.
        x_star = brent(objective, -10.0, 10.0, tol=1e-12, max_iter=200)

        # Start the decomposition sum.
        price = 0.0

        # All coupons except the last one contribute K * delta_i * bond-put.
        for t_i, delta_i in zip(times[:-1], accruals[:-1]):
            # Bond strike induced by x*.
            k_i = self.discount_bond(swaption.expiry, t_i, x_star)

            # Add the corresponding bond-put contribution.
            price += swaption.strike * delta_i * self.zcb_put(swaption.expiry, t_i, k_i)

        # The last cashflow also contains the notional repayment term.
        last_t = times[-1]
        last_delta = accruals[-1]
        k_last = self.discount_bond(swaption.expiry, last_t, x_star)
        price += (1.0 + swaption.strike * last_delta) * self.zcb_put(swaption.expiry, last_t, k_last)

        # Receiver handling is not fully specialized here; the present project is
        # focused on payer Bermudans.
        return price if swaption.payer else price
