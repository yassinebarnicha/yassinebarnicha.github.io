from __future__ import annotations

# Standard-library imports.
import math
from typing import Callable


# A genuine Brent solver, used for the Jamshidian root and for bootstrap.
def brent(f: Callable[[float], float], a: float, b: float, tol: float = 1e-10, max_iter: int = 100) -> float:
    # Evaluate the function at the bracket endpoints.
    fa = f(a)
    fb = f(b)

    # If one endpoint is already a root, return it immediately.
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b

    # Brent requires a bracket.
    if fa * fb > 0.0:
        raise ValueError("root must be bracketed")

    # Initialize the third point and step sizes.
    c, fc = a, fa
    d = e = b - a

    # Main Brent iteration loop.
    for _ in range(max_iter):
        # If b and c are on the same side, reset c to a.
        if fb * fc > 0.0:
            c, fc = a, fa
            d = e = b - a

        # Keep |fb| <= |fc| so that b is the best current approximation.
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb

        # Numerical tolerance for the current iterate.
        tol1 = 2.0 * math.ulp(1.0) * abs(b) + 0.5 * tol

        # Midpoint of the current bracket.
        xm = 0.5 * (c - b)

        # Stopping criterion.
        if abs(xm) <= tol1 or fb == 0.0:
            return b

        # Attempt inverse quadratic interpolation / secant if safe.
        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                # Secant step.
                p = 2.0 * xm * s
                q = 1.0 - s
            else:
                # Inverse quadratic interpolation step.
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)

            # Adjust signs so that p/q gives the proposed step.
            if p > 0.0:
                q = -q
            p = abs(p)

            # Accept interpolation only if it stays conservative.
            min1 = 3.0 * xm * q - abs(tol1 * q)
            min2 = abs(e * q)
            if 2.0 * p < min(min1, min2):
                e = d
                d = p / q
            else:
                # Fallback to bisection.
                d = xm
                e = d
        else:
            # No safe interpolation available -> bisection.
            d = xm
            e = d

        # Move the old best point into a.
        a, fa = b, fb

        # Update b by the chosen step, while respecting tolerance.
        if abs(d) > tol1:
            b += d
        else:
            b += tol1 if xm >= 0.0 else -tol1

        # Re-evaluate the function.
        fb = f(b)

    # If we reach this point, the iteration budget was not enough.
    raise RuntimeError("Brent did not converge")
