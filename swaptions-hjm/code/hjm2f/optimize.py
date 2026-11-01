from __future__ import annotations

from typing import Callable, List, Sequence, Tuple


def nelder_mead(f: Callable[[Sequence[float]], float], x0: Sequence[float], step: float = 0.5,
                max_fev: int = 300, xtol: float = 1e-4, ftol: float = 1e-4) -> Tuple[List[float], float, int]:
    """Nelder-Mead simplex (reflection 1, expansion 2, contraction 1/2, shrink 1/2). Returns (x, f(x), evaluations)."""
    n = len(x0)
    simplex = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += step
        simplex.append(p)
    values = [f(p) for p in simplex]
    nfev = len(simplex)
    while nfev < max_fev:
        order = sorted(range(n + 1), key=lambda i: values[i])
        simplex = [simplex[i] for i in order]
        values = [values[i] for i in order]
        size = max(abs(simplex[i][k] - simplex[0][k]) for i in range(1, n + 1) for k in range(n))
        if size < xtol and abs(values[-1] - values[0]) < ftol:
            break
        centre = [sum(simplex[i][k] for i in range(n)) / n for k in range(n)]
        worst = simplex[-1]
        refl = [centre[k] + (centre[k] - worst[k]) for k in range(n)]
        f_refl = f(refl)
        nfev += 1
        if f_refl < values[0]:
            expd = [centre[k] + 2.0 * (centre[k] - worst[k]) for k in range(n)]
            f_exp = f(expd)
            nfev += 1
            if f_exp < f_refl:
                simplex[-1], values[-1] = expd, f_exp
            else:
                simplex[-1], values[-1] = refl, f_refl
        elif f_refl < values[-2]:
            simplex[-1], values[-1] = refl, f_refl
        else:
            outside = f_refl < values[-1]
            if outside:
                point = [centre[k] + 0.5 * (refl[k] - centre[k]) for k in range(n)]
            else:
                point = [centre[k] + 0.5 * (worst[k] - centre[k]) for k in range(n)]
            f_point = f(point)
            nfev += 1
            if f_point < min(f_refl, values[-1]):
                simplex[-1], values[-1] = point, f_point
            else:
                for i in range(1, n + 1):
                    simplex[i] = [simplex[0][k] + 0.5 * (simplex[i][k] - simplex[0][k]) for k in range(n)]
                    values[i] = f(simplex[i])
                    nfev += 1
    best = min(range(n + 1), key=lambda i: values[i])
    return simplex[best], values[best], nfev
