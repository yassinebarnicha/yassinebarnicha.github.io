from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import numpy as np

from hjm1f.calibration_1f import CalibrationInstrument, CoterminalBootstrap1F
from hjm1f.curve import ZeroCurve
from hjm1f.hjm_1f import HJM1FModel
from hjm1f.instruments import CoterminalSwaption

from .calibration_2f import CoterminalBootstrap2F, implied_normal_vol
from .hjm_2f import HJM2FModel
from .optimize import nelder_mead

# order of the structural parameters: lambda_1, lambda_2, rho, sigma_2 / sigma_1
LAM1_MAX, LAM2_MAX = 0.5, 3.0


def to_params(x: Sequence[float]) -> Tuple[Tuple[float, float], float, float]:
    """Unconstrained coordinates -> (lambda_1, lambda_2), rho, ratio (logistic / tanh / exp transforms)."""
    lam1 = LAM1_MAX / (1.0 + math.exp(-x[0]))
    lam2 = LAM2_MAX / (1.0 + math.exp(-x[1]))
    return (lam1, lam2), math.tanh(x[2]), math.exp(x[3])


def from_params(lam: Tuple[float, float], rho: float, ratio: float) -> List[float]:
    logit = lambda p, top: math.log((p / top) / (1.0 - p / top))
    eps = 1e-9
    return [logit(min(max(lam[0], eps), LAM1_MAX - eps), LAM1_MAX), logit(min(max(lam[1], eps), LAM2_MAX - eps), LAM2_MAX),
            math.atanh(min(max(rho, -0.999999), 0.999999)), math.log(ratio)]


@dataclass
class MatrixFit:
    lam: Tuple[float, float]
    rho: float
    ratio: float
    rmse_bp: float
    errors_bp: List[float]
    sigma1: List[float]
    evaluations: int
    history: List[float] = field(default_factory=list)


def _errors_2f(curve, lam, rho, ratio, exercise_times, instruments, targets, n_grid) -> Tuple[List[float], List[float]]:
    res = CoterminalBootstrap2F(curve, lam, rho, ratio, exercise_times).calibrate(instruments, n_grid=n_grid)
    model = HJM2FModel(curve, lam, res.sigma1, res.sigma2, rho)
    errs = []
    for sw, vol in targets:
        price = model.european_swaption_price(sw, n_grid=n_grid)
        errs.append((implied_normal_vol(curve, sw, price) - vol) * 1e4 if price > 0 else -vol * 1e4)
    return errs, list(res.sigma1.values)


def fit_two_factor(curve: ZeroCurve, exercise_times: Sequence[float], instruments: Sequence[CalibrationInstrument],
                   targets: Sequence[Tuple[CoterminalSwaption, float]], x0: Sequence[float],
                   free: Sequence[bool] = (True, True, True, True), n_grid: int = 61, max_fev: int = 200) -> MatrixFit:
    """
    Global fit of (lambda_1, lambda_2, rho, sigma_2/sigma_1) to a set of European swaptions `targets`
    (list of (swaption, normal vol)). For every candidate the volatility level sigma_1(t) is re-bootstrapped
    exactly on the coterminal `instruments` (the ones of the Bermudan), so that the Bermudan's own Europeans
    are always matched and the free parameters only have to explain the other swaptions.
    """
    x0 = list(x0)
    idx = [i for i, f in enumerate(free) if f]
    history: List[float] = []

    def objective(y: Sequence[float]) -> float:
        x = list(x0)
        for k, i in enumerate(idx):
            x[i] = y[k]
        lam, rho, ratio = to_params(x)
        try:
            errs, _ = _errors_2f(curve, lam, rho, ratio, exercise_times, instruments, targets, n_grid)
        except Exception:
            return 1e3
        rmse = float(np.sqrt(np.mean(np.square(errs))))
        history.append(rmse)
        return rmse

    y, fy, nfev = nelder_mead(objective, [x0[i] for i in idx], max_fev=max_fev)
    x = list(x0)
    for k, i in enumerate(idx):
        x[i] = y[k]
    lam, rho, ratio = to_params(x)
    errs, sig = _errors_2f(curve, lam, rho, ratio, exercise_times, instruments, targets, n_grid)
    return MatrixFit(lam, rho, ratio, float(np.sqrt(np.mean(np.square(errs)))), errs, sig, nfev, history)


def fit_one_factor(curve: ZeroCurve, exercise_times: Sequence[float], instruments: Sequence[CalibrationInstrument],
                   targets: Sequence[Tuple[CoterminalSwaption, float]], lam_max: float = 0.5) -> Tuple[float, float, List[float]]:
    """
    One-factor counterpart: the best mean reversion lambda in [0, lam_max] (golden-section search), the volatility
    bootstrapped on the coterminals for each lambda. Returns (lambda, rmse in bp, errors in bp).
    """
    def errors(lam: float) -> List[float]:
        calib = CoterminalBootstrap1F(curve, lam, exercise_times).calibrate(instruments)
        model = HJM1FModel(curve, lam, calib.sigma_curve)
        errs = []
        for sw, vol in targets:
            price = model.european_swaption_price_jamshidian(sw)
            errs.append((implied_normal_vol(curve, sw, price) - vol) * 1e4 if price > 0 else -vol * 1e4)
        return errs

    rmse = lambda lam: float(np.sqrt(np.mean(np.square(errors(lam)))))
    a, b = 0.0, lam_max
    g = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = rmse(c), rmse(d)
    for _ in range(30):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = rmse(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = rmse(d)
    best = min([(rmse(a), a), (fc, c), (fd, d), (rmse(b), b)])
    return best[1], best[0], errors(best[1])
