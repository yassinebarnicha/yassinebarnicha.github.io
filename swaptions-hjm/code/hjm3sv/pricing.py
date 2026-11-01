"""
Calibration and pricing for the stochastic-volatility model, always with the exact 2F price as a control variate:

    price_SV = price_2F(exact) + ( MC_SV - MC_{xi=0} )      (same random numbers in both Monte Carlo runs)

so that the Monte Carlo noise only has to resolve the (small) effect of the stochastic volatility.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from hjm1f.instruments import BermudanCoterminalSwaption, CoterminalSwaption
from hjm1f.market_quotes import forward_swap_rate
from hjm2f.calibration_2f import atm_swaption, implied_normal_vol
from hjm2f.caps import atm_cap_strike, cap_fixings, implied_cap_vol, model_cap_price
from hjm2f.hjm_2f import HJM2FModel
from hjm2f.market_data import MarketData
from hjm2f.optimize import nelder_mead

from .model import SVParams, State, bond, europeans, simulate, swap_value

ZERO = SVParams(1.0, 0.0, 0.0)
KAPPA_RANGE, XI_RANGE = (0.02, 3.0), (0.02, 2.5)
RHO_V_MAX = 0.95
NAN_PENALTY_BP = 100.0   # a strike whose model price is not positive (no implied vol) counts as a 100 bp miss


def sv_from_x(x: Sequence[float]) -> SVParams:
    lg = lambda v, lo, hi: lo + (hi - lo) / (1.0 + math.exp(-v))
    return SVParams(lg(x[0], *KAPPA_RANGE), lg(x[1], *XI_RANGE), RHO_V_MAX * math.tanh(x[2]))


def x_from_sv(sv: SVParams) -> List[float]:
    lt = lambda p, lo, hi: math.log(min(max((p - lo) / (hi - lo), 1e-6), 1 - 1e-6) / (1 - min(max((p - lo) / (hi - lo), 1e-6), 1 - 1e-6)))
    return [lt(sv.kappa, *KAPPA_RANGE), lt(sv.xi, *XI_RANGE), math.atanh(min(max(sv.rho_v / RHO_V_MAX, -0.999), 0.999))]


# ---------------------------------------------------------------- smiles
@dataclass
class SmilePoint:
    expiry: float
    tenor: int
    offset_bp: float          # 0 for the ATM anchor
    swaption: CoterminalSwaption
    market_vol: float         # nan for the anchor
    atm_vol: float            # market ATM vol of the smile


def smile_points(md: MarketData, pairs: Sequence[Tuple[float, int]], z_max: float = 3.2) -> List[SmilePoint]:
    """Strikes more than z_max ATM standard deviations away are left out: their prices are below the Monte Carlo resolution."""
    pts: List[SmilePoint] = []
    for p in md.raw["swaptionSmile"]["pairs"]:
        e, t = float(p["expiry"]), int(p["tenor"])
        if (e, t) not in pairs or not md.has(e, t):
            continue
        atm = md.swaption_vol(e, t)
        pay = [e + i for i in range(1, t + 1)]
        acc = [1.0] * t
        fwd = forward_swap_rate(md.curve, e, pay[-1], pay, acc)
        pts.append(SmilePoint(e, t, 0.0, CoterminalSwaption(e, pay[-1], fwd, pay, acc, True), float("nan"), atm))
        for off, vol in zip(md.raw["swaptionSmile"]["offsetsBp"], p["vol"]):
            if vol is None:
                continue
            if abs(off) / (atm * 1e4 * math.sqrt(e)) > z_max:
                continue
            k = fwd + off / 1e4
            pts.append(SmilePoint(e, t, float(off), CoterminalSwaption(e, pay[-1], k, pay, acc, off >= 0), float(vol), atm))
    return pts


def _safe_iv(curve, sw, price) -> float:
    try:
        return implied_normal_vol(curve, sw, price) if price > 0 else float("nan")
    except Exception:
        return float("nan")


class SmileEngine:
    """Model smiles of a 2F base model with a stochastic-volatility multiplier, by control-variate Monte Carlo."""

    def __init__(self, base: HJM2FModel, points: Sequence[SmilePoint], t_star: float, n_paths: int = 40000, seed: int = 21, n_grid: int = 121):
        self.base, self.points, self.t_star, self.n, self.seed = base, list(points), t_star, n_paths, seed
        sws = [p.swaption for p in self.points]
        self.exact = np.array([base.european_swaption_price(s, n_grid=n_grid) for s in sws])
        self.mc0 = europeans(base, ZERO, sws, t_star, n_paths, seed)

    def vols(self, sv: SVParams) -> np.ndarray:
        sws = [p.swaption for p in self.points]
        mc = europeans(self.base, sv, sws, self.t_star, self.n, self.seed)
        price = self.exact + (mc - self.mc0)
        return np.array([_safe_iv(self.base.curve, s, p) for s, p in zip(sws, price)])

    def errors_bp(self, sv: SVParams, shape_only: bool = True) -> List[Tuple[SmilePoint, float]]:
        v = self.vols(sv)
        atm_model: Dict[Tuple[float, int], float] = {(p.expiry, p.tenor): v[i] for i, p in enumerate(self.points) if p.offset_bp == 0.0}
        out = []
        for i, p in enumerate(self.points):
            if p.offset_bp == 0.0:
                continue
            m = v[i]
            if shape_only:
                m = p.atm_vol * v[i] / atm_model[(p.expiry, p.tenor)]
            err = (m - p.market_vol) * 1e4
            out.append((p, err if not math.isnan(err) else NAN_PENALTY_BP))
        return out

    def rmse(self, sv: SVParams) -> float:
        e = [x for _, x in self.errors_bp(sv)]
        return float(np.sqrt(np.mean(np.square(e))))


def flat_smile_errors(points: Sequence[SmilePoint]) -> List[Tuple[SmilePoint, float]]:
    """What a Gaussian model (1F or 2F, flat normal vol in the strike) gets wrong: the ATM vol at every strike."""
    return [(p, (p.atm_vol - p.market_vol) * 1e4) for p in points if p.offset_bp != 0.0]


def fit_sv(engine: SmileEngine, start: SVParams, max_fev: int = 120) -> Tuple[SVParams, float, int]:
    x, f, n = nelder_mead(lambda y: engine.rmse(sv_from_x(y)), x_from_sv(start), step=0.7, max_fev=max_fev, xtol=1e-3, ftol=1e-3)
    return sv_from_x(x), f, n


# ---------------------------------------------------------------- caps
def caps_vols(base: HJM2FModel, sv: SVParams, maturities: Sequence[float], n_paths: int, seed: int, t_star: float, delta: float = 0.5) -> List[float]:
    """Model ATM cap vols with the control variate (exact 2F cap price + MC difference)."""
    curve = base.curve
    strikes = {m: atm_cap_strike(curve, m, delta) for m in maturities}
    fixings = sorted({t for m in maturities for t in cap_fixings(m, delta)})
    sums = {s: np.zeros(len(maturities)) for s in ("sv", "zero")}
    for tag, par in (("sv", sv), ("zero", ZERO)):
        acc = np.zeros(len(maturities))

        def observer(st: State, acc=acc) -> None:
            tt = round(st.t, 10)
            if tt not in fx:
                return
            pn = bond(base, st.t, t_star, st)
            p = bond(base, st.t, st.t + delta, st)
            for j, m in enumerate(maturities):
                if tt in fixset[m]:
                    k = 1.0 / (1.0 + delta * strikes[m])
                    acc[j] += (1.0 + delta * strikes[m]) * curve.discount(t_star) * float(np.mean(np.maximum(k - p, 0.0) / pn))

        fx = set(round(f, 10) for f in fixings)
        fixset = {m: set(round(f, 10) for f in cap_fixings(m, delta)) for m in maturities}
        simulate(base, par, t_star, fixings, n_paths, seed, observer)
        sums[tag] = acc
    out = []
    for j, m in enumerate(maturities):
        price = model_cap_price(base, m, strikes[m], delta) + (sums["sv"][j] - sums["zero"][j])
        out.append(implied_cap_vol(curve, m, strikes[m], price, delta) if price > 0 else 0.0)
    return out


# ---------------------------------------------------------------- Bermudan
def bermudan_sv(base: HJM2FModel, sv: SVParams, berm: BermudanCoterminalSwaption, exact_price: float, n_paths: int = 100000, seed: int = 5,
                use_v: bool = True) -> Tuple[float, float, float]:
    """(price with control variate, LSMC_SV, LSMC_0). use_v: add the volatility multiplier to the regression basis of the SV run."""
    out = []
    for par in (sv, ZERO):
        ex = list(berm.exercise_times)
        tn = berm.terminal
        pay: Dict[float, np.ndarray] = {}
        reg: Dict[float, np.ndarray] = {}

        def observer(st: State) -> None:
            pn = bond(base, st.t, tn, st)
            und = berm.underlying_at(st.t)
            pay[round(st.t, 10)] = np.maximum(swap_value(base, und, st, pn), 0.0) / pn
            reg[round(st.t, 10)] = np.column_stack([st.x[:, 0], st.x[:, 1], st.v])

        simulate(base, par, tn, ex, n_paths, seed, observer)
        keys = [round(t, 10) for t in ex]
        value = pay[keys[-1]].copy()
        for k in reversed(keys[:-1]):
            itm = pay[k] > 1e-14
            if not np.any(itm):
                continue
            r = reg[k][itm]
            s = r.std(axis=0) + 1e-12
            u, v, w = r[:, 0] / s[0], r[:, 1] / s[1], (r[:, 2] - 1.0) / s[2] if s[2] > 1e-9 else np.zeros(len(r))
            cols = [np.ones_like(u), u, v, u * u, u * v, v * v, u ** 3, u * u * v, u * v * v, v ** 3]
            if use_v:
                cols += [w, w * u, w * v, w * w]
            basis = np.column_stack(cols)
            beta, *_ = np.linalg.lstsq(basis, value[itm], rcond=1e-10)
            idx = np.where(itm)[0][pay[k][itm] > basis @ beta]
            value[idx] = pay[k][idx]
        out.append(base.curve.discount(tn) * float(np.mean(value)))
    return exact_price + (out[0] - out[1]), out[0], out[1]
