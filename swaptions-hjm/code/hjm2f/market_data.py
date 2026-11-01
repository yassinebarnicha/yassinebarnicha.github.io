from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from hjm1f.calibration_1f import CalibrationInstrument
from hjm1f.curve import ZeroCurve
from hjm1f.instruments import CoterminalSwaption
from hjm1f.market_quotes import forward_swap_rate

from .calibration_2f import atm_swaption, implied_normal_vol
from .calibration_global import Targets, model_errors_bp
from .caps import model_cap_vol
from .hjm_2f import HJM2FModel

DEFAULT_PATH = Path(__file__).resolve().parent.parent / "data" / "eur_2016-03-30.json"
MAX_HORIZON = 30.0  # the curve stops at 30 years: a swaption must end before


@dataclass
class SmileSurface:
    """
    Smile of normal volatility, interpolated between the 13 expiry x tenor smiles of the data.

    For each tabulated smile, vol(offset) / vol(ATM) is a piecewise-linear function of the strike offset (bp, strike minus
    forward swap rate; flat beyond the tabulated range). The ratio at a given (expiry, tenor, offset) is the average of these
    ratios weighted by a Gaussian kernel in (log expiry, log tenor).
    """

    offsets: List[float]
    pairs: List[Tuple[float, float, List[Optional[float]], float]]  # expiry, tenor, vols, ATM vol
    bandwidth: float = 0.8

    def _ratio_one(self, vols: Sequence[Optional[float]], atm: float, offset_bp: float) -> float:
        pts = [(o, v / atm) for o, v in zip(self.offsets, vols) if v is not None]
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        return float(np.interp(offset_bp, xs, ys))

    def ratio(self, expiry: float, tenor: float, offset_bp: float) -> float:
        num = den = 0.0
        for e, t, vols, atm in self.pairs:
            d2 = (math.log(expiry / e) ** 2 + math.log(tenor / t) ** 2) / (self.bandwidth ** 2)
            w = math.exp(-0.5 * d2)
            num += w * self._ratio_one(vols, atm, offset_bp)
            den += w
        return num / den if den > 1e-300 else 1.0


class MarketData:
    """The real market of the course data: EUR curve, ATM swaption matrix, ATM caps, smiles (30 March 2016)."""

    def __init__(self, path: Path = DEFAULT_PATH):
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        self.raw = d
        self.curve = ZeroCurve([c["t"] for c in d["curve"]], [c["df"] for c in d["curve"]])
        m = d["swaptionNormalVolATM"]
        self.expiries: List[float] = m["expiries"]
        self.tenors: List[int] = m["tenors"]
        self.matrix = np.array(m["vol"], dtype=float)
        self.cap_maturities: List[float] = d["capsATM"]["maturities"]
        self.cap_vols: List[float] = d["capsATM"]["atmVol"]
        sm = d["swaptionSmile"]
        pairs = []
        for p in sm["pairs"]:
            atm = self.swaption_vol(p["expiry"], int(p["tenor"]))
            pairs.append((p["expiry"], p["tenor"], p["vol"], atm))
        self.smile = SmileSurface([float(o) for o in sm["offsetsBp"]], pairs)

    # ---- quotes
    def swaption_vol(self, expiry: float, tenor: int) -> float:
        return float(self.matrix[self.expiries.index(float(expiry))][self.tenors.index(int(tenor))])

    def has(self, expiry: float, tenor: int) -> bool:
        return float(expiry) in self.expiries and int(tenor) in self.tenors and expiry + tenor <= MAX_HORIZON + 1e-9

    def cap_vol(self, maturity: float) -> float:
        return self.cap_vols[self.cap_maturities.index(float(maturity))]

    # ---- the Bermudan
    def bermudan(self, terminal: int, exercise: Sequence[int]):
        """Strike = par rate of the swap to `terminal`; payment dates annual; returns (strike, payments, accruals)."""
        pay = [float(i) for i in range(1, terminal + 1)]
        acc = [1.0] * terminal
        return forward_swap_rate(self.curve, 0.0, float(terminal), pay, acc), pay, acc

    def coterminal_instruments(self, terminal: int, exercise: Sequence[int], use_smile: bool) -> List[CalibrationInstrument]:
        """Co-terminal swaptions struck at the Bermudan's rate; their quote is the ATM vol, moved along the smile if asked."""
        strike, pay, acc = self.bermudan(terminal, exercise)
        out = []
        for e in exercise:
            p = [x for x in pay if x > e]
            fwd = forward_swap_rate(self.curve, float(e), float(terminal), p, [1.0] * len(p))
            vol = self.swaption_vol(e, terminal - e)
            if use_smile:
                vol *= self.smile.ratio(float(e), float(terminal - e), (strike - fwd) * 1e4)
            out.append(CalibrationInstrument(CoterminalSwaption(float(e), float(terminal), strike, pay, acc, True), vol))
        return out

    def bucket_plan(self, terminal: int, exercise: Sequence[int], use_smile: bool):
        """
        Volatility segments at EVERY expiry of the matrix up to the last exercise date (1M, 2M, 3M, 6M, 9M, 1Y, 18M, 2Y, ...), each fitted
        exactly: on the co-terminal swaption quoted at the Bermudan's strike for an exercise date, on the ATM swaption of tenor
        closest to terminal - expiry for the other expiries. Returns (bucket_ends, instruments).
        """
        window = float(max(exercise))
        ends = [e for e in self.expiries if e <= window + 1e-9]
        coterminal = {i.swaption.expiry: i for i in self.coterminal_instruments(terminal, exercise, use_smile)}
        instruments = []
        for e in ends:
            if e in coterminal:
                instruments.append(coterminal[e])
            else:
                tenor = min(max(int(math.floor(terminal - e + 0.5)), 1), 10)
                instruments.append(CalibrationInstrument(atm_swaption(self.curve, e, float(tenor)), self.swaption_vol(e, tenor)))
        return ends, instruments
    # ---- everything the model is asked to reproduce
    def targets(self, window: float, tenors: Sequence[int], with_caps: bool) -> Targets:
        """Every ATM swaption whose expiry is within the window, and (if asked) every ATM cap whose last fixing is."""
        sw = [(atm_swaption(self.curve, e, float(t)), self.swaption_vol(e, t)) for e in self.expiries if e <= window + 1e-9
              for t in tenors if self.has(e, t)]
        caps = [(m, v) for m, v in zip(self.cap_maturities, self.cap_vols) if with_caps and m - 0.5 <= window + 1e-9]
        return Targets(sw, caps)
    # ---- report card on the whole market
    def report_card(self, model: HJM2FModel, window: float, tenors: Sequence[int], with_caps: bool, n_grid: int = 121) -> Dict:
        """
        Model minus market (bp of normal vol) on EVERY ATM swaption of the matrix that ends within the curve, and on every ATM cap.
        group = 'fitted' (in the calibration set), 'window' (expires within the exercise window but was not fitted: sub-annual
        expiries, long tenors) or 'beyond' (expires after the last exercise date: the last volatility segment is extrapolated).
        """
        rows = []
        for i, e in enumerate(self.expiries):
            for j, t in enumerate(self.tenors):
                if e + t > MAX_HORIZON + 1e-9:
                    continue
                sw = atm_swaption(self.curve, e, float(t))
                price = model.european_swaption_price(sw, n_grid=n_grid)
                mv = implied_normal_vol(self.curve, sw, price) if price > 0 else 0.0
                group = "fitted" if (e <= window + 1e-9 and t in tenors) else ("window" if e <= window + 1e-9 else "beyond")
                rows.append({"expiry": e, "tenor": t, "market": self.matrix[i][j], "model": mv, "err_bp": (mv - self.matrix[i][j]) * 1e4, "group": group})
        caps = []
        for m, v in zip(self.cap_maturities, self.cap_vols):
            mv = model_cap_vol(model, m)
            caps.append({"maturity": m, "market": v, "model": mv, "err_bp": (mv - v) * 1e4, "group": "fitted" if (with_caps and m - 0.5 <= window + 1e-9) else ("window" if m - 0.5 <= window + 1e-9 else "beyond")})
        return {"swaptions": rows, "caps": caps}


def rmse(values: Sequence[float]) -> float:
    return float(np.sqrt(np.mean(np.square(values)))) if len(values) else float("nan")


def summarise(card: Dict) -> Dict[str, float]:
    out = {}
    for kind in ("swaptions", "caps"):
        for g in ("fitted", "window", "beyond"):
            errs = [r["err_bp"] for r in card[kind] if r["group"] == g]
            out[f"{kind}_{g}"] = rmse(errs)
            out[f"n_{kind}_{g}"] = len(errs)
    return out