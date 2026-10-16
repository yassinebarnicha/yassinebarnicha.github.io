/* lookback.js — JavaScript port of the project's unified DLL (code/unified/unified_dll.cpp), extended.

   priceLookback(method, ttm, optionType, spot, rate, vol, paths, steps, seed, withStdErrors, useBGK, options)
   returns the 48 values of the C++ function price_lookback: for the Monte Carlo then the closed form, 6 quantities
   (price, Δ, Γ, Θ, ρ, vega) × (estimate, standard error, 95 % CI low, 95 % CI high).
     method 1: Barnicha engine (exact one-step Brownian bridge)   — code/Barnicha/src/lookback_dll.cpp
     method 0: Chivet engine (time stepping, BGK correction)      — code/Chivet/src/LookbackDLL.cpp
   Without "options", the computation is the C++ one (see the differences below).

   Generator: std::mt19937 and libstdc++'s uniform_real_distribution / normal_distribution reproduced exactly: with
   the same seed, the Barnicha engine gives the numbers of the C++ compiled with g++ (checked to 1e-15).

   Deliberate differences with the C++ (defects found in the review, reported on the page):
   - Chivet: the seed entered is used (the C++ takes the clock) and the normal law restarts at each price (the C++
     keeps the reserved normal from one price to the next: the common draws of the Greeks are shifted);
   - Chivet: the standard error of rho is divided by 2·bump·r like the estimate (the C++ forgets the factor r);
   - Chivet, closed form: Γ divided by h² (the C++ divides by (2h)²: Γ four times too small).

   Extensions (options): continuous dividend q; extremum already observed (seasoned option); choice of the discrete
   monitoring correction (none, BGK, or exact Brownian bridge on each step); bump sizes; antithetic variable;
   centred Θ; BGK coefficient β.
   Note: at inception (no extremum given) a spot bump also moves the starting extremum, so the Greeks are those of the
   inception price V(S, S) (Γ = 0); with an extremum given, the extremum stays fixed (Greeks of a held contract). */
(function (racine) {
  "use strict";

  // ------------------------------------------------------------------ std::mt19937
  function MT19937(graine) {
    const mt = new Uint32Array(624);
    let idx = 624;
    mt[0] = graine >>> 0;
    for (let i = 1; i < 624; i++) {
      const p = mt[i - 1] ^ (mt[i - 1] >>> 30);
      mt[i] = (((((p & 0xffff0000) >>> 16) * 1812433253) << 16) + (p & 0x0000ffff) * 1812433253 + i) >>> 0;
    }
    function regenerer() {
      for (let i = 0; i < 624; i++) {
        const y = (mt[i] & 0x80000000) | (mt[(i + 1) % 624] & 0x7fffffff);
        let v = mt[(i + 397) % 624] ^ (y >>> 1);
        if (y & 1) v ^= 0x9908b0df;
        mt[i] = v >>> 0;
      }
      idx = 0;
    }
    return function () {
      if (idx >= 624) regenerer();
      let y = mt[idx++];
      y ^= y >>> 11; y ^= (y << 7) & 0x9d2c5680; y ^= (y << 15) & 0xefc60000; y ^= y >>> 18;
      return y >>> 0;
    };
  }
  // libstdc++'s std::generate_canonical<double, 53>(mt19937): two 32-bit draws
  const DEUX32 = 4294967296, UN_MOINS = 1 - Math.pow(2, -53);
  function canonique(g) {
    const r = (g() + g() * DEUX32) / (DEUX32 * DEUX32);
    return r >= 1 ? UN_MOINS : r;
  }
  // libstdc++'s std::normal_distribution<double>: polar method, one value kept in reserve
  function LoiNormale() {
    let reserve = 0, dispo = false;
    return function (g) {
      if (dispo) { dispo = false; return reserve; }
      let x, y, r2;
      do { x = 2 * canonique(g) - 1; y = 2 * canonique(g) - 1; r2 = x * x + y * y; } while (r2 > 1 || r2 === 0);
      const mult = Math.sqrt(-2 * Math.log(r2) / r2);
      reserve = x * mult; dispo = true;
      return y * mult;
    };
  }

  // ------------------------------------------------------------------ normal law
  function erfc(x) { // ~1e-15: series (|x| < 3) or Lentz continued fraction
    if (x < 0) return 2 - erfc(-x);
    if (x < 3) {
      let t = x, s = x;
      for (let n = 1; n < 200; n++) { t *= 2 * x * x / (2 * n + 1); s += t; if (t < s * 1e-17) break; }
      return 1 - 2 / Math.sqrt(Math.PI) * Math.exp(-x * x) * s;
    }
    let f = x, C = x, D = 0;
    for (let n = 1; n < 300; n++) { const a = n / 2; D = x + a * D; D = 1 / (Math.abs(D) < 1e-300 ? 1e-300 : D); C = x + a / C; if (Math.abs(C) < 1e-300) C = 1e-300; const d = C * D; f *= d; if (Math.abs(d - 1) < 1e-16) break; }
    return Math.exp(-x * x) / Math.sqrt(Math.PI) / f;
  }
  const N = x => 0.5 * erfc(-x / Math.SQRT2);
  const n = x => 0.3989422804014327 * Math.exp(-0.5 * x * x);

  // ------------------------------------------------------------------ closed forms
  // Barnicha (analytic_lookback.cpp): inception (extremum = spot), no dividend, exact limit as r → 0
  function prixFermeBarnicha(put, S, r, sig, T) {
    if (Math.abs(r) < 1e-10) {
      const d = 0.5 * sig * Math.sqrt(T), pd = n(d), cd = N(d);
      return put ? S * ((1 + d * d) * (2 * cd - 1) + 2 * d * pd + d * d) : S * ((1 + d * d) * (2 * cd - 1) + 2 * d * pd - d * d);
    }
    const sq = Math.sqrt(T), a1 = ((r + 0.5 * sig * sig) * T) / (sig * sq), a2 = a1 - sig * sq, a3 = a1 - 2 * r * sq / sig;
    const disc = Math.exp(-r * T), mult = S * sig * sig / (2 * r);
    if (!put) return S * N(a1) - S * disc * N(a2) - mult * (N(-a1) - disc * N(-a3));
    return -S * N(-a1) + S * disc * N(-a2) + mult * (N(a1) - disc * N(a3));
  }
  // Chivet (AnalyticLookback.cpp): the same price written differently; undefined for r = 0 (division by 2r)
  function prixFermeChivet(put, S0, T, sig, r) {
    const dp = (Math.log(1) + (r + 0.5 * sig * sig) * T) / (sig * Math.sqrt(T)), dm = (Math.log(1) + (r - 0.5 * sig * sig) * T) / (sig * Math.sqrt(T));
    if (!put) return S0 * N(dp) - S0 * Math.exp(-r * T) * N(dm) + Math.exp(-r * T) * (sig * sig * S0 / (2 * r)) * N(dm) - N(-dp) * S0 * (sig * sig / (2 * r));
    return S0 * Math.exp(-r * T) * N(-dm) + S0 * N(dp) * (1 + (sig * sig) / (2 * r)) - Math.exp(-r * T) * S0 * (sig * sig / (2 * r)) * N(-dm) - S0;
  }
  // general Goldman-Sosin-Gatto: carry b = r − q, extremum already observed m (min for the call, max for the put)
  function gsg(put, S, m, r, b, sig, T) {
    const sq = Math.sqrt(T), v2 = sig * sig;
    if (!put) {
      const a1 = (Math.log(S / m) + (b + 0.5 * v2) * T) / (sig * sq), a2 = a1 - sig * sq;
      return S * Math.exp((b - r) * T) * N(a1) - m * Math.exp(-r * T) * N(a2)
        + S * Math.exp(-r * T) * v2 / (2 * b) * (Math.pow(S / m, -2 * b / v2) * N(-a1 + 2 * b * sq / sig) - Math.exp(b * T) * N(-a1));
    }
    const b1 = (Math.log(S / m) + (b + 0.5 * v2) * T) / (sig * sq), b2 = b1 - sig * sq;
    return -S * Math.exp((b - r) * T) * N(-b1) + m * Math.exp(-r * T) * N(-b2)
      + S * Math.exp(-r * T) * v2 / (2 * b) * (-Math.pow(S / m, -2 * b / v2) * N(b1 - 2 * b * sq / sig) + Math.exp(b * T) * N(b1));
  }
  function prixFermeGeneral(put, S, m, r, q, sig, T) {
    const b = r - q;
    if (Math.abs(b) < 1e-6) return 0.5 * (gsg(put, S, m, r, 1e-6, sig, T) + gsg(put, S, m, r, -1e-6, sig, T)); // limit b → 0
    return gsg(put, S, m, r, b, sig, T);
  }
  // closed form used by an engine: its own code's formula while within its scope (inception, q = 0)
  function prixFerme(moteur, m) {
    const ext = m.ext == null ? m.S : m.ext;
    if (m.q === 0 && ext === m.S) return moteur === 1 ? prixFermeBarnicha(m.put, m.S, m.r, m.sig, m.T) : prixFermeChivet(m.put, m.S, m.T, m.sig, m.r);
    return prixFermeGeneral(m.put, m.S, ext, m.r, m.q, m.sig, m.T);
  }

  // ------------------------------------------------------------------ statistics (pricing_engine.cpp)
  function stats(xs) {
    const k = xs.length;
    if (!k) return [0, 0, 0, 0];
    let s = 0; for (let i = 0; i < k; i++) s += xs[i];
    const moy = s / k; let v = 0; for (let i = 0; i < k; i++) { const d = xs[i] - moy; v += d * d; }
    const se = k > 1 ? Math.sqrt(v / (k - 1)) / Math.sqrt(k) : 0;
    return [moy, se, moy - 1.96 * se, moy + 1.96 * se];
  }
  const avec = (m, o) => Object.assign({}, m, o);

  // ================================================================== Barnicha engine
  // draws (z, u_max, u_min): OneStepBridgeSampler::generate, in the C++ order
  function tiragesPont(nChemins, graine) {
    const g = MT19937(graine), loi = LoiNormale();
    const z = new Float64Array(nChemins), umax = new Float64Array(nChemins), umin = new Float64Array(nChemins);
    for (let i = 0; i < nChemins; i++) {
      z[i] = loi(g);
      let u = canonique(g); umax[i] = u === 0 ? 1e-16 : u;
      u = canonique(g); umin[i] = u === 0 ? 1e-16 : u;
    }
    // logarithms of the uniforms computed once (reused for every bumped market of the Greeks)
    const lmax = umax.map(Math.log), lmin = umin.map(Math.log);
    return {z, umax, umin, lmax, lmin, n: nChemins};
  }
  // discounted payoff per path (monte_carlo_engine.cpp + lookback_payoff_evaluator.cpp)
  function gainsPont(t, m, antithetique) {
    const out = new Float64Array(t.n), x0 = Math.log(m.S), derive = (m.r - m.q - 0.5 * m.sig * m.sig) * m.T, vol = m.sig * Math.sqrt(m.T);
    const disc = Math.exp(-m.r * m.T), k = -0.5 * m.sig * m.sig * m.T, put = m.put;
    const ext = m.ext == null ? null : m.ext; // extremum already observed (otherwise: the spot)
    const gain = (z, lumax, lumin) => {
      const xT = x0 + derive + vol * z, dx = (x0 - xT) * (x0 - xT);
      if (!put) { let mn = Math.exp(0.5 * (x0 + xT - Math.sqrt(dx + 4 * (k * lumin)))); if (ext !== null && ext < mn) mn = ext; return Math.exp(xT) - mn; }
      let mx = Math.exp(0.5 * (x0 + xT + Math.sqrt(dx + 4 * (k * lumax)))); if (ext !== null && ext > mx) mx = ext; return mx - Math.exp(xT);
    };
    for (let i = 0; i < t.n; i++) {
      const p = gain(t.z[i], t.lmax[i], t.lmin[i]);
      out[i] = disc * (antithetique ? 0.5 * (p + gain(-t.z[i], t.lmax[i], t.lmin[i])) : p);
    }
    return out;
  }
  // path-by-path detail (no antithetic) for the descriptive charts:
  // log(S_T/S), log(bridge extremum / S) and discounted payoff
  function detailPont(t, m) {
    const k = t.n, lx = new Float64Array(k), le = new Float64Array(k), gain = new Float64Array(k);
    const derive = (m.r - m.q - 0.5 * m.sig * m.sig) * m.T, vol = m.sig * Math.sqrt(m.T), c = -0.5 * m.sig * m.sig * m.T, disc = Math.exp(-m.r * m.T);
    for (let i = 0; i < k; i++) {
      const xT = derive + vol * t.z[i], d = xT * xT;
      const e = m.put ? 0.5 * (xT + Math.sqrt(d + 4 * (c * t.lmax[i]))) : 0.5 * (xT - Math.sqrt(d + 4 * (c * t.lmin[i])));
      lx[i] = xT; le[i] = e;
      let ext = m.S * Math.exp(e);
      if (m.ext != null) ext = m.put ? Math.max(ext, m.ext) : Math.min(ext, m.ext);
      gain[i] = disc * (m.put ? ext - m.S * Math.exp(xT) : m.S * Math.exp(xT) - ext);
    }
    return {lx, le, gain};
  }
  // exact distribution function of the extremum of X_t = μt + σW_t on [0, T] (in log-return)
  //   max: P(M ≤ y) = N((y − μT)/(σ√T)) − e^{2μy/σ²} N((−y − μT)/(σ√T)), y ≥ 0
  //   min: P(m ≥ y) = N((−y + μT)/(σ√T)) − e^{2μy/σ²} N((y + μT)/(σ√T)), y ≤ 0
  function repartitionExtremum(put, y, mu, sig, T) {
    const s = sig * Math.sqrt(T);
    if (put) return y <= 0 ? 0 : N((y - mu * T) / s) - Math.exp(2 * mu * y / (sig * sig)) * N((-y - mu * T) / s);
    return y >= 0 ? 0 : N((-y + mu * T) / s) - Math.exp(2 * mu * y / (sig * sig)) * N((y + mu * T) / s);
  }
  // BumpAndRevalueGreeks::compute (default bumps of bump_and_revalue_greeks.h: 1e-4)
  function grecquesBarnicha(gains, m, bp, thetaCentre) {
    const base = gains(m), k = base.length, o = new Array(24);
    const fixer = (j, arr) => { const s = stats(arr); for (let q = 0; q < 4; q++) o[j * 4 + q] = s[q]; };
    const tab = f => { const a = new Float64Array(k); for (let i = 0; i < k; i++) a[i] = f(i); return a; };
    fixer(0, base);
    const hS = Math.max(1e-12, m.S * bp.S), up = gains(avec(m, {S: m.S + hS})), dn = gains(avec(m, {S: m.S - hS}));
    fixer(1, tab(i => (up[i] - dn[i]) / (2 * hS)));
    fixer(2, tab(i => (up[i] - 2 * base[i] + dn[i]) / (hS * hS)));
    const hT = Math.max(1e-12, bp.T), Tdn = Math.max(1e-10, m.T - hT), tD = gains(avec(m, {T: Tdn}));
    if (thetaCentre) { const tU = gains(avec(m, {T: m.T + hT})); fixer(3, tab(i => -(tU[i] - tD[i]) / (m.T + hT - Tdn))); }
    else fixer(3, tab(i => (tD[i] - base[i]) / hT));
    const hR = Math.max(1e-12, bp.r), rU = gains(avec(m, {r: m.r + hR})), rD = gains(avec(m, {r: m.r - hR}));
    fixer(4, tab(i => (rU[i] - rD[i]) / (2 * hR)));
    const hV = Math.max(1e-12, bp.sig), vU = gains(avec(m, {sig: m.sig + hV})), vD = gains(avec(m, {sig: m.sig - hV}));
    fixer(5, tab(i => (vU[i] - vD[i]) / (2 * hV)));
    return o;
  }
  function barnicha(m, paths, seed, o) {
    const t = tiragesPont(paths, seed), bp = Object.assign({S: 1e-4, r: 1e-4, sig: 1e-4, T: 1e-4}, o.bumpsBarnicha);
    const mc = grecquesBarnicha(x => gainsPont(t, x, o.antithetique !== false), m, bp, !!o.thetaCentre);
    const an = grecquesBarnicha(x => new Float64Array([prixFerme(1, x)]), m, bp, !!o.thetaCentre);
    return mc.concat(an);
  }

  // ================================================================== Chivet engine
  const BGK_BETA1 = 0.5826; // −ζ(1/2)/√(2π)
  // computeRawPrice: the same draws (fixed seed) for every bumped market.
  // correction: "aucune" (none), "bgk" (C++), or "pont" (bridge): exact Brownian-bridge extremum on each step (unbiased)
  // Several bumped markets evaluated on the SAME draws: each path draws its normals (and its uniforms for the bridge)
  // once, in the C++ order, then replays them for every market. Identical to one simulation per market reseeded with
  // the same seed (what the C++ does), in a single pass of the generator.
  // Computed in logs: the extremum of S is the exponential of the extremum of log S.
  function gainsChivetMulti(seed, paths, steps, correction, beta, marches) {
    const g = MT19937(seed), loi = LoiNormale(), pont = correction === "pont", bgk = correction === "bgk";
    const z = new Float64Array(steps), lu = pont ? new Float64Array(steps) : null;
    const P = marches.map(m => {
      const dt = m.T / steps, s = m.sig * Math.sqrt(dt);
      return {put: m.put, x0: Math.log(m.S), derive: (m.r - m.q - 0.5 * m.sig * m.sig) * dt, s, v2dt: 2 * s * s, disc: Math.exp(-m.r * m.T),
        lcorr: bgk ? beta * m.sig * Math.sqrt(dt) : 0, ext: m.ext == null ? null : m.ext, out: new Float64Array(paths)};
    });
    for (let i = 0; i < paths; i++) {
      for (let k = 0; k < steps; k++) { z[k] = loi(g); if (pont) lu[k] = Math.log(canonique(g) || 1e-16); }
      for (const q of P) {
        let x = q.x0, xmax = x, xmin = x;
        const derive = q.derive, s = q.s;
        if (pont) {
          // extremum of a Brownian bridge between x and x1: P(max ≤ y) = 1 − exp(−2 (y − x)(y − x1) / (σ² Δt))
          for (let k = 0; k < steps; k++) {
            const x1 = x + derive + s * z[k], dx = x1 - x, w = Math.sqrt(dx * dx - q.v2dt * lu[k]);
            if (q.put) { const e = 0.5 * (x + x1 + w); if (e > xmax) xmax = e; } else { const e = 0.5 * (x + x1 - w); if (e < xmin) xmin = e; }
            x = x1;
          }
        } else {
          for (let k = 0; k < steps; k++) { x += derive + s * z[k]; if (x > xmax) xmax = x; else if (x < xmin) xmin = x; }
        }
        let smax = Math.exp(xmax + q.lcorr), smin = Math.exp(xmin - q.lcorr);
        if (q.ext !== null) { if (q.ext < smin) smin = q.ext; if (q.ext > smax) smax = q.ext; }
        const ST = Math.exp(x);
        q.out[i] = q.disc * (q.put ? smax - ST : ST - smin);
      }
    }
    return P.map(q => q.out);
  }
  const gainsChivet = (seed, paths, steps, correction, beta, m) => gainsChivetMulti(seed, paths, steps, correction, beta, [m])[0];
  function chivet(m, paths, steps, seed, stdErr, correction, o) {
    const bc = Object.assign({S: 0.001, sig: 0.001, r: 0.0001, T: 1 / 365}, o.bumpsChivet), beta = o.beta == null ? BGK_BETA1 : o.beta;
    const moy = a => { let s = 0; for (let i = 0; i < a.length; i++) s += a[i]; return s / a.length; };
    const ecart = a => stats(a)[1];
    // the 9 markets of Pricer::compute, on the same draws
    const [base, up, dn, tA, tB, rU, rD, vU, vD] = gainsChivetMulti(seed, paths, steps, correction, beta, [
      m, avec(m, {S: m.S * (1 + bc.S)}), avec(m, {S: m.S * (1 - bc.S)}), avec(m, {T: m.T * (1 + bc.T)}), avec(m, {T: m.T * (1 - bc.T)}),
      avec(m, {r: m.r * (1 + bc.r)}), avec(m, {r: m.r * (1 - bc.r)}), avec(m, {sig: m.sig * (1 + bc.sig)}), avec(m, {sig: m.sig * (1 - bc.sig)})]);
    const k = base.length, tab = f => { const a = new Float64Array(k); for (let i = 0; i < k; i++) a[i] = f(i); return a; };
    const hS = m.S * bc.S, hV = m.sig * bc.sig;
    const p0 = moy(base), pU = moy(up), pD = moy(dn), dT = 2 * m.T * bc.T, dR = 2 * bc.r * m.r;
    const est = [p0, (pU - pD) / (2 * hS), (pU - 2 * p0 + pD) / (hS * hS), (moy(tA) - moy(tB)) / dT, (moy(rU) - moy(rD)) / dR, (moy(vU) - moy(vD)) / (2 * hV)];
    const se = [ecart(base), 0, 0, 0, 0, 0];
    if (stdErr) {
      se[1] = ecart(tab(i => (up[i] - dn[i]) / (2 * hS)));
      se[2] = ecart(tab(i => (up[i] - 2 * base[i] + dn[i]) / (hS * hS)));
      se[3] = ecart(tab(i => (tA[i] - tB[i]) / dT));
      se[4] = ecart(tab(i => (rU[i] - rD[i]) / dR));
      se[5] = ecart(tab(i => (vU[i] - vD[i]) / (2 * hV)));
    }
    const out = [];
    for (let j = 0; j < 6; j++) out.push(est[j], se[j], est[j] - 1.96 * se[j], est[j] + 1.96 * se[j]);
    // closed form and Greeks by relative bumps of 1e-4 (AnalyticLookback.cpp), Γ corrected
    const F = x => prixFerme(0, avec(m, x)), V0 = F({}), Vp = F({S: m.S * 1.0001}), Vm = F({S: m.S * 0.9999});
    const ex = v => [v, 0, v, v];
    out.push(...ex(V0), ...ex((Vp - Vm) / (m.S * 0.0002)), ...ex((Vp - 2 * V0 + Vm) / Math.pow(m.S * 0.0001, 2)),
      ...ex((F({T: m.T * 1.0001}) - F({T: m.T * 0.9999})) / (m.T * 0.0002)),
      ...ex((F({r: m.r * 1.0001}) - F({r: m.r * 0.9999})) / (m.r * 0.0002)),
      ...ex((F({sig: m.sig * 1.0001}) - F({sig: m.sig * 0.9999})) / (m.sig * 0.0002)));
    return out;
  }

  // ================================================================== entry point (unified_dll.cpp)
  function priceLookback(method, ttm, optionType, spot, rate, vol, paths, steps, seed, withStdErrors, useBGK, options) {
    const o = options || {};
    if (!(ttm > 0)) throw new Error("zero or negative maturity (DLL code -2)");
    if (!(paths > 0)) throw new Error("zero number of paths (DLL code -3)");
    const put = optionType === 1, ext = o.extremum == null ? null : +o.extremum;
    if (ext !== null && (put ? ext < spot : ext > spot)) throw new Error(put ? "the maximum already observed must be ≥ the spot" : "the minimum already observed must be ≤ the spot");
    const m = {put, S: spot, r: rate, q: o.q || 0, sig: vol, T: ttm, ext};
    if (method === 1) return barnicha(m, Math.round(paths), seed >>> 0, o);
    if (method === 0) return chivet(m, Math.round(paths), Math.max(1, Math.round(steps)), seed >>> 0, !!withStdErrors, o.correction || (useBGK ? "bgk" : "aucune"), o);
    throw new Error("unknown method (DLL code -999)");
  }

  const API = {priceLookback, MT19937, canonique, LoiNormale, tiragesPont, gainsPont, gainsChivet, prixFermeBarnicha, prixFermeChivet,
    prixFermeGeneral, detailPont, repartitionExtremum, gainsChivetMulti, BGK_BETA1};
  if (typeof module !== "undefined" && module.exports) module.exports = API; else racine.LB = API;
})(typeof self !== "undefined" ? self : this);
