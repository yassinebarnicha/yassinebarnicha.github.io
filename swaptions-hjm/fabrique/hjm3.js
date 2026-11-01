/* hjm3.js — JavaScript port of the Python package hjm3sv/ (two-factor Gaussian HJM with a stochastic-volatility multiplier),
   built on hjm.js (window.HJM) and hjm2.js (window.HJM2). Runs in the browser (window.HJM3) and under Node for the tests.

   Both forward-rate volatilities are scaled by m_t = exp(y_t - v(t)), V_t = m_t^2:
       dy = -kappa y dt + xi dZ,   v(t) = xi^2 (1 - exp(-2 kappa t)) / (2 kappa),   dZ = rho_v dW^1 + sqrt(1 - rho_v^2) dB
   (E[V_t] = 1 under the spot measure). The no-arbitrage drift holds path by path, so the bond formula of the 2F model survives with
   RANDOM covariances A_kl(t), G_kl(t). Everything is simulated under the T*-forward measure, where X has the drift
       -V_t sigma_k e^{lambda_k t} sum_l rho_kl sigma_l B_l(T* - t)
   and y the extra drift  -xi rho_v m_t (sigma_1 B_1(T* - t) + rho sigma_2 B_2(T* - t)).  With xi = 0 the simulation is the exact
   2F transition, so a run with xi = 0 on the same random numbers is a control variate for the exact 2F price:
       price_SV = price_2F(exact) + (MC_SV - MC_{xi=0}).                                                                       */
(function (racine) {
  "use strict";
  const node = typeof module !== "undefined" && module.exports;
  const H = node ? require("./hjm.js") : racine.HJM, H2 = node ? require("./hjm2.js") : racine.HJM2;
  const B = (l, t) => Math.abs(l) < 1e-12 ? t : -Math.expm1(-l * t) / l;
  const ZERO = {kappa: 1, xi: 0, rhoV: 0};
  const KAPPA = [0.02, 3], XI = [0.02, 2.5], RHOV_MAX = 0.95, NAN_PENALTY_BP = 100;

  // time grid: 1/24 up to 5 years, 1/12 beyond, always containing the observation times
  function grille(obs, fineUntil = 5, dtF = 1 / 24, dtC = 1 / 12) {
    const last = Math.max(...obs), pts = new Set(obs.map(t => +t.toFixed(10)));
    let t = 0; pts.add(0);
    while (t < last - 1e-12) { t += t < fineUntil - 1e-12 ? dtF : dtC; pts.add(+Math.min(t, last).toFixed(10)); }
    return Array.from(pts).sort((a, b) => a - b);
  }

  // P(t,T) on every path
  function bondArr(M, t, T, S, out) {
    const lam = M.lam, bt0 = B(lam[0], t), bt1 = B(lam[1], t), bT0 = B(lam[0], T), bT1 = B(lam[1], T), d0 = bT0 - bt0, d1 = bT1 - bt1;
    const c0 = Math.log(M.courbe.discount(T) / M.courbe.discount(t)), k00 = bT0 * bT0 - bt0 * bt0, k01 = 2 * (bT0 * bT1 - bt0 * bt1), k11 = bT1 * bT1 - bt1 * bt1;
    const n = S.n;
    for (let j = 0; j < n; j++) {
      const det = -0.5 * (k00 * S.a11[j] + k01 * S.a12[j] + k11 * S.a22[j] - 2 * (d0 * S.g11[j] + d0 * S.g12[j] + d1 * S.g21[j] + d1 * S.g22[j]));
      out[j] = Math.exp(c0 + det - d0 * S.x1[j] - d1 * S.x2[j]);
    }
    return out;
  }

  // Simulate (X, y, A, G) under the T*-forward measure; call observer(S) at each observation time (S.t).
  function simuler(M, sv, tStar, obs, nPaths, graine, observer) {
    const g = grille(obs), obsSet = new Set(obs.map(t => +t.toFixed(10))), rng = H.generateur(graine), n = nPaths - (nPaths % 2), rho = M.rho;
    const S = {n, t: 0, x1: new Float64Array(n), x2: new Float64Array(n), y: new Float64Array(n), v: new Float64Array(n).fill(1),
      a11: new Float64Array(n), a12: new Float64Array(n), a22: new Float64Array(n), g11: new Float64Array(n), g12: new Float64Array(n), g21: new Float64Array(n), g22: new Float64Array(n)};
    const sV = Math.sqrt(Math.max(1 - sv.rhoV * sv.rhoV, 0)), half = n / 2;
    const u1 = new Float64Array(n), u2 = new Float64Array(n), u3 = new Float64Array(n), vNew = new Float64Array(n), yNew = new Float64Array(n);
    let aP = [[0, 0], [0, 0]], gP = [[0, 0], [0, 0]], mP = [0, 0];
    for (let i = 1; i < g.length; i++) {
      const t0 = g[i - 1], t1 = g[i], dt = t1 - t0, {A, G} = M.AG(t1), m1 = M.moyenneForward(t1, tStar);
      const da00 = A[0][0] - aP[0][0], da01 = A[0][1] - aP[0][1], da11 = A[1][1] - aP[1][1];
      const dg00 = G[0][0] - gP[0][0], dg01 = G[0][1] - gP[0][1], dg10 = G[1][0] - gP[1][0], dg11 = G[1][1] - gP[1][1];
      const dm0 = m1[0] - mP[0], dm1 = m1[1] - mP[1];
      aP = A; gP = G; mP = m1;
      const l11 = Math.sqrt(Math.max(da00, 0)), l21 = l11 > 1e-150 ? da01 / l11 : 0, l22 = Math.sqrt(Math.max(da11 - l21 * l21, 0));
      for (let j = 0; j < half; j++) {
        const a = rng.normale(), b = rng.normale(), c = rng.normale();
        u1[j] = a; u2[j] = b; u3[j] = c; u1[j + half] = -a; u2[j + half] = -b; u3[j + half] = -c;
      }
      let c = 0, decay = 1, sd = 0, vt = 0;
      const stoch = sv.xi !== 0;
      if (stoch) {
        const tm = 0.5 * (t0 + t1), s = M.sigmaA(tm);
        c = s[0] * B(M.lam[0], tStar - tm) + rho * s[1] * B(M.lam[1], tStar - tm);
        decay = Math.exp(-sv.kappa * dt); sd = sv.xi * Math.sqrt((1 - Math.exp(-2 * sv.kappa * dt)) / (2 * sv.kappa));
        vt = sv.xi * sv.xi * (1 - Math.exp(-2 * sv.kappa * t1)) / (2 * sv.kappa);
      }
      for (let j = 0; j < n; j++) {
        const v = S.v[j], sq = Math.sqrt(v);
        S.x1[j] += dm0 * v + sq * l11 * u1[j];
        S.x2[j] += dm1 * v + sq * (l21 * u1[j] + l22 * u2[j]);
        let vn = v;
        if (stoch) {
          const zv = sv.rhoV * u1[j] + sV * u3[j];
          let yn = S.y[j] * decay + sd * zv - sv.xi * sv.rhoV * sq * c * dt;
          yn = yn > 4 ? 4 : yn < -4 ? -4 : yn;
          S.y[j] = yn; vn = Math.exp(2 * (yn - vt));
        }
        const vb = 0.5 * (v + vn);
        S.a11[j] += vb * da00; S.a12[j] += vb * da01; S.a22[j] += vb * da11;
        S.g11[j] += vb * dg00; S.g12[j] += vb * dg01; S.g21[j] += vb * dg10; S.g22[j] += vb * dg11;
        S.v[j] = vn;
      }
      S.t = t1;
      if (obsSet.has(+t1.toFixed(10))) observer(S);
    }
  }

  // Monte Carlo prices of European swaptions (deflated payoff, T* measure): one simulation for all
  function europeens(M, sv, swaptions, tStar, nPaths, graine) {
    const out = new Float64Array(swaptions.length), parTemps = new Map();
    swaptions.forEach((s, i) => { const k = +s.expiry.toFixed(10); if (!parTemps.has(k)) parTemps.set(k, []); parTemps.get(k).push(i); });
    let pn = null, pT = null, tmp = null, ann = null;
    simuler(M, sv, tStar, Array.from(parTemps.keys()), nPaths, graine, S => {
      const idx = parTemps.get(+S.t.toFixed(10)); if (!idx) return;
      if (!pn || pn.length !== S.n) { pn = new Float64Array(S.n); pT = new Float64Array(S.n); tmp = new Float64Array(S.n); ann = new Float64Array(S.n); }
      bondArr(M, S.t, tStar, S, pn);
      for (const i of idx) {
        const sw = swaptions[i];
        bondArr(M, S.t, sw.terminal, S, pT); ann.fill(0);
        sw.rest.forEach((T, q) => { bondArr(M, S.t, T, S, tmp); const a = sw.acc[q]; for (let j = 0; j < S.n; j++) ann[j] += a * tmp[j]; });
        let sum = 0;
        for (let j = 0; j < S.n; j++) { const p = 1 - pT[j] - sw.strike * ann[j], val = sw.payer ? p : -p; if (val > 0) sum += val / pn[j]; }
        out[i] = M.courbe.discount(tStar) * sum / S.n;
      }
    });
    return out;
  }

  const volSegura = (courbe, sw, p) => { if (!(p > 0)) return NaN; const v = H.volNImplicite(courbe, sw, p); return v; };

  // implied normal vols of the SV model with the control variate (zero: the xi = 0 prices, computed once and cached by the caller)
  function volsSV(M, sv, swaptions, tStar, nPaths, graine, exact, mc0) {
    exact = exact || swaptions.map(s => M.prixEuropeen(s, 121)); mc0 = mc0 || europeens(M, ZERO, swaptions, tStar, nPaths, graine);
    const mc = europeens(M, sv, swaptions, tStar, nPaths, graine);
    return swaptions.map((s, i) => volSegura(M.courbe, s, exact[i] + mc[i] - mc0[i]));
  }

  // ------------------------------------------------------------------ smiles
  // strikes more than zMax ATM standard deviations away are left out (their prices are below the Monte Carlo resolution)
  function pointsSmile(courbe, D, paires, zMax = 3.2) {
    const pts = [], off = D.smile.offsetsBp;
    for (const p of D.smile.pairs) {
      if (!paires.some(q => q[0] === p.expiry && q[1] === p.tenor) || !H2.aSwaption(D, p.expiry, p.tenor)) continue;
      const e = p.expiry, t = p.tenor, atm = H2.volATM(D, e, t), pay = Array.from({length: t}, (_, i) => e + i + 1), acc = pay.map(() => 1);
      const fwd = H.tauxSwapForward(courbe, H.Swaption(e, pay[t - 1], 0, pay, acc, true));
      pts.push({e, t, off: 0, sw: H.Swaption(e, pay[t - 1], fwd, pay, acc, true), marche: NaN, atm});
      p.vol.forEach((v, k) => {
        if (v == null || Math.abs(off[k]) / (atm * 1e4 * Math.sqrt(e)) > zMax) return;
        pts.push({e, t, off: off[k], sw: H.Swaption(e, pay[t - 1], fwd + off[k] / 1e4, pay, acc, off[k] >= 0), marche: v, atm});
      });
    }
    return pts;
  }
  function erreursFlat(pts) { return pts.filter(p => p.off !== 0).map(p => Object.assign({}, p, {err: (p.atm - p.marche) * 1e4})); }

  function MoteurSmile(M, pts, tStar, nPaths = 40000, graine = 21) {
    const sws = pts.map(p => p.sw), exact = sws.map(s => M.prixEuropeen(s, 121)), mc0 = europeens(M, ZERO, sws, tStar, nPaths, graine);
    function erreurs(sv, formeSeule = true) {
      const v = volsSV(M, sv, sws, tStar, nPaths, graine, exact, mc0), atm = new Map();
      pts.forEach((p, i) => { if (p.off === 0) atm.set(p.e + "|" + p.t, v[i]); });
      const out = [];
      pts.forEach((p, i) => {
        if (p.off === 0) return;
        let m = formeSeule ? p.atm * v[i] / atm.get(p.e + "|" + p.t) : v[i];
        let err = (m - p.marche) * 1e4; if (Number.isNaN(err)) err = NAN_PENALTY_BP;
        out.push(Object.assign({}, p, {modele: m, err}));
      });
      return out;
    }
    return {erreurs, rmse: sv => H2.rmse(erreurs(sv).map(r => r.err))};
  }

  const logistique = (x, lo, hi) => lo + (hi - lo) / (1 + Math.exp(-x));
  const logit = (p, lo, hi) => { const q = Math.min(Math.max((p - lo) / (hi - lo), 1e-6), 1 - 1e-6); return Math.log(q / (1 - q)); };
  const svDepuisX = x => ({kappa: logistique(x[0], KAPPA[0], KAPPA[1]), xi: logistique(x[1], XI[0], XI[1]), rhoV: RHOV_MAX * Math.tanh(x[2])});
  const xDepuisSV = sv => [logit(sv.kappa, KAPPA[0], KAPPA[1]), logit(sv.xi, XI[0], XI[1]), Math.atanh(Math.min(Math.max(sv.rhoV / RHOV_MAX, -0.999), 0.999))];
  function ajusterSV(moteur, depart, maxFev = 100, progres) {
    const r = H2.nelderMead(y => moteur.rmse(svDepuisX(y)), xDepuisSV(depart), 0.7, maxFev, 1e-3, 1e-3, progres);
    return {sv: svDepuisX(r.x), rmse: r.f, evaluations: r.nfev};
  }

  // ------------------------------------------------------------------ caps (exact 2F cap price + MC difference)
  function volsCaps(M, sv, maturites, tStar, nPaths, graine, delta = 0.5) {
    const courbe = M.courbe, strikes = maturites.map(m => H2.strikeATMCap(courbe, m, delta));
    const fixes = maturites.map(m => new Set(H2.capFixings(m, delta).map(t => +t.toFixed(10))));
    const tous = Array.from(new Set([].concat(...maturites.map(m => H2.capFixings(m, delta).map(t => +t.toFixed(10)))))).sort((a, b) => a - b);
    const somme = par => {
      const acc = new Float64Array(maturites.length); let pn = null, p = null;
      simuler(M, par, tStar, tous, nPaths, graine, S => {
        const tt = +S.t.toFixed(10); if (!tous.includes(tt)) return;
        if (!pn || pn.length !== S.n) { pn = new Float64Array(S.n); p = new Float64Array(S.n); }
        bondArr(M, S.t, tStar, S, pn); bondArr(M, S.t, S.t + delta, S, p);
        maturites.forEach((m, k) => {
          if (!fixes[k].has(tt)) return;
          const kb = 1 / (1 + delta * strikes[k]); let s = 0;
          for (let j = 0; j < S.n; j++) { const q = kb - p[j]; if (q > 0) s += q / pn[j]; }
          acc[k] += (1 + delta * strikes[k]) * courbe.discount(tStar) * s / S.n;
        });
      });
      return acc;
    };
    const a = somme(sv), z = somme(ZERO);
    return maturites.map((m, k) => {
      const prix = H2.prixCapModele(M, m, strikes[k], delta) + (a[k] - z[k]);
      return prix > 0 ? H2.volCapImplicite(courbe, m, strikes[k], prix, delta) : 0;
    });
  }

  // ------------------------------------------------------------------ re-bootstrap of the volatility level under the SV dynamics
  // proportional fixed-point iteration: sigma_j <- sigma_j * quote_j / model_j (the price is almost proportional to the level)
  function rebootstrap(courbe, lam, rho, ratio, fins, instruments, sv, sigma0, nPaths = 60000, iterations = 14, progres) {
    let s = sigma0.slice(), pire = NaN, M;
    const sws = instruments.map(i => i.swaption), cible = instruments.map(i => i.volN);
    for (let it = 0; it < iterations; it++) {
      M = H2.Modele2(courbe, lam, fins, s, s.map(v => ratio * v), rho);
      const v = volsSV(M, sv, sws, 10, nPaths, 17), err = v.map((x, i) => (x - cible[i]) * 1e4);
      pire = Math.max(...err.map(Math.abs)); if (progres) progres((it + 1) / iterations);
      if (!(pire > 0.05)) break;
      s = s.map((x, i) => x * Math.pow(cible[i] / v[i], 1.6)); // over-relaxed: the segments are coupled
    }
    return {sigmas: s, modele: M, erreurMax: pire};
  }

  // ------------------------------------------------------------------ Bermudan: Longstaff-Schwartz with the exact 2F price as control variate
  function lsmcSV(M, sv, berm, nPaths, graine) {
    const ex = berm.exerciseTimes, tn = berm.terminal, pay = [], reg = [];
    let pn = null, pT = null, tmp = null, ann = null;
    simuler(M, sv, tn, ex, nPaths, graine, S => {
      if (!pn || pn.length !== S.n) { pn = new Float64Array(S.n); pT = new Float64Array(S.n); tmp = new Float64Array(S.n); ann = new Float64Array(S.n); }
      bondArr(M, S.t, tn, S, pn);
      const sw = berm.sousJacente(S.t); ann.fill(0);
      sw.rest.forEach((T, q) => { bondArr(M, S.t, T, S, tmp); const a = sw.acc[q]; for (let j = 0; j < S.n; j++) ann[j] += a * tmp[j]; });
      const p = new Float64Array(S.n);
      for (let j = 0; j < S.n; j++) { const v = 1 - pn[j] - sw.strike * ann[j]; p[j] = Math.max(sw.payer ? v : -v, 0) / pn[j]; }
      pay.push(p); reg.push({x1: Float64Array.from(S.x1), x2: Float64Array.from(S.x2), v: Float64Array.from(S.v)});
    });
    const n = pay.length, N = pay[0].length, valeur = Float64Array.from(pay[n - 1]), stoch = sv.xi !== 0;
    for (let j = n - 2; j >= 0; j--) {
      const itm = []; for (let i = 0; i < N; i++) if (pay[j][i] > 1e-14) itm.push(i);
      const p = stoch ? 14 : 10; if (itm.length < p) continue;
      const m = itm.length, r = reg[j];
      let m1 = 0, m2 = 0, m3 = 0; itm.forEach(i => { m1 += r.x1[i]; m2 += r.x2[i]; m3 += r.v[i]; }); m1 /= m; m2 /= m; m3 /= m;
      let s1 = 0, s2 = 0, s3 = 0; itm.forEach(i => { s1 += (r.x1[i] - m1) ** 2; s2 += (r.x2[i] - m2) ** 2; s3 += (r.v[i] - m3) ** 2; });
      s1 = Math.sqrt(s1 / m) + 1e-12; s2 = Math.sqrt(s2 / m) + 1e-12; s3 = Math.sqrt(s3 / m) + 1e-12;
      const A = new Float64Array(m * p), y = new Float64Array(m);
      itm.forEach((i, q) => {
        const u = r.x1[i] / s1, v = r.x2[i] / s2, w = (r.v[i] - 1) / s3, row = [1, u, v, u * u, u * v, v * v, u ** 3, u * u * v, u * v * v, v ** 3];
        if (stoch) row.push(w, w * u, w * v, w * w);
        for (let k = 0; k < p; k++) A[q * p + k] = row[k];
        y[q] = valeur[i];
      });
      const beta = H.moindresCarres(A, y, m, p);
      itm.forEach((i, q) => { let c = 0; for (let k = 0; k < p; k++) c += beta[k] * A[q * p + k]; if (pay[j][i] > c) valeur[i] = pay[j][i]; });
    }
    let s = 0; for (let i = 0; i < N; i++) s += valeur[i];
    return M.courbe.discount(tn) * s / N;
  }
  function prixBermudeenSV(M, sv, berm, prixExact, nPaths = 100000, graine = 5) {
    const a = lsmcSV(M, sv, berm, nPaths, graine), z = lsmcSV(M, ZERO, berm, nPaths, graine);
    return {prix: prixExact + (a - z), lsmcSV: a, lsmc0: z};
  }

  const API = {ZERO, grille, simuler, bondArr, europeens, volsSV, pointsSmile, erreursFlat, MoteurSmile, svDepuisX, xDepuisSV, ajusterSV, volsCaps, rebootstrap, lsmcSV, prixBermudeenSV, KAPPA, XI, RHOV_MAX};
  if (node) module.exports = API; else racine.HJM3 = API;
})(typeof self !== "undefined" ? self : this);
