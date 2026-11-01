/* hjm2.js — JavaScript port of the Python package hjm2f/ (two-factor Gaussian HJM), built on hjm.js (window.HJM).
   Runs in the browser (window.HJM2) and under Node (module.exports) for the tests.

   Model: forward-rate volatilities sigma_k(t) exp(-lambda_k (T - t)), k = 1, 2, correlation rho between the two
   Brownian motions. State variables X^k_t = int_0^t sigma_k(u) exp(lambda_k u) dW^k_u are driftless Gaussian and

       P(t,T) = P(0,T)/P(0,t) exp( -D_1(t,T) X^1_t - D_2(t,T) X^2_t + det(t,T) ).

   What the second factor removes: with one factor the swap value is monotone in the factor, which gives Jamshidian's
   decomposition into zero-coupon-bond options. With two factors the exercise boundary is a curve of the plane (x1, x2):
   European swaptions are priced by a 2-D integral, Bermudans by backward induction on a 2-D grid. */
(function (racine) {
  "use strict";
  const H = (typeof module !== "undefined" && module.exports) ? require("./hjm.js") : racine.HJM;

  // ------------------------------------------------------------------ Gauss-Legendre nodes (Newton)
  function gaussLegendre(n) {
    const x = new Float64Array(n), w = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), pp = 1;
      for (let it = 0; it < 100; it++) {
        let p1 = 1, p2 = 0;
        for (let j = 1; j <= n; j++) { const p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
        pp = n * (z * p1 - p2) / (z * z - 1);
        const z1 = z; z = z1 - p1 / pp;
        if (Math.abs(z - z1) < 1e-15) break;
      }
      x[i] = z; w[i] = 2 / ((1 - z * z) * pp * pp);
    }
    return {n, x, w};
  }
  const GL = gaussLegendre(24);
  const B = (l, t) => Math.abs(l) < 1e-12 ? t : -Math.expm1(-l * t) / l;

  // L with L L' = M for a symmetric positive semi-definite 2x2 matrix [[a, b], [b, c]] (eigen-decomposition)
  function racine2x2(a, b, c) {
    const th = 0.5 * Math.atan2(2 * b, a - c), cs = Math.cos(th), sn = Math.sin(th);
    const l1 = Math.max(a * cs * cs + 2 * b * cs * sn + c * sn * sn, 0), l2 = Math.max(a * sn * sn - 2 * b * cs * sn + c * cs * cs, 0);
    const s1 = Math.sqrt(l1), s2 = Math.sqrt(l2);
    return [[cs * s1, -sn * s2], [sn * s1, cs * s2]]; // columns: sqrt(l1) v1, sqrt(l2) v2
  }
  // uniform grid on [-demi, demi] with normalised normal weights (a fine uniform rule integrates a kinked payoff
  // better than Gauss-Hermite, whose nodes are far apart near the centre)
  function grilleNormale(n, demi = 6) {
    const z = new Float64Array(n), w = new Float64Array(n);
    let s = 0;
    for (let i = 0; i < n; i++) { z[i] = -demi + 2 * demi * i / (n - 1); w[i] = Math.exp(-0.5 * z[i] * z[i]); s += w[i]; }
    for (let i = 0; i < n; i++) w[i] /= s;
    return {z, w};
  }

  // ------------------------------------------------------------------ two-factor model (hjm_2f.py)
  function Modele2(courbe, lam, bucketEnds, sig1, sig2, rho) {
    if (bucketEnds.length !== sig1.length || bucketEnds.length !== sig2.length) throw new Error("one volatility per segment and per factor expected");
    for (let i = 1; i < bucketEnds.length; i++) if (!(bucketEnds[i] > bucketEnds[i - 1])) throw new Error("segment ends must be strictly increasing");
    if (!(rho >= -1 && rho <= 1)) throw new Error("rho must lie in [-1, 1]");
    const rk = (k, l) => k === l ? 1 : rho, cache = new Map();
    const n = bucketEnds.length;
    function segments(t) { // pieces of [0, t]; the last segment's volatility extends beyond its end
      const out = []; let left = 0;
      for (let i = 0; i < n; i++) {
        const end = bucketEnds[i], right = (i === n - 1 || end >= t) ? t : end;
        if (right > left) out.push([left, right, [sig1[i], sig2[i]]]);
        if (right >= t) break;
        left = right;
      }
      return out;
    }
    function integrer(t, fn) {
      let tot = 0;
      for (const [a, b, s] of segments(t)) {
        const demi = 0.5 * (b - a), mil = 0.5 * (a + b);
        for (let i = 0; i < GL.n; i++) tot += demi * GL.w[i] * fn(mil + demi * GL.x[i], s);
      }
      return tot;
    }
    function AG(t) {
      let r = cache.get("AG" + t);
      if (!r) {
        const A = [[0, 0], [0, 0]], G = [[0, 0], [0, 0]];
        for (let k = 0; k < 2; k++) for (let l = 0; l < 2; l++) {
          const lk = lam[k], ll = lam[l], c = rk(k, l);
          A[k][l] = c * integrer(t, (u, s) => s[k] * s[l] * Math.exp((lk + ll) * u));
          G[k][l] = c * integrer(t, (u, s) => s[k] * s[l] * Math.exp((lk + ll) * u) * B(ll, u));
        }
        r = {A, G}; cache.set("AG" + t, r);
      }
      return r;
    }
    const cov = t => AG(t).A;
    // P(t,T;x) = exp(c - d1 x1 - d2 x2)
    function coeffs(t, T) {
      const key = "C" + t + "|" + T; let r = cache.get(key);
      if (!r) {
        const {A, G} = AG(t), bt = [B(lam[0], t), B(lam[1], t)], bT = [B(lam[0], T), B(lam[1], T)], d = [bT[0] - bt[0], bT[1] - bt[1]];
        let det = 0;
        for (let k = 0; k < 2; k++) for (let l = 0; l < 2; l++) det += (bT[k] * bT[l] - bt[k] * bt[l]) * A[k][l] - d[k] * G[k][l] - d[l] * G[l][k];
        det *= -0.5;
        r = {c: Math.log(courbe.discount(T) / courbe.discount(t)) + det, d1: d[0], d2: d[1]};
        cache.set(key, r);
      }
      return r;
    }
    const zc = (t, T, x1, x2) => { if (T <= t) return 1; const q = coeffs(t, T); return Math.exp(q.c - q.d1 * x1 - q.d2 * x2); };
    // r_t = f(0,t) + sum_k b_k X^k + (convexity)  (Gaussian: used for the short-rate fan)
    function tauxCourtFixe(t) {
      const {A, G} = AG(t), b = [Math.exp(-lam[0] * t), Math.exp(-lam[1] * t)], bb = [B(lam[0], t), B(lam[1], t)];
      let conv = 0;
      for (let k = 0; k < 2; k++) for (let l = 0; l < 2; l++) conv += (b[k] * bb[l] + bb[k] * b[l]) * A[k][l] - b[k] * G[k][l] - b[l] * G[l][k];
      return courbe.instFwd(t) + 0.5 * conv;
    }
    // E[X_t] under the T-forward measure (deterministic drift)
    function moyenneForward(t, T) {
      const key = "M" + t + "|" + T; let r = cache.get(key);
      if (!r) {
        r = [0, 0];
        for (let k = 0; k < 2; k++) for (let l = 0; l < 2; l++) {
          const lk = lam[k], ll = lam[l];
          r[k] -= rk(k, l) * integrer(t, (u, s) => s[k] * s[l] * Math.exp(lk * u) * B(ll, T - u));
        }
        cache.set(key, r);
      }
      return r;
    }
    // value at expiry of the underlying swap, on arrays of states
    function valeursSwap(sw, X1, X2) {
      const m = X1.length, out = new Float64Array(m);
      if (!sw.rest.length) return out;
      const pn = coeffs(sw.expiry, sw.terminal), qi = sw.rest.map(T => coeffs(sw.expiry, T));
      for (let j = 0; j < m; j++) {
        const x1 = X1[j], x2 = X2[j];
        let ann = 0;
        for (let i = 0; i < qi.length; i++) ann += sw.acc[i] * Math.exp(qi[i].c - qi[i].d1 * x1 - qi[i].d2 * x2);
        const p = 1 - Math.exp(pn.c - pn.d1 * x1 - pn.d2 * x2) - sw.strike * ann;
        out[j] = sw.payer ? p : -p;
      }
      return out;
    }
    // European swaption: P(0,T) E^T[(swap value)^+] by 2-D quadrature of the Gaussian X_T (this replaces Jamshidian)
    function prixEuropeen(sw, nGrille = 241) {
      if (!sw.rest.length) return 0;
      const T = sw.expiry, m = moyenneForward(T, T), C = cov(T), L = racine2x2(C[0][0], C[0][1], C[1][1]);
      const g = grilleNormale(nGrille), X1 = new Float64Array(nGrille * nGrille), X2 = new Float64Array(nGrille * nGrille), W = new Float64Array(nGrille * nGrille);
      for (let i = 0; i < nGrille; i++) for (let j = 0; j < nGrille; j++) {
        const q = i * nGrille + j, z1 = g.z[i], z2 = g.z[j];
        X1[q] = m[0] + L[0][0] * z1 + L[0][1] * z2; X2[q] = m[1] + L[1][0] * z1 + L[1][1] * z2; W[q] = g.w[i] * g.w[j];
      }
      const v = valeursSwap(sw, X1, X2);
      let s = 0;
      for (let q = 0; q < v.length; q++) if (v[q] > 0) s += W[q] * v[q];
      return courbe.discount(T) * s;
    }
    // E^T[P(t,b)/P(t,T)] - P(0,b)/P(0,T): zero if the bond formula and the drift are consistent (no arbitrage)
    function defautMartingale(t, b, T) {
      const qb = coeffs(t, b), qT = coeffs(t, T), dd = [qb.d1 - qT.d1, qb.d2 - qT.d2], m = moyenneForward(t, T), C = cov(t);
      const expo = (qb.c - qT.c) - (dd[0] * m[0] + dd[1] * m[1]) + 0.5 * (dd[0] * dd[0] * C[0][0] + 2 * dd[0] * dd[1] * C[0][1] + dd[1] * dd[1] * C[1][1]);
      return Math.exp(expo) - courbe.discount(b) / courbe.discount(T);
    }
    // volatility of the forward rate f(t, T) at time t: two factor loadings
    const sigmaA = t => {
      if (t <= 0) return [sig1[0], sig2[0]];
      let lo = 0, hi = n; while (lo < hi) { const mm = (lo + hi) >> 1; if (bucketEnds[mm] < t) lo = mm + 1; else hi = mm; }
      const i = Math.min(lo, n - 1); return [sig1[i], sig2[i]];
    };
    function chargements(t, T) { const s = sigmaA(t); return [s[0] * Math.exp(-lam[0] * (T - t)), s[1] * Math.exp(-lam[1] * (T - t))]; }
    function correlationTaux(t, T1, T2) { // correlation of the instantaneous forward rates f(t,T1), f(t,T2)
      const a = chargements(t, T1), b = chargements(t, T2), c = (u, v) => u[0] * v[0] + rho * (u[0] * v[1] + u[1] * v[0]) + u[1] * v[1];
      const den = Math.sqrt(c(a, a) * c(b, b));
      return den > 0 ? c(a, b) / den : 1;
    }
    return {courbe, lam, rho, bucketEnds, sig1, sig2, cov, coeffs, zc, tauxCourtFixe, moyenneForward, valeursSwap, prixEuropeen, defautMartingale,
      sigmaA, chargements, correlationTaux};
  }

  // ------------------------------------------------------------------ swaption helpers
  function swaptionATM(courbe, expiry, tenor) { // annual-pay payer swaption struck at the forward swap rate
    const n = Math.round(tenor), pay = [], acc = [];
    for (let i = 1; i <= n; i++) { pay.push(expiry + i); acc.push(1); }
    const k = H.tauxSwapForward(courbe, H.Swaption(expiry, pay[n - 1], 0, pay, acc, true));
    return H.Swaption(expiry, pay[n - 1], k, pay, acc, true);
  }

  // ------------------------------------------------------------------ co-terminal bootstrap (calibration_2f.py)
  // sigma_2(t) = ratio * sigma_1(t); the volatility level is fitted exactly on the co-terminal European swaptions.
  function bootstrap2(courbe, lam, rho, ratio, bucketEnds, instruments, lower = 1e-6, upper = 0.2, nGrille = 241) {
    const sig = bucketEnds.map(() => 0.01), modeles = [], marche = [];
    const M = s => Modele2(courbe, lam, bucketEnds, s, s.map(v => ratio * v), rho);
    instruments.forEach((inst, i) => {
      const pm = H.prixBachelier(courbe, inst.swaption, inst.volN);
      const obj = s => { const v = sig.slice(); v[i] = s; return M(v).prixEuropeen(inst.swaption, nGrille) - pm; };
      try { sig[i] = H.brent(obj, lower, upper, 1e-10, 200); }
      catch (e) { throw new Error(`instrument ${i + 1} (expiry ${inst.swaption.expiry}): no volatility in [${lower}, ${upper}] reproduces the market price (${e.message})`); }
      modeles.push(M(sig).prixEuropeen(inst.swaption, nGrille)); marche.push(pm);
    });
    return {sigmas: sig, sigmas2: sig.map(v => ratio * v), modeles, marche, modele: M(sig)};
  }

  // ------------------------------------------------------------------ Bermudan on a 2-D grid (quad_2f.py)
  // Backward induction at the exercise dates under the T_N-forward measure, on values deflated by P(t, T_N): the
  // transition between two dates is exactly Gaussian, so only the grid and the quadrature are approximations.
  function prixQuad(M, berm, opts = {}) {
    const ng = opts.nGrille || 61, nStd = opts.nStd || 4.5, nq = opts.nQuad || 81, profil = !!opts.profil;
    const ex = berm.exerciseTimes, N = ex.length, tn = berm.terminal;
    const moy = ex.map(t => M.moyenneForward(t, tn)), cov = ex.map(t => M.cov(t));
    const axes = moy.map((m, j) => {
      const s1 = Math.sqrt(Math.max(cov[j][0][0], 1e-18)), s2 = Math.sqrt(Math.max(cov[j][1][1], 1e-18));
      return {a1: H.linspace(m[0] - nStd * s1, m[0] + nStd * s1, ng), a2: H.linspace(m[1] - nStd * s2, m[1] + nStd * s2, ng)};
    });
    const g = grilleNormale(nq), nq2 = nq * nq, wq = new Float64Array(nq2), z1 = new Float64Array(nq2), z2 = new Float64Array(nq2);
    for (let i = 0; i < nq; i++) for (let j = 0; j < nq; j++) { const q = i * nq + j; wq[q] = g.w[i] * g.w[j]; z1[q] = g.z[i]; z2[q] = g.z[j]; }
    const payoff = j => {
      const X1 = new Float64Array(ng * ng), X2 = new Float64Array(ng * ng), {a1, a2} = axes[j];
      for (let i = 0; i < ng; i++) for (let k = 0; k < ng; k++) { X1[i * ng + k] = a1[i]; X2[i * ng + k] = a2[k]; }
      const v = M.valeursSwap(berm.sousJacente(ex[j]), X1, X2), q = M.coeffs(ex[j], tn);
      for (let r = 0; r < v.length; r++) v[r] = Math.max(v[r], 0) / Math.exp(q.c - q.d1 * X1[r] - q.d2 * X2[r]);
      return v;
    };
    let U = payoff(N - 1);
    const tranches = [];
    if (profil) tranches[N - 1] = {t: ex[N - 1], a1: axes[N - 1].a1, a2: axes[N - 1].a2, exo: Uint8Array.from(U, v => v > 0 ? 1 : 0)};
    // interpolation of U (on the axes of date j+1) at the points a + increment
    const interp = (U1, ax, p1, p2) => { // p1, p2: arrays of coordinates -> index/fraction pairs
      const h1 = ax.a1[1] - ax.a1[0], h2 = ax.a2[1] - ax.a2[0], n = p1.length, i1 = new Int32Array(n), f1 = new Float64Array(n), i2 = new Int32Array(n), f2 = new Float64Array(n);
      for (let r = 0; r < n; r++) {
        let u = (p1[r] - ax.a1[0]) / h1, v = (p2[r] - ax.a2[0]) / h2;
        u = Math.min(Math.max(u, 0), ng - 1 - 1e-9); v = Math.min(Math.max(v, 0), ng - 1 - 1e-9);
        const iu = Math.floor(u), iv = Math.floor(v); i1[r] = iu; f1[r] = u - iu; i2[r] = iv; f2[r] = v - iv;
      }
      return {i1, f1, i2, f2};
    };
    for (let j = N - 2; j >= 0; j--) {
      const C = [[cov[j + 1][0][0] - cov[j][0][0], cov[j + 1][0][1] - cov[j][0][1]], [0, cov[j + 1][1][1] - cov[j][1][1]]];
      const L = racine2x2(C[0][0], C[0][1], C[1][1]), dm = [moy[j + 1][0] - moy[j][0], moy[j + 1][1] - moy[j][1]];
      const inc1 = new Float64Array(nq2), inc2 = new Float64Array(nq2);
      for (let q = 0; q < nq2; q++) { inc1[q] = dm[0] + L[0][0] * z1[q] + L[0][1] * z2[q]; inc2[q] = dm[1] + L[1][0] * z1[q] + L[1][1] * z2[q]; }
      const {a1, a2} = axes[j], ax = axes[j + 1];
      // the x1 part (row i) and the x2 part (column k) of the interpolation are independent: tabulate them once
      const P1 = new Float64Array(ng * nq2), P2 = new Float64Array(ng * nq2);
      for (let i = 0; i < ng; i++) for (let q = 0; q < nq2; q++) { P1[i * nq2 + q] = a1[i] + inc1[q]; P2[i * nq2 + q] = a2[i] + inc2[q]; }
      const T1 = interp(U, ax, P1, new Float64Array(P1.length)), T2 = interp(U, ax, new Float64Array(P2.length), P2);
      const cont = new Float64Array(ng * ng), pay = payoff(j), Unew = new Float64Array(ng * ng);
      for (let i = 0; i < ng; i++) for (let k = 0; k < ng; k++) {
        let s = 0; const o1 = i * nq2, o2 = k * nq2;
        for (let q = 0; q < nq2; q++) {
          const iu = T1.i1[o1 + q], fu = T1.f1[o1 + q], iv = T2.i2[o2 + q], fv = T2.f2[o2 + q], b = iu * ng + iv;
          s += wq[q] * (U[b] * (1 - fu) * (1 - fv) + U[b + ng] * fu * (1 - fv) + U[b + 1] * (1 - fu) * fv + U[b + ng + 1] * fu * fv);
        }
        cont[i * ng + k] = s;
      }
      const exo = profil ? new Uint8Array(ng * ng) : null;
      for (let r = 0; r < Unew.length; r++) { Unew[r] = Math.max(pay[r], cont[r]); if (exo) exo[r] = pay[r] > cont[r] ? 1 : 0; }
      if (profil) tranches[j] = {t: ex[j], a1: axes[j].a1, a2: axes[j].a2, exo};
      U = Unew;
    }
    const C0 = moy.length ? M.cov(ex[0]) : null, L0 = racine2x2(C0[0][0], C0[0][1], C0[1][1]);
    const P1 = new Float64Array(nq2), P2 = new Float64Array(nq2);
    for (let q = 0; q < nq2; q++) { P1[q] = moy[0][0] + L0[0][0] * z1[q] + L0[0][1] * z2[q]; P2[q] = moy[0][1] + L0[1][0] * z1[q] + L0[1][1] * z2[q]; }
    const T = interp(U, axes[0], P1, P2);
    let v0 = 0;
    for (let q = 0; q < nq2; q++) {
      const iu = T.i1[q], fu = T.f1[q], iv = T.i2[q], fv = T.f2[q], b = iu * ng + iv;
      v0 += wq[q] * (U[b] * (1 - fu) * (1 - fv) + U[b + ng] * fu * (1 - fv) + U[b + 1] * (1 - fu) * fv + U[b + ng + 1] * fu * fv);
    }
    const prix = M.courbe.discount(tn) * v0;
    return profil ? {prix, tranches} : prix;
  }

  // ------------------------------------------------------------------ Longstaff-Schwartz in two factors (mc_2f.py)
  // Under the T_N-forward measure, values deflated by P(t, T_N) are martingales: no path-dependent discounting.
  function prixLSMC2(M, berm, nChemins = 20000, graine = 1234, opts = {}) {
    const ex = berm.exerciseTimes, n = ex.length, tn = berm.terminal, N = nChemins - (nChemins % 2), rng = H.generateur(graine);
    const degre = Math.max(1, Math.min(4, Math.round(opts.degre || 2)));
    const monomes = []; for (let a = 0; a <= degre; a++) for (let b = 0; a + b <= degre; b++) monomes.push([a, b]);
    const p = monomes.length;
    const moy = ex.map(t => M.moyenneForward(t, tn)), cov = ex.map(t => M.cov(t));
    const ecart = cov.map(c => [Math.sqrt(Math.max(c[0][0], 1e-300)), Math.sqrt(Math.max(c[1][1], 1e-300))]);
    function simuler() {
      const x1 = new Float64Array(N * n), x2 = new Float64Array(N * n);
      let pm = [0, 0], pc = [[0, 0], [0, 0]];
      const cur1 = new Float64Array(N), cur2 = new Float64Array(N);
      for (let j = 0; j < n; j++) {
        const dC = [cov[j][0][0] - pc[0][0], cov[j][0][1] - pc[0][1], cov[j][1][1] - pc[1][1]], L = racine2x2(dC[0], dC[1], dC[2]);
        for (let i = 0; i < N / 2; i++) {
          const a = rng.normale(), b = rng.normale();
          for (const [idx, s] of [[i, 1], [i + N / 2, -1]]) {
            cur1[idx] += (moy[j][0] - pm[0]) + s * (L[0][0] * a + L[0][1] * b);
            cur2[idx] += (moy[j][1] - pm[1]) + s * (L[1][0] * a + L[1][1] * b);
          }
        }
        for (let i = 0; i < N; i++) { x1[i * n + j] = cur1[i]; x2[i * n + j] = cur2[i]; }
        pm = moy[j]; pc = cov[j];
      }
      return {x1, x2};
    }
    const col = (X, j) => { const c = new Float64Array(N); for (let i = 0; i < N; i++) c[i] = X[i * n + j]; return c; };
    const payoffs = S => ex.map((t, j) => {
      const X1 = col(S.x1, j), X2 = col(S.x2, j), v = M.valeursSwap(berm.sousJacente(t), X1, X2), q = M.coeffs(t, tn);
      for (let i = 0; i < N; i++) v[i] = Math.max(v[i], 0) / Math.exp(q.c - q.d1 * X1[i] - q.d2 * X2[i]);
      return v;
    });
    const base = (u, v, k) => Math.pow(u, monomes[k][0]) * Math.pow(v, monomes[k][1]);
    const S = simuler(), pay = payoffs(S), valeur = Float64Array.from(pay[n - 1]), betas = new Array(n).fill(null);
    for (let j = n - 2; j >= 0; j--) {
      const x1 = col(S.x1, j), x2 = col(S.x2, j), itm = [];
      for (let i = 0; i < N; i++) if (pay[j][i] > 1e-14) itm.push(i);
      if (itm.length < p) continue;
      const m = itm.length, A = new Float64Array(m * p), y = new Float64Array(m);
      itm.forEach((i, q) => { for (let k = 0; k < p; k++) A[q * p + k] = base(x1[i] / ecart[j][0], x2[i] / ecart[j][1], k); y[q] = valeur[i]; });
      const beta = H.moindresCarres(A, y, m, p); betas[j] = beta;
      itm.forEach((i, q) => { let c = 0; for (let k = 0; k < p; k++) c += beta[k] * A[q * p + k]; if (pay[j][i] > c) valeur[i] = pay[j][i]; });
    }
    let SP = S, valP = valeur, payP = pay;
    const nbExerces = new Array(n).fill(0);
    if (opts.independants) {
      SP = simuler(); payP = payoffs(SP); valP = Float64Array.from(payP[n - 1]);
      const fait = new Uint8Array(N);
      for (let j = 0; j < n - 1; j++) {
        if (!betas[j]) continue;
        const x1 = col(SP.x1, j), x2 = col(SP.x2, j);
        for (let i = 0; i < N; i++) {
          if (fait[i] || !(payP[j][i] > 1e-14)) continue;
          let c = 0; for (let k = 0; k < p; k++) c += betas[j][k] * base(x1[i] / ecart[j][0], x2[i] / ecart[j][1], k);
          if (payP[j][i] > c) { valP[i] = payP[j][i]; fait[i] = 1; nbExerces[j]++; }
        }
      }
      for (let i = 0; i < N; i++) if (!fait[i] && payP[n - 1][i] > 0) nbExerces[n - 1]++;
    } else {
      // exercise date of each path = the last date where the strategy value equals the payoff
      for (let i = 0; i < N; i++) { for (let j = 0; j < n; j++) if (valP[i] > 0 && Math.abs(valP[i] - payP[j][i]) < 1e-15 * Math.max(1, valP[i])) { nbExerces[j]++; break; } }
    }
    let s = 0, s2 = 0;
    for (let i = 0; i < N; i++) { s += valP[i]; s2 += valP[i] * valP[i]; }
    const P0 = M.courbe.discount(tn), moyenne = s / N, v = (s2 - N * moyenne * moyenne) / (N - 1);
    // martingale check at the first exercise date: E^{T_N}[P(t0, Ti)/P(t0, T_N)] = P(0, Ti)/P(0, T_N)
    const t0 = ex[0], X1 = col(SP.x1, 0), X2 = col(SP.x2, 0), qn = M.coeffs(t0, tn), martingale = [];
    berm.payTimes.filter(Ti => Ti > t0 && Ti < tn).forEach(Ti => {
      const q = M.coeffs(t0, Ti); let a = 0, a2 = 0;
      for (let i = 0; i < N; i++) { const r = Math.exp((q.c - qn.c) - (q.d1 - qn.d1) * X1[i] - (q.d2 - qn.d2) * X2[i]); a += r; a2 += r * r; }
      const mm = a / N;
      martingale.push({t: Ti, mc: mm, se: Math.sqrt(Math.max(a2 / N - mm * mm, 0) / N), courbe: M.courbe.discount(Ti) / M.courbe.discount(tn)});
    });
    return {prix: P0 * moyenne, erreur: P0 * Math.sqrt(Math.max(v, 0) / N), nbExerces, martingale, t0};
  }

  // ------------------------------------------------------------------ Nelder-Mead and global fit (optimize.py, calibration_global.py)
  function nelderMead(f, x0, step = 0.5, maxFev = 300, xtol = 1e-4, ftol = 1e-4, progres) {
    const n = x0.length; let simplex = [x0.slice()];
    for (let i = 0; i < n; i++) { const p = x0.slice(); p[i] += step; simplex.push(p); }
    let values = simplex.map(f), nfev = simplex.length;
    while (nfev < maxFev) {
      const ordre = values.map((v, i) => i).sort((a, b) => values[a] - values[b]);
      simplex = ordre.map(i => simplex[i]); values = ordre.map(i => values[i]);
      let taille = 0; for (let i = 1; i <= n; i++) for (let k = 0; k < n; k++) taille = Math.max(taille, Math.abs(simplex[i][k] - simplex[0][k]));
      if (taille < xtol && Math.abs(values[n] - values[0]) < ftol) break;
      if (progres) progres(Math.min(nfev / maxFev, 0.999), values[0]);
      const centre = Array.from({length: n}, (_, k) => { let s = 0; for (let i = 0; i < n; i++) s += simplex[i][k]; return s / n; }), pire = simplex[n];
      const refl = centre.map((c, k) => c + (c - pire[k])), fr = f(refl); nfev++;
      if (fr < values[0]) {
        const exp = centre.map((c, k) => c + 2 * (c - pire[k])), fe = f(exp); nfev++;
        if (fe < fr) { simplex[n] = exp; values[n] = fe; } else { simplex[n] = refl; values[n] = fr; }
      } else if (fr < values[n - 1]) { simplex[n] = refl; values[n] = fr; }
      else {
        const dehors = fr < values[n];
        const pt = dehors ? centre.map((c, k) => c + 0.5 * (refl[k] - c)) : centre.map((c, k) => c + 0.5 * (pire[k] - c)), fp = f(pt); nfev++;
        if (fp < Math.min(fr, values[n])) { simplex[n] = pt; values[n] = fp; }
        else for (let i = 1; i <= n; i++) { simplex[i] = simplex[i].map((v, k) => simplex[0][k] + 0.5 * (v - simplex[0][k])); values[i] = f(simplex[i]); nfev++; }
      }
    }
    let b = 0; for (let i = 1; i <= n; i++) if (values[i] < values[b]) b = i;
    return {x: simplex[b], f: values[b], nfev};
  }
  const LAM1_MAX = 0.5, LAM2_MAX = 3;
  const versParams = x => ({lam: [LAM1_MAX / (1 + Math.exp(-x[0])), LAM2_MAX / (1 + Math.exp(-x[1]))], rho: Math.tanh(x[2]), ratio: Math.exp(x[3])});
  const depuisParams = (lam, rho, ratio) => {
    const logit = (p, top) => { const q = Math.min(Math.max(p, 1e-9), top - 1e-9) / top; return Math.log(q / (1 - q)); };
    return [logit(lam[0], LAM1_MAX), logit(lam[1], LAM2_MAX), Math.atanh(Math.min(Math.max(rho, -0.999999), 0.999999)), Math.log(ratio)];
  };
  // errors (model - market, bp of normal vol) on the target swaptions [{sw, volN}] for given structural parameters
  function erreurs2(courbe, lam, rho, ratio, ex, instruments, cibles, nGrille) {
    const cal = bootstrap2(courbe, lam, rho, ratio, ex, instruments, 1e-6, 0.2, nGrille);
    const err = cibles.map(c => { const p = cal.modele.prixEuropeen(c.sw, nGrille); const v = p > 0 ? H.volNImplicite(courbe, c.sw, p) : 0; return (v - c.volN) * 1e4; });
    return {err, cal};
  }
  const rmse = a => Math.sqrt(a.reduce((s, v) => s + v * v, 0) / a.length);
  function ajuster2F(courbe, ex, instruments, cibles, x0, libres, nGrille = 61, maxFev = 200, progres) {
    const idx = libres.map((f, i) => f ? i : -1).filter(i => i >= 0), historique = [];
    const compose = y => { const x = x0.slice(); idx.forEach((i, k) => { x[i] = y[k]; }); return x; };
    const obj = y => {
      const q = versParams(compose(y));
      try { const r = rmse(erreurs2(courbe, q.lam, q.rho, q.ratio, ex, instruments, cibles, nGrille).err); historique.push(r); return r; } catch (e) { return 1e3; }
    };
    const res = nelderMead(obj, idx.map(i => x0[i]), 0.5, maxFev, 1e-4, 1e-4, progres);
    const q = versParams(compose(res.x)), fin = erreurs2(courbe, q.lam, q.rho, q.ratio, ex, instruments, cibles, nGrille);
    return {lam: q.lam, rho: q.rho, ratio: q.ratio, rmse: rmse(fin.err), erreurs: fin.err, sigmas: fin.cal.sigmas, evaluations: res.nfev, historique};
  }
  function ajuster1F(courbe, ex, instruments, cibles, lamMax = 0.5) {
    const erreurs1 = lam => {
      const cal = H.bootstrap(courbe, lam, ex, instruments, 1e-6, 0.2, "corrigee");
      return cibles.map(c => { const p = cal.modele.jamshidian(c.sw); const v = p > 0 ? H.volNImplicite(courbe, c.sw, p) : 0; return (v - c.volN) * 1e4; });
    };
    const f = lam => { try { return rmse(erreurs1(lam)); } catch (e) { return 1e3; } };
    let a = 0, b = lamMax; const g = (Math.sqrt(5) - 1) / 2;
    let c = b - g * (b - a), d = a + g * (b - a), fc = f(c), fd = f(d);
    for (let it = 0; it < 30; it++) {
      if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = f(c); } else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = f(d); }
    }
    const cand = [[f(a), a], [fc, c], [fd, d], [f(b), b]].sort((u, v) => u[0] - v[0])[0];
    return {lam: cand[1], rmse: cand[0], erreurs: erreurs1(cand[1])};
  }

  const API = {gaussLegendre, racine2x2, grilleNormale, Modele2, swaptionATM, bootstrap2, prixQuad, prixLSMC2, nelderMead, versParams, depuisParams,
    erreurs2, ajuster2F, ajuster1F, rmse};
  if (typeof module !== "undefined" && module.exports) module.exports = API; else racine.HJM2 = API;
})(typeof self !== "undefined" ? self : this);
