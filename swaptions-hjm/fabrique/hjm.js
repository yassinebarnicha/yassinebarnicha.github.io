/* hjm.js — faithful JavaScript port of the Python package hjm1f/ (same model, same formulas, same schemes).
   Runs in the browser (window.HJM) and under Node (module.exports) for the tests.

   The only deliberate differences with the Python:
   - the normal law is computed with erf/erfc by series and continued fraction (accuracy ~1e-15);
   - the Monte Carlo has its own generator (xoshiro128** + Box-Muller), so not the same draws as NumPy;
   - the European receiver swaption is priced with calls on zero-coupon bonds (the Python returns the payer formula
     in that case: the code itself says the receiver is not specialised there).
   Options added after the model review: corrected convexity sign, exact Monte Carlo discounting, Crank–Nicolson with
   Rannacher steps, regression degree and independent regression paths (see the Theory tab and the report). */
(function (racine) {
  "use strict";

  // ------------------------------------------------------------------ normal law
  const SQRTPI = Math.sqrt(Math.PI);
  function erfc(x) {
    if (x < 0) return 2 - erfc(-x);
    if (x < 3) return 1 - erfPos(x);
    // continued fraction (Lentz): erfc(x) = exp(-x²)/√π · 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + …))))
    const tiny = 1e-300;
    let f = x, C = x, D = 0;
    for (let n = 1; n < 300; n++) {
      const a = n / 2;
      D = x + a * D; if (Math.abs(D) < tiny) D = tiny; D = 1 / D;
      C = x + a / C; if (Math.abs(C) < tiny) C = tiny;
      const delta = C * D; f *= delta;
      if (Math.abs(delta - 1) < 1e-16) break;
    }
    return Math.exp(-x * x) / SQRTPI / f;
  }
  // erf(x) = 2/√π e^{-x²} Σ 2^n x^{2n+1} / (1·3·…·(2n+1)): all terms positive, no cancellation
  function erfPos(x) {
    let terme = x, somme = x;
    for (let n = 1; n < 200; n++) {
      terme *= 2 * x * x / (2 * n + 1);
      somme += terme;
      if (terme < somme * 1e-17) break;
    }
    return 2 / SQRTPI * Math.exp(-x * x) * somme;
  }
  const erf = x => x < 0 ? -erf(-x) : (x < 3 ? erfPos(x) : 1 - erfc(x));
  const phi = x => 0.5 * erfc(-x / Math.SQRT2);
  const pdf = x => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);

  // ------------------------------------------------------------------ Brent solver (solvers.py)
  const ULP1 = Math.pow(2, -52);
  function ulp(x) { // math.ulp(1.0) in the Python: a constant
    return ULP1 * (x === x ? 1 : 1);
  }
  function brent(f, a, b, tol = 1e-10, maxIter = 100) {
    let fa = f(a), fb = f(b);
    if (fa === 0) return a;
    if (fb === 0) return b;
    if (fa * fb > 0) throw new Error("root not bracketed: f(" + a + ") and f(" + b + ") have the same sign");
    let c = a, fc = fa, d = b - a, e = d;
    for (let it = 0; it < maxIter; it++) {
      if (fb * fc > 0) { c = a; fc = fa; d = e = b - a; }
      if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
      const tol1 = 2 * ulp(1) * Math.abs(b) + 0.5 * tol;
      const xm = 0.5 * (c - b);
      if (Math.abs(xm) <= tol1 || fb === 0) return b;
      if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
        const s = fb / fa;
        let p, q;
        if (a === c) { p = 2 * xm * s; q = 1 - s; }
        else { const qq = fa / fc, r = fb / fc; p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1)); q = (qq - 1) * (r - 1) * (s - 1); }
        if (p > 0) q = -q;
        p = Math.abs(p);
        const min1 = 3 * xm * q - Math.abs(tol1 * q), min2 = Math.abs(e * q);
        if (2 * p < Math.min(min1, min2)) { e = d; d = p / q; } else { d = xm; e = d; }
      } else { d = xm; e = d; }
      a = b; fa = fb;
      if (Math.abs(d) > tol1) b += d; else b += xm >= 0 ? tol1 : -tol1;
      fb = f(b);
    }
    throw new Error("Brent did not converge");
  }

  // ------------------------------------------------------------------ zero-coupon curve (curve.py)
  function interpLin(x, xs, ys) {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
    let lo = 0, hi = xs.length - 1; // bisect_right - 1
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (xs[m] <= x) lo = m; else hi = m; }
    const w = (x - xs[lo]) / (xs[lo + 1] - xs[lo]);
    return ys[lo] * (1 - w) + ys[lo + 1] * w;
  }
  function Courbe(times, discounts) {
    times = times.slice(); discounts = discounts.slice();
    if (times.length !== discounts.length) throw new Error("maturities and discount factors have different lengths");
    if (times.length < 1) throw new Error("empty curve");
    for (let i = 1; i < times.length; i++) if (!(times[i] > times[i - 1])) throw new Error("the curve maturities must be strictly increasing");
    if (times[0] > 0) { times.unshift(0); discounts.unshift(1); }
    if (discounts[0] !== 1) throw new Error("the discount factor at 0 must be 1");
    for (const d of discounts) if (!(d > 0)) throw new Error("discount factors must be strictly positive");
    const logs = discounts.map(d => Math.log(Math.max(d, 1e-16)));
    const discount = t => t <= 0 ? 1 : Math.exp(interpLin(t, times, logs));
    const instFwd = (t, h = 1e-5) => { const t0 = Math.max(0, t - h), t1 = t + h; return -(Math.log(discount(t1)) - Math.log(discount(t0))) / (t1 - t0); };
    return {times, discounts, discount, instFwd};
  }

  // ------------------------------------------------------------------ instruments (instruments.py)
  function Swaption(expiry, terminal, strike, payTimes, accruals, payer = true) {
    const rest = [], acc = [];
    payTimes.forEach((t, i) => { if (t > expiry) { rest.push(t); acc.push(accruals[i]); } });
    return {expiry, terminal, strike, payTimes, accruals, payer, rest, acc};
  }
  function Bermudeenne(exerciseTimes, terminal, strike, payTimes, accruals, payer = true) {
    return {exerciseTimes, terminal, strike, payTimes, accruals, payer,
      sousJacente: t => Swaption(t, terminal, strike, payTimes, accruals, payer)};
  }

  // ------------------------------------------------------------------ market quotes (market_quotes.py)
  const annuite = (courbe, times, acc) => times.reduce((s, t, i) => s + acc[i] * courbe.discount(t), 0);
  function tauxSwapForward(courbe, sw) {
    const ann = annuite(courbe, sw.rest, sw.acc);
    if (!(ann > 0)) throw new Error("zero annuity: no coupon after the expiry " + sw.expiry);
    return (courbe.discount(sw.expiry) - courbe.discount(sw.terminal)) / ann;
  }
  function prixBachelier(courbe, sw, volN) {
    const ann = annuite(courbe, sw.rest, sw.acc), fwd = tauxSwapForward(courbe, sw), k = sw.strike;
    const tau = Math.max(sw.expiry, 1e-12), s = Math.max(volN, 0) * Math.sqrt(tau);
    if (s < 1e-16) return ann * (sw.payer ? Math.max(fwd - k, 0) : Math.max(k - fwd, 0));
    const d = (fwd - k) / s;
    return sw.payer ? ann * ((fwd - k) * phi(d) + s * pdf(d)) : ann * ((k - fwd) * phi(-d) + s * pdf(d));
  }

  // ------------------------------------------------------------------ one-factor HJM model (hjm_1f.py)
  // convexite: "code" reproduces hjm_1f.py; "corrigee" (corrected) flips the sign of the convexity term c(t,T) in
  // P(t,T) = P(0,T)/P(0,t) · exp(−D(t,T) (c(t,T) + X_t)) and therefore in the short rate. Only the corrected version makes
  // discounted zero-coupon bonds martingales (no arbitrage): see the Theory tab.
  function Modele(courbe, lam, bucketEnds, sigmas, convexite = "code") {
    const sg = convexite === "corrigee" ? -1 : 1;
    if (bucketEnds.length !== sigmas.length) throw new Error("one volatility per segment expected");
    for (let i = 1; i < bucketEnds.length; i++) if (!(bucketEnds[i] > bucketEnds[i - 1])) throw new Error("segment ends must be strictly increasing");
    const petit = Math.abs(lam) < 1e-14;
    const sigma = t => {
      if (t <= 0) return sigmas[0];
      let lo = 0, hi = bucketEnds.length; // bisect_left
      while (lo < hi) { const m = (lo + hi) >> 1; if (bucketEnds[m] < t) lo = m + 1; else hi = m; }
      return sigmas[Math.min(lo, sigmas.length - 1)];
    };
    const Bbar = t => petit ? t : (1 - Math.exp(-lam * t)) / lam;
    const b = t => Math.exp(-lam * t);
    const D = (t, T) => Bbar(T) - Bbar(t);
    function alpha(t) {
      if (t <= 0) return 0;
      let tot = 0, left = 0;
      for (let i = 0; i < bucketEnds.length; i++) {
        const end = bucketEnds[i], s = sigmas[i], right = Math.min(end, t);
        if (right > left) tot += petit ? s * s * (right - left) : s * s * (Math.exp(2 * lam * right) - Math.exp(2 * lam * left)) / (2 * lam);
        if (end >= t) break;
        left = end;
      }
      return tot;
    }
    function gamma(t) {
      if (t <= 0) return 0;
      let tot = 0, left = 0;
      for (let i = 0; i < bucketEnds.length; i++) {
        const end = bucketEnds[i], s = sigmas[i], right = Math.min(end, t);
        if (right > left) tot += petit ? 0.5 * s * s * (right * right - left * left)
          : (s * s / lam) * ((Math.exp(2 * lam * right) - Math.exp(2 * lam * left)) / (2 * lam) - (Math.exp(lam * right) - Math.exp(lam * left)) / lam);
        if (end >= t) break;
        left = end;
      }
      return tot;
    }
    const eta = t => sigma(t) * Math.exp(lam * t);
    function zc(t, T, x) {
      if (T <= t) return 1;
      const d = D(t, T), c = 0.5 * (Bbar(T) + Bbar(t)) * alpha(t) - gamma(t);
      return courbe.discount(T) / courbe.discount(t) * Math.exp(d * (sg * c - x));
    }
    // r_t = f(0,t) - b(t)(B(t) α(t) - γ(t) - X_t)  (original sign); "partieFixe" gathers what does not depend on X
    const partieFixe = t => courbe.instFwd(t) - sg * b(t) * (Bbar(t) * alpha(t) - gamma(t));
    const tauxCourt = (t, x) => sg === 1 ? courbe.instFwd(t) - b(t) * (Bbar(t) * alpha(t) - gamma(t) - x) : partieFixe(t) + b(t) * x;
    function valeurSwap(sw, x) {
      if (!sw.rest.length) return 0;
      const pn = zc(sw.expiry, sw.terminal, x);
      let ann = 0;
      for (let i = 0; i < sw.rest.length; i++) ann += sw.acc[i] * zc(sw.expiry, sw.rest[i], x);
      const payeur = 1 - pn - sw.strike * ann;
      return sw.payer ? payeur : -payeur;
    }
    // fast version on a whole grid of x (same computation as valeurSwap, deterministic terms factored out)
    function valeursSwapGrille(sw, xs) {
      const out = new Float64Array(xs.length);
      if (!sw.rest.length) return out;
      const t = sw.expiry, at = alpha(t), gt = gamma(t), Bt = Bbar(t), Pt = courbe.discount(t);
      const facteur = T => { if (T <= t) return null; const d = D(t, T); return {d, k: courbe.discount(T) / Pt * Math.exp(d * sg * (0.5 * (Bbar(T) + Bt) * at - gt))}; };
      const fn = facteur(sw.terminal), fi = sw.rest.map(facteur);
      for (let j = 0; j < xs.length; j++) {
        const x = xs[j], pn = fn ? fn.k * Math.exp(-fn.d * x) : 1;
        let ann = 0;
        for (let i = 0; i < fi.length; i++) ann += sw.acc[i] * (fi[i] ? fi[i].k * Math.exp(-fi[i].d * x) : 1);
        const p = 1 - pn - sw.strike * ann;
        out[j] = sw.payer ? p : -p;
      }
      return out;
    }
    const volOptionZC = (te, tm) => Math.abs(D(te, tm)) * Math.sqrt(Math.max(alpha(te), 0));
    function putZC(te, tm, K) {
      if (te <= 0) return Math.max(K * courbe.discount(te) - courbe.discount(tm), 0);
      const p0t = courbe.discount(te), p0T = courbe.discount(tm), v = volOptionZC(te, tm);
      if (v < 1e-14) return Math.max(K * p0t - p0T, 0);
      const h = Math.log(p0T / (K * p0t)) / v + 0.5 * v;
      return K * p0t * phi(v - h) - p0T * phi(-h);
    }
    function callZC(te, tm, K) {
      if (te <= 0) return Math.max(courbe.discount(tm) - K * courbe.discount(te), 0);
      const p0t = courbe.discount(te), p0T = courbe.discount(tm), v = volOptionZC(te, tm);
      if (v < 1e-14) return Math.max(p0T - K * p0t, 0);
      const h = Math.log(p0T / (K * p0t)) / v + 0.5 * v;
      return p0T * phi(h) - K * p0t * phi(h - v);
    }
    function jamshidian(sw) {
      if (!sw.rest.length) return 0;
      const xStar = brent(x => valeurSwap(sw, x), -10, 10, 1e-12, 200);
      const opt = sw.payer ? putZC : callZC;
      let prix = 0;
      const n = sw.rest.length;
      for (let i = 0; i < n - 1; i++) prix += sw.strike * sw.acc[i] * opt(sw.expiry, sw.rest[i], zc(sw.expiry, sw.rest[i], xStar));
      prix += (1 + sw.strike * sw.acc[n - 1]) * opt(sw.expiry, sw.rest[n - 1], zc(sw.expiry, sw.rest[n - 1], xStar));
      return prix;
    }
    return {courbe, lam, bucketEnds, sigmas, convexite, sigma, Bbar, b, D, alpha, gamma, eta, zc, tauxCourt, partieFixe,
      valeurSwap, valeursSwapGrille, putZC, callZC, jamshidian};
  }

  // ------------------------------------------------------------------ co-terminal bootstrap (calibration_1f.py)
  function bootstrap(courbe, lam, bucketEnds, instruments, lower = 1e-6, upper = 0.20, convexite = "code") {
    const sig = bucketEnds.map(() => 0.01);
    const modeles = [], marche = [], parBrent = [];
    instruments.forEach((inst, i) => {
      const prixMarche = prixBachelier(courbe, inst.swaption, inst.volN);
      const obj = s => { const v = sig.slice(); v[i] = s; return Modele(courbe, lam, bucketEnds, v, convexite).jamshidian(inst.swaption) - prixMarche; };
      try { sig[i] = brent(obj, lower, upper, 1e-10, 200); }
      catch (e) { throw new Error(`instrument ${i + 1} (expiry ${inst.swaption.expiry}): no volatility in [${lower}, ${upper}] reproduces the market price (${e.message})`); }
      modeles.push(Modele(courbe, lam, bucketEnds, sig, convexite).jamshidian(inst.swaption));
      marche.push(prixMarche);
    });
    return {sigmas: sig, modeles, marche, modele: Modele(courbe, lam, bucketEnds, sig, convexite)};
  }

  // ------------------------------------------------------------------ PDE (pde_1f.py)
  function arrondiBanquier(x) { // Python 3 round()
    const f = Math.floor(x), r = x - f;
    if (r > 0.5) return f + 1;
    if (r < 0.5) return f;
    return f % 2 === 0 ? f : f + 1;
  }
  function linspace(a, b, n) {
    const out = new Float64Array(n);
    if (n === 1) { out[0] = a; return out; }
    const pas = (b - a) / (n - 1);
    for (let k = 0; k < n; k++) out[k] = a + k * pas;
    out[n - 1] = b;
    return out;
  }
  function grilleEcartsTypes(modele, tMax, nStd = 5, nx = 301, nt = 500) {
    const std = Math.sqrt(Math.max(modele.alpha(tMax), 1e-16));
    return {xMin: -nStd * std, xMax: nStd * std, nx, nt, tMax};
  }
  function grilleTemps(berm, grille) {
    const cles = [...new Set([0, grille.tMax, ...berm.exerciseTimes])].sort((a, b) => a - b);
    if (cles.length < 2) return [0, grille.tMax];
    const base = (cles[cles.length - 1] - cles[0]) / Math.max(grille.nt, 1);
    const ts = [cles[0]];
    for (let k = 0; k < cles.length - 1; k++) {
      const t0 = cles[k], t1 = cles[k + 1], n = Math.max(1, arrondiBanquier((t1 - t0) / base));
      const loc = linspace(t0, t1, n + 1);
      for (let j = 1; j < loc.length; j++) ts.push(loc[j]);
    }
    return ts;
  }
  const r12 = t => Math.round(t * 1e12) / 1e12;
  function thomas(a, b, c, d) {
    const n = d.length, cp = new Float64Array(n - 1), dp = new Float64Array(n);
    cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
    for (let i = 1; i < n - 1; i++) { const den = b[i] - a[i - 1] * cp[i - 1]; cp[i] = c[i] / den; dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / den; }
    dp[n - 1] = (d[n - 1] - a[n - 2] * dp[n - 2]) / (b[n - 1] - a[n - 2] * cp[n - 2]);
    const x = new Float64Array(n);
    x[n - 1] = dp[n - 1];
    for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
    return x;
  }
  // θ-scheme: θ = 1 implicit (the Python code), θ = ½ Crank–Nicolson. Rannacher: after the final date and after each
  // exercise date (non-smooth payoff), a few implicit steps damp the Crank–Nicolson oscillations.
  function prixEDP(modele, berm, grille, avecProfil = false, opts = {}) {
    if (grille.nx < 4) throw new Error("at least 4 space points");
    const theta0 = opts.theta == null ? 1 : opts.theta, rannacher = Math.max(0, Math.round(opts.rannacher || 0));
    const xs = linspace(grille.xMin, grille.xMax, grille.nx), dx = xs[1] - xs[0];
    const ts = grilleTemps(berm, grille);
    const exSet = new Set(berm.exerciseTimes.map(r12));
    const sousJ = new Map(berm.exerciseTimes.map(t => [t, berm.sousJacente(t)]));
    const exercice = t => { const v = modele.valeursSwapGrille(sousJ.get(t) || berm.sousJacente(t), xs); for (let j = 0; j < v.length; j++) v[j] = Math.max(v[j], 0); return v; };
    const dernier = berm.exerciseTimes[berm.exerciseTimes.length - 1];
    let v = exercice(dernier);
    const tranches = [];
    const m = grille.nx - 2, lower = new Float64Array(m - 1), upper = new Float64Array(m - 1), diag = new Float64Array(m), rhs = new Float64Array(m);
    // intrinsic value at both edges for each exercise date (does not depend on the time step)
    const bordsEx = berm.exerciseTimes.map(ex => {
      const sw = sousJ.get(ex);
      return {ex, g: Math.max(modele.valeurSwap(sw, xs[0]), 0), d: Math.max(modele.valeurSwap(sw, xs[grille.nx - 1]), 0)};
    });
    let lisser = rannacher;
    for (let n = ts.length - 2; n >= 0; n--) {
      const th = lisser > 0 ? 1 : theta0;
      if (lisser > 0) lisser--;
      const t = ts[n], dt = ts[n + 1] - ts[n], A = modele.partieFixe(t), bt = modele.b(t);
      let diff;
      if (theta0 === 1) { // the Python code's expressions, unchanged (η taken at the start of the step)
        const eta = modele.eta(t);
        diff = 0.5 * eta * eta / (dx * dx);
        for (let i = 0; i < m - 1; i++) { lower[i] = -dt * diff; upper[i] = -dt * diff; }
        for (let i = 0; i < m; i++) { diag[i] = 1 + 2 * dt * diff + dt * (A + bt * xs[i + 1]); rhs[i] = v[i + 1]; }
      } else {
        // θ-scheme with every coefficient at its own end of the step: η(t_n) in the implicit part, η(t_{n+1}) in the
        // explicit part (η(t) = σ(t) e^{λt} moves within a step; one value for both ends would make the scheme first
        // order). σ is read at the middle of the step, i.e. on the segment the step belongs to.
        const A1 = modele.partieFixe(ts[n + 1]), b1 = modele.b(ts[n + 1]), s = modele.sigma(0.5 * (t + ts[n + 1]));
        const e0 = s * Math.exp(modele.lam * t), e1 = s * Math.exp(modele.lam * ts[n + 1]);
        diff = 0.5 * e0 * e0 / (dx * dx);
        const diff1 = 0.5 * e1 * e1 / (dx * dx);
        for (let i = 0; i < m - 1; i++) { lower[i] = -th * dt * diff; upper[i] = -th * dt * diff; }
        for (let i = 0; i < m; i++) {
          diag[i] = 1 + th * dt * (2 * diff + A + bt * xs[i + 1]);
          rhs[i] = v[i + 1] - (1 - th) * dt * ((2 * diff1 + A1 + b1 * xs[i + 1]) * v[i + 1] - diff1 * (v[i] + v[i + 2]));
        }
      }
      let gauche = v[0], droite = v[grille.nx - 1];
      let ig = 0, id = 0, aDesFuturs = false;
      for (const be of bordsEx) if (be.ex >= t - 1e-12) { aDesFuturs = true; ig = Math.max(ig, be.g); id = Math.max(id, be.d); }
      if (aDesFuturs) { gauche = Math.max(gauche, ig); droite = Math.max(droite, id); }
      rhs[0] += th * dt * diff * gauche; rhs[m - 1] += th * dt * diff * droite;
      const inner = thomas(lower, diag, upper, rhs);
      const vNew = new Float64Array(grille.nx);
      vNew[0] = gauche; vNew[grille.nx - 1] = droite; vNew.set(inner, 1);
      if (exSet.has(r12(t))) {
        const tEx = berm.exerciseTimes.find(u => r12(u) === r12(t)), ex = exercice(tEx);
        if (avecProfil) tranches.push({t: tEx, continuation: Float64Array.from(vNew), exercice: ex});
        for (let j = 0; j < vNew.length; j++) vNew[j] = Math.max(vNew[j], ex[j]);
        lisser = rannacher;
      }
      v = vNew;
    }
    if (avecProfil) {
      // at the last date the continuation value is zero: the option is exercised as soon as it is in the money
      if (!tranches.some(tr => r12(tr.t) === r12(dernier)))
        tranches.push({t: dernier, continuation: new Float64Array(grille.nx), exercice: exercice(dernier)});
      tranches.sort((a, b) => a.t - b.t);
    }
    // np.interp(0, xs, v)
    let prix;
    if (0 <= xs[0]) prix = v[0];
    else if (0 >= xs[xs.length - 1]) prix = v[v.length - 1];
    else { let k = 0; while (xs[k + 1] < 0) k++; const w = (0 - xs[k]) / (xs[k + 1] - xs[k]); prix = v[k] * (1 - w) + v[k + 1] * w; }
    return avecProfil ? {prix, xs, v, nPasTemps: ts.length - 1, tranches} : prix;
  }

  // ------------------------------------------------------------------ Longstaff-Schwartz Monte Carlo (mc_1f.py)
  function generateur(graine) { // xoshiro128** seeded by splitmix32, normals by Box-Muller
    let s = (graine >>> 0) || 1;
    const sm = () => { s = (s + 0x9e3779b9) >>> 0; let z = s; z = Math.imul(z ^ (z >>> 16), 0x85ebca6b) >>> 0; z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35) >>> 0; return (z ^ (z >>> 16)) >>> 0; };
    let a = sm(), b = sm(), c = sm(), d = sm();
    const u32 = () => { const r = Math.imul(((Math.imul(b, 5) << 7) | (Math.imul(b, 5) >>> 25)), 9) >>> 0; const t = b << 9;
      c ^= a; d ^= b; b ^= c; a ^= d; c ^= t; d = (d << 11) | (d >>> 21); return r; };
    const unif = () => ((u32() >>> 5) * 67108864 + (u32() >>> 6)) / 9007199254740992 || 1e-300;
    let reserve = null;
    const normale = () => {
      if (reserve !== null) { const r = reserve; reserve = null; return r; }
      const u1 = unif(), u2 = unif(), r = Math.sqrt(-2 * Math.log(u1)), th = 2 * Math.PI * u2;
      reserve = r * Math.sin(th); return r * Math.cos(th);
    };
    return {unif, normale};
  }
  // least squares by Householder QR (like np.linalg.lstsq at full rank)
  function moindresCarres(A, y, n, p) {
    const Q = Float64Array.from(A), z = Float64Array.from(y);
    for (let k = 0; k < p; k++) {
      let nrm = 0; for (let i = k; i < n; i++) nrm += Q[i * p + k] * Q[i * p + k];
      nrm = Math.sqrt(nrm); if (nrm === 0) continue;
      const alpha = Q[k * p + k] > 0 ? -nrm : nrm;
      const v = new Float64Array(n); for (let i = k; i < n; i++) v[i] = Q[i * p + k]; v[k] -= alpha;
      let vv = 0; for (let i = k; i < n; i++) vv += v[i] * v[i]; if (vv === 0) continue;
      for (let j = k; j < p; j++) { let s = 0; for (let i = k; i < n; i++) s += v[i] * Q[i * p + j]; s = 2 * s / vv; for (let i = k; i < n; i++) Q[i * p + j] -= s * v[i]; }
      let s = 0; for (let i = k; i < n; i++) s += v[i] * z[i]; s = 2 * s / vv; for (let i = k; i < n; i++) z[i] -= s * v[i];
    }
    const beta = new Float64Array(p);
    for (let k = p - 1; k >= 0; k--) { let s = z[k]; for (let j = k + 1; j < p; j++) s -= Q[k * p + j] * beta[j]; const dkk = Q[k * p + k]; beta[k] = Math.abs(dkk) > 1e-300 ? s / dkk : 0; }
    return beta;
  }
  // Paths of the factor at the exercise dates.
  //  "code": X drawn exactly at the exercise dates only (as in mc_1f.py);
  //  "exacte" (exact): X drawn on a fine grid, and the integral I = ∫ r ds accumulated by trapezoids (r = A(t) + b(t) X_t).
  function simulerChemins(modele, ex, N, rng, exacte, sousPasParAn) {
    const n = ex.length, x = new Float64Array(N * n);
    if (!exacte) {
      let aPrev = 0;
      for (let j = 0; j < n; j++) {
        const at = modele.alpha(ex[j]), sd = Math.sqrt(Math.max(at - aPrev, 0));
        for (let i = 0; i < N; i++) x[i * n + j] = (j > 0 ? x[i * n + j - 1] : 0) + sd * rng.normale();
        aPrev = at;
      }
      return {x, I: null};
    }
    const grille = [0], marque = [];
    let t0 = 0;
    for (let j = 0; j < n; j++) {
      const k = Math.max(1, Math.ceil((ex[j] - t0) * sousPasParAn));
      for (let q = 1; q <= k; q++) grille.push(q === k ? ex[j] : t0 + (ex[j] - t0) * q / k);
      marque.push(grille.length - 1); t0 = ex[j];
    }
    const G = grille.length, sd = new Float64Array(G), A = new Float64Array(G), B = new Float64Array(G);
    let aPrev = 0;
    for (let g = 0; g < G; g++) { const at = modele.alpha(grille[g]); sd[g] = Math.sqrt(Math.max(at - aPrev, 0)); aPrev = at; A[g] = modele.partieFixe(grille[g]); B[g] = modele.b(grille[g]); }
    const I = new Float64Array(N * n);
    for (let i = 0; i < N; i++) {
      let xc = 0, ic = 0, rPrev = A[0], j = 0;
      for (let g = 1; g < G; g++) {
        xc += sd[g] * rng.normale();
        const r = A[g] + B[g] * xc;
        ic += 0.5 * (rPrev + r) * (grille[g] - grille[g - 1]);
        rPrev = r;
        if (g === marque[j]) { x[i * n + j] = xc; I[i * n + j] = ic; j++; }
      }
    }
    return {x, I};
  }
  // Longstaff-Schwartz. Options:
  //  actualisation (discounting): "code" (mc_1f.py: short rate frozen between two dates, then e^{−r(0) τ} back to today)
  //                  or "exacte" (e^{−∫r} simulated path by path);
  //  degre: degree of the regression polynomial (the code: 2, basis [1, X, X²]);
  //  independants: regression on one set of paths, pricing on another (removes the in-sample high bias);
  //  detail: returns, for each date, a sample (X, discounted continuation value) and the regression.
  function prixLSMC(modele, berm, nChemins = 20000, graine = 1234, opts = {}) {
    const ex = berm.exerciseTimes.slice(), n = ex.length, N = nChemins, rng = generateur(graine);
    const exacte = opts.actualisation === "exacte", degre = Math.max(1, Math.min(6, Math.round(opts.degre || 2))), p = degre + 1;
    const sousPas = Math.max(1, opts.sousPasParAn || 52);
    const S = simulerChemins(modele, ex, N, rng, exacte, sousPas);
    const col = (X, j, M) => { const c = new Float64Array(M); for (let i = 0; i < M; i++) c[i] = X[i * n + j]; return c; };
    const ecartX = ex.map(t => Math.sqrt(Math.max(modele.alpha(t), 1e-300)));
    const base = (xx, j, k) => Math.pow(xx / ecartX[j], k); // normalised X: same fit, better conditioned
    const r0 = modele.tauxCourt(0, 0);
    // discounting from date j to exercise date k, path i
    const actu = (Sx, i, j, k, xx) => {
      if (exacte) return Math.exp(-(Sx.I[i * n + k] - Sx.I[i * n + j]));
      return Math.exp(-(modele.partieFixe(ex[j]) + modele.b(ex[j]) * xx) * (ex[k] - ex[j]));
    };
    const actu0 = (Sx, i, k) => exacte ? Math.exp(-Sx.I[i * n + k]) : Math.exp(-r0 * ex[k]);
    // backward pass on the training set
    const flux = new Float64Array(N), kEx = new Int32Array(N).fill(n - 1);
    const valeurs = X => ex.map((t, j) => { const v = modele.valeursSwapGrille(berm.sousJacente(t), col(X, j, N)); for (let i = 0; i < N; i++) v[i] = Math.max(v[i], 0); return v; });
    const vals = valeurs(S.x);
    for (let i = 0; i < N; i++) flux[i] = vals[n - 1][i];
    const betas = new Array(n).fill(null), details = [];
    for (let j = n - 2; j >= 0; j--) {
      const xj = col(S.x, j, N), vEx = vals[j], itm = [];
      for (let i = 0; i < N; i++) if (vEx[i] > 1e-14) itm.push(i);
      if (!itm.length) continue;
      const m = itm.length, M = new Float64Array(m * p), y = new Float64Array(m);
      itm.forEach((i, q) => { for (let k = 0; k < p; k++) M[q * p + k] = base(xj[i], j, k); y[q] = flux[i] * actu(S, i, j, kEx[i], xj[i]); });
      const beta = moindresCarres(M, y, m, p);
      betas[j] = beta;
      if (opts.detail) {
        const pas = Math.max(1, Math.floor(m / 1200)), dx = [], dy = [], dv = [];
        for (let q = 0; q < m; q += pas) { dx.push(xj[itm[q]]); dy.push(y[q]); dv.push(vEx[itm[q]]); }
        details[j] = {t: ex[j], x: dx, y: dy, v: dv, beta: Array.from(beta), sd: ecartX[j], nItm: m};
      }
      itm.forEach(i => { let c = 0; for (let k = 0; k < p; k++) c += beta[k] * base(xj[i], j, k);
        if (vEx[i] > c) { flux[i] = vEx[i]; kEx[i] = j; } });
    }
    let Sp = S, fluxP = flux, kP = kEx;
    if (opts.independants) {
      // new paths: the exercise rule learnt above is applied forwards
      Sp = simulerChemins(modele, ex, N, rng, exacte, sousPas);
      const valsP = valeurs(Sp.x);
      fluxP = new Float64Array(N); kP = new Int32Array(N).fill(n - 1);
      for (let i = 0; i < N; i++) {
        let fait = false;
        for (let j = 0; j < n - 1 && !fait; j++) {
          const v = valsP[j][i]; if (!(v > 1e-14) || !betas[j]) continue;
          const xx = Sp.x[i * n + j]; let c = 0; for (let k = 0; k < p; k++) c += betas[j][k] * base(xx, j, k);
          if (v > c) { fluxP[i] = v; kP[i] = j; fait = true; }
        }
        if (!fait) fluxP[i] = valsP[n - 1][i];
      }
    }
    let s = 0, s2 = 0;
    const nbExerces = new Array(n).fill(0);
    for (let i = 0; i < N; i++) { const val = fluxP[i] * actu0(Sp, i, kP[i]); s += val; s2 += val * val; if (fluxP[i] > 0) nbExerces[kP[i]]++; }
    const moy = s / N, var_ = N > 1 ? (s2 - N * moy * moy) / (N - 1) : 0;
    // martingale check: Monte Carlo prices of the zero-coupon bonds P(0, T_k) with the chosen discounting
    const martingale = ex.map((t, k) => {
      let a = 0, a2 = 0;
      for (let i = 0; i < N; i++) { const d = actu0(Sp, i, k); a += d; a2 += d * d; }
      const mm = a / N;
      return {t, mc: mm, se: Math.sqrt(Math.max(a2 / N - mm * mm, 0) / N), courbe: modele.courbe.discount(t)};
    });
    return {prix: moy, erreur: Math.sqrt(Math.max(var_, 0) / N), nbExerces, details, martingale};
  }
  // implied normal (Bachelier) volatility of a swaption price
  function volNImplicite(courbe, sw, prix) {
    const f = v => prixBachelier(courbe, sw, v) - prix;
    if (f(1e-9) > 0) return NaN;
    try { return brent(f, 1e-9, 1, 1e-14, 200); } catch (e) { return NaN; }
  }
  // short-rate fan: quantiles of r_t = A(t) + b(t) X_t, X_t ~ N(0, α(t))
  function eventailTaux(modele, tMax, zs, nPts = 120) {
    const ts = [], q = zs.map(() => []);
    for (let k = 0; k <= nPts; k++) {
      const t = tMax * k / nPts, A = modele.partieFixe(t), b = modele.b(t), sd = Math.sqrt(Math.max(modele.alpha(t), 0));
      ts.push(t); zs.forEach((z, i) => q[i].push(A + b * z * sd));
    }
    return {ts, q};
  }

  const API = {erf, erfc, phi, pdf, brent, Courbe, Swaption, Bermudeenne, annuite, tauxSwapForward, prixBachelier,
    Modele, bootstrap, grilleEcartsTypes, grilleTemps, prixEDP, prixLSMC, simulerChemins, volNImplicite, eventailTaux, generateur, moindresCarres, linspace};
  if (typeof module !== "undefined" && module.exports) module.exports = API; else racine.HJM = API;
})(typeof self !== "undefined" ? self : this);
