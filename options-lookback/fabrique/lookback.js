/* lookback.js — portage JavaScript de la DLL unifiée du projet (code/unified/unified_dll.cpp), étendu.

   priceLookback(method, ttm, optionType, spot, rate, vol, paths, steps, seed, withStdErrors, useBGK, options)
   renvoie les 48 valeurs de la fonction C++ price_lookback : pour le Monte-Carlo puis la formule fermée, 6 grandeurs
   (prix, Δ, Γ, Θ, ρ, vega) × (estimation, erreur standard, IC 95 % bas, IC 95 % haut).
     method 1 : moteur Barnicha (pont brownien exact en un pas)   — code/Barnicha/src/lookback_dll.cpp
     method 0 : moteur Chivet (schéma en temps, correction BGK)  — code/Chivet/src/LookbackDLL.cpp
   Sans « options », le calcul est celui du C++ (voir les écarts ci-dessous).

   Générateur : std::mt19937 et les lois uniform_real_distribution / normal_distribution de libstdc++ reproduits à
   l'identique : avec la même graine, le moteur Barnicha donne les nombres du C++ compilé avec g++ (vérifié à 1e-15).

   Écarts volontaires avec le C++ (défauts relevés à la relecture, signalés dans la page) :
   - Chivet : la graine saisie est utilisée (le C++ prend l'horloge) et la loi normale repart de zéro à chaque prix
     (le C++ garde d'un prix à l'autre la normale mise en réserve : les tirages communs des grecques sont décalés) ;
   - Chivet : l'erreur standard du rho est divisée par 2·bump·r comme l'estimation (le C++ oublie le facteur r) ;
   - Chivet, formule fermée : Γ divisé par h² (le C++ divise par (2h)² : Γ quatre fois trop petit).

   Extensions (options) : dividende continu q ; extremum déjà observé (option en cours de vie) ; correction de la
   surveillance discrète au choix (aucune, BGK, ou pont brownien exact sur chaque pas) ; tailles des bumps ;
   variable antithétique ; Θ centré ; coefficient β de BGK. */
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
  // std::generate_canonical<double, 53>(mt19937) de libstdc++ : deux tirages de 32 bits
  const DEUX32 = 4294967296, UN_MOINS = 1 - Math.pow(2, -53);
  function canonique(g) {
    const r = (g() + g() * DEUX32) / (DEUX32 * DEUX32);
    return r >= 1 ? UN_MOINS : r;
  }
  // std::normal_distribution<double> de libstdc++ : méthode polaire, une valeur gardée en réserve
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

  // ------------------------------------------------------------------ loi normale
  function erfc(x) { // ~1e-15 : série (|x| < 3) ou fraction continue de Lentz
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

  // ------------------------------------------------------------------ formules fermées
  // Barnicha (analytic_lookback.cpp) : émission (extremum = spot), sans dividende, limite exacte en r → 0
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
  // Chivet (AnalyticLookback.cpp) : même prix écrit autrement ; non défini pour r = 0 (division par 2r)
  function prixFermeChivet(put, S0, T, sig, r) {
    const dp = (Math.log(1) + (r + 0.5 * sig * sig) * T) / (sig * Math.sqrt(T)), dm = (Math.log(1) + (r - 0.5 * sig * sig) * T) / (sig * Math.sqrt(T));
    if (!put) return S0 * N(dp) - S0 * Math.exp(-r * T) * N(dm) + Math.exp(-r * T) * (sig * sig * S0 / (2 * r)) * N(dm) - N(-dp) * S0 * (sig * sig / (2 * r));
    return S0 * Math.exp(-r * T) * N(-dm) + S0 * N(dp) * (1 + (sig * sig) / (2 * r)) - Math.exp(-r * T) * S0 * (sig * sig / (2 * r)) * N(-dm) - S0;
  }
  // Goldman-Sosin-Gatto général : portage b = r − q, extremum déjà observé m (min pour le call, max pour le put)
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
    if (Math.abs(b) < 1e-6) return 0.5 * (gsg(put, S, m, r, 1e-6, sig, T) + gsg(put, S, m, r, -1e-6, sig, T)); // limite b → 0
    return gsg(put, S, m, r, b, sig, T);
  }
  // prix fermé utilisé par un moteur : celui de son code tant qu'on reste dans son cadre (émission, q = 0)
  function prixFerme(moteur, m) {
    const ext = m.ext == null ? m.S : m.ext;
    if (m.q === 0 && ext === m.S) return moteur === 1 ? prixFermeBarnicha(m.put, m.S, m.r, m.sig, m.T) : prixFermeChivet(m.put, m.S, m.T, m.sig, m.r);
    return prixFermeGeneral(m.put, m.S, ext, m.r, m.q, m.sig, m.T);
  }

  // ------------------------------------------------------------------ statistiques (pricing_engine.cpp)
  function stats(xs) {
    const k = xs.length;
    if (!k) return [0, 0, 0, 0];
    let s = 0; for (let i = 0; i < k; i++) s += xs[i];
    const moy = s / k; let v = 0; for (let i = 0; i < k; i++) { const d = xs[i] - moy; v += d * d; }
    const se = k > 1 ? Math.sqrt(v / (k - 1)) / Math.sqrt(k) : 0;
    return [moy, se, moy - 1.96 * se, moy + 1.96 * se];
  }
  const avec = (m, o) => Object.assign({}, m, o);

  // ================================================================== moteur Barnicha
  // tirages (z, u_max, u_min) : OneStepBridgeSampler::generate, dans l'ordre du C++
  function tiragesPont(nChemins, graine) {
    const g = MT19937(graine), loi = LoiNormale();
    const z = new Float64Array(nChemins), umax = new Float64Array(nChemins), umin = new Float64Array(nChemins);
    for (let i = 0; i < nChemins; i++) {
      z[i] = loi(g);
      let u = canonique(g); umax[i] = u === 0 ? 1e-16 : u;
      u = canonique(g); umin[i] = u === 0 ? 1e-16 : u;
    }
    // logarithmes des uniformes calculés une fois (ils resservent pour chaque marché décalé des grecques)
    const lmax = umax.map(Math.log), lmin = umin.map(Math.log);
    return {z, umax, umin, lmax, lmin, n: nChemins};
  }
  // gains actualisés par chemin (monte_carlo_engine.cpp + lookback_payoff_evaluator.cpp)
  function gainsPont(t, m, antithetique) {
    const out = new Float64Array(t.n), x0 = Math.log(m.S), derive = (m.r - m.q - 0.5 * m.sig * m.sig) * m.T, vol = m.sig * Math.sqrt(m.T);
    const disc = Math.exp(-m.r * m.T), k = -0.5 * m.sig * m.sig * m.T, put = m.put;
    const ext = m.ext == null ? null : m.ext; // extremum déjà observé (sinon : le spot)
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
  // détail chemin par chemin (sans antithétique) pour les graphes descriptifs :
  // log(S_T/S), log(extremum du pont / S) et gain actualisé
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
  // fonction de répartition exacte de l'extremum de X_t = μt + σW_t sur [0, T] (en log-rendement)
  //   max : P(M ≤ y) = N((y − μT)/(σ√T)) − e^{2μy/σ²} N((−y − μT)/(σ√T)), y ≥ 0
  //   min : P(m ≥ y) = N((−y + μT)/(σ√T)) − e^{2μy/σ²} N((y + μT)/(σ√T)), y ≤ 0
  function repartitionExtremum(put, y, mu, sig, T) {
    const s = sig * Math.sqrt(T);
    if (put) return y <= 0 ? 0 : N((y - mu * T) / s) - Math.exp(2 * mu * y / (sig * sig)) * N((-y - mu * T) / s);
    return y >= 0 ? 0 : N((-y + mu * T) / s) - Math.exp(2 * mu * y / (sig * sig)) * N((y + mu * T) / s);
  }
  // BumpAndRevalueGreeks::compute (bumps par défaut de bump_and_revalue_greeks.h : 1e-4)
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

  // ================================================================== moteur Chivet
  const BGK_BETA1 = 0.5826; // −ζ(1/2)/√(2π)
  // computeRawPrice : mêmes tirages (graine fixe) pour chaque marché décalé.
  // correction : "aucune", "bgk" (C++), ou "pont" : extremum exact du pont brownien sur chaque pas (sans biais)
  // Plusieurs marchés décalés évalués sur les MÊMES tirages : chaque chemin tire ses normales (et ses uniformes pour le
  // pont) une seule fois, dans l'ordre du C++, puis les rejoue pour chaque marché. Résultat identique à une simulation
  // séparée par marché réamorcée à la même graine (ce que fait le C++), en une seule passe du générateur.
  // Calcul en log : l'extremum de S est l'exponentielle de l'extremum de log S.
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
          // extremum d'un pont brownien entre x et x1 : P(max ≤ y) = 1 − exp(−2 (y − x)(y − x1) / (σ² Δt))
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
    // les 9 marchés de Pricer::compute, sur les mêmes tirages
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
    // formule fermée et grecques par bumps relatifs de 1e-4 (AnalyticLookback.cpp), Γ corrigé
    const F = x => prixFerme(0, avec(m, x)), V0 = F({}), Vp = F({S: m.S * 1.0001}), Vm = F({S: m.S * 0.9999});
    const ex = v => [v, 0, v, v];
    out.push(...ex(V0), ...ex((Vp - Vm) / (m.S * 0.0002)), ...ex((Vp - 2 * V0 + Vm) / Math.pow(m.S * 0.0001, 2)),
      ...ex((F({T: m.T * 1.0001}) - F({T: m.T * 0.9999})) / (m.T * 0.0002)),
      ...ex((F({r: m.r * 1.0001}) - F({r: m.r * 0.9999})) / (m.r * 0.0002)),
      ...ex((F({sig: m.sig * 1.0001}) - F({sig: m.sig * 0.9999})) / (m.sig * 0.0002)));
    return out;
  }

  // ================================================================== point d'entrée (unified_dll.cpp)
  function priceLookback(method, ttm, optionType, spot, rate, vol, paths, steps, seed, withStdErrors, useBGK, options) {
    const o = options || {};
    if (!(ttm > 0)) throw new Error("maturité nulle ou négative (code -2 de la DLL)");
    if (!(paths > 0)) throw new Error("nombre de chemins nul (code -3 de la DLL)");
    const put = optionType === 1, ext = o.extremum == null ? null : +o.extremum;
    if (ext !== null && (put ? ext < spot : ext > spot)) throw new Error(put ? "le maximum déjà observé doit être ≥ au spot" : "le minimum déjà observé doit être ≤ au spot");
    const m = {put, S: spot, r: rate, q: o.q || 0, sig: vol, T: ttm, ext};
    if (method === 1) return barnicha(m, Math.round(paths), seed >>> 0, o);
    if (method === 0) return chivet(m, Math.round(paths), Math.max(1, Math.round(steps)), seed >>> 0, !!withStdErrors, o.correction || (useBGK ? "bgk" : "aucune"), o);
    throw new Error("méthode inconnue (code -999 de la DLL)");
  }

  const API = {priceLookback, MT19937, canonique, LoiNormale, tiragesPont, gainsPont, gainsChivet, prixFermeBarnicha, prixFermeChivet,
    prixFermeGeneral, detailPont, repartitionExtremum, gainsChivetMulti, BGK_BETA1};
  if (typeof module !== "undefined" && module.exports) module.exports = API; else racine.LB = API;
})(typeof self !== "undefined" ? self : this);
