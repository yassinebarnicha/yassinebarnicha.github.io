/* graphe.js — petits graphiques sur canevas, sans bibliothèque.
   tracer(conteneur, {series, xlab, ylab, xfmt, yfmt, xmin, xmax, ymin, ymax, zero})
     series : [{nom, x:[], y:[], type:'ligne'|'points'|'marches'|'barres', couleur:'c1'..'c8', tirets, epaisseur, bande:[ybas,yhaut]}]
   carte(conteneur, {nx, ny, x0, x1, y0, y1, valeur(i,j)->0..1 ou NaN, couleurs:[bas,haut], xlab, ylab, points:[{x,y,couleur}]})
   Le graphe se redessine seul quand son conteneur change de taille ou que le thème change. */
(function () {
  "use strict";
  const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const couleur = c => (c && /^c\d$/.test(c)) ? css("--" + c) : (c && c.startsWith("--") ? css(c) : c || css("--c1"));
  const fmtDefaut = v => {
    if (!isFinite(v)) return "—";
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e6)) return v.toExponential(2);
    if (a >= 100) return v.toFixed(a >= 1e4 ? 0 : 1);
    return +v.toPrecision(4) + "";
  };
  // graduations « rondes »
  function graduations(a, b, n) {
    if (!(b > a)) { const m = a || 1; a = m - Math.abs(m) * 0.1 - 1e-12; b = m + Math.abs(m) * 0.1 + 1e-12; }
    const brut = (b - a) / Math.max(n, 1), p = Math.pow(10, Math.floor(Math.log10(brut))), r = brut / p;
    const pas = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * p;
    const t = [], d = Math.ceil(a / pas - 1e-9) * pas;
    for (let v = d; v <= b + pas * 1e-9; v += pas) t.push(Math.abs(v) < pas * 1e-9 ? 0 : v);
    return t;
  }
  const tous = new Set();
  function surveiller(el, dessin) {
    el.__dessin = dessin;
    if (!el.__obs) {
      el.__obs = new ResizeObserver(() => el.__dessin && el.__dessin());
      el.__obs.observe(el);
      tous.add(el);
    }
    dessin();
  }
  const mq = window.matchMedia("(prefers-color-scheme: dark)");
  const redessinerTout = () => tous.forEach(el => el.isConnected ? el.__dessin && el.__dessin() : tous.delete(el));
  if (mq.addEventListener) mq.addEventListener("change", redessinerTout);
  new MutationObserver(redessinerTout).observe(document.documentElement, {attributes: true, attributeFilter: ["data-theme"]});

  function preparer(el) {
    let cv = el.querySelector("canvas"), bulle = el.querySelector(".bulle");
    if (!cv) { cv = document.createElement("canvas"); el.appendChild(cv); }
    if (!bulle) { bulle = document.createElement("div"); bulle.className = "bulle"; el.appendChild(bulle); }
    const r = el.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
    const W = Math.max(50, el.clientWidth), H = Math.max(50, el.clientHeight);
    cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
    const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0);
    return {cv, bulle, g, W, H, r};
  }
  function legende(el, series) {
    let lg = el.nextElementSibling;
    if (!lg || !lg.classList || !lg.classList.contains("legende")) {
      lg = document.createElement("div"); lg.className = "legende"; el.after(lg);
    }
    const vis = series.filter(s => s.nom && !s.sansLegende);
    lg.innerHTML = vis.length < 2 ? "" : vis.map(s => {
      const c = couleur(s.couleur), pt = s.type === "points";
      return `<span><i class="${pt ? "pt" : ""}" style="background:${c}${s.tirets ? ";background:repeating-linear-gradient(90deg," + c + " 0 4px,transparent 4px 7px)" : ""}"></i>${s.nom}</span>`;
    }).join("");
  }

  function tracer(el, o) {
    o = Object.assign({series: [], xfmt: fmtDefaut, yfmt: fmtDefaut}, o);
    legende(el, o.series);
    surveiller(el, () => {
      const {cv, bulle, g, W, H} = preparer(el);
      const encre = css("--ink"), encre2 = css("--ink2"), grille = css("--grid"), fond = css("--panel");
      const mono = css("--mono") || "monospace";
      // domaine
      let xa = Infinity, xb = -Infinity, ya = Infinity, yb = -Infinity;
      for (const s of o.series) for (let i = 0; i < s.x.length; i++) {
        const x = s.x[i], y = s.y[i];
        if (!isFinite(x) || !isFinite(y)) continue;
        if (x < xa) xa = x; if (x > xb) xb = x;
        if (y < ya) ya = y; if (y > yb) yb = y;
        if (s.bande) { const lo = s.bande[0][i], hi = s.bande[1][i]; if (isFinite(lo) && lo < ya) ya = lo; if (isFinite(hi) && hi > yb) yb = hi; }
      }
      if (!isFinite(xa)) { g.fillStyle = encre2; g.font = "12px " + mono; g.fillText("Aucune donnée", 12, 24); return; }
      if (o.series.some(s => s.type === "barres")) {
        ya = Math.min(ya, 0); yb = Math.max(yb, 0);
        // une demi-largeur de barre de marge de chaque côté
        let pasMin = Infinity;
        for (const s of o.series) if (s.type === "barres") for (let i = 1; i < s.x.length; i++) pasMin = Math.min(pasMin, Math.abs(s.x[i] - s.x[i - 1]));
        if (!isFinite(pasMin)) pasMin = 1;
        xa -= pasMin / 2; xb += pasMin / 2;
      }
      if (o.zero) { ya = Math.min(ya, 0); yb = Math.max(yb, 0); }
      if (o.xmin != null) xa = o.xmin; if (o.xmax != null) xb = o.xmax;
      if (o.ymin != null) ya = o.ymin; if (o.ymax != null) yb = o.ymax;
      if (yb - ya < 1e-300) { const m = Math.abs(ya) || 1; ya -= m * 0.05; yb += m * 0.05; }
      else { const pad = (yb - ya) * 0.06; if (o.ymin == null) ya -= pad; if (o.ymax == null) yb += pad; }
      if (xb - xa < 1e-300) { xa -= 0.5; xb += 0.5; }
      const ty = graduations(ya, yb, Math.max(3, Math.round(H / 60))), tx = graduations(xa, xb, Math.max(3, Math.round(W / 95)));
      g.font = "11px " + mono;
      const lyMax = Math.max(...ty.map(v => g.measureText(o.yfmt(v)).width));
      const gauche = Math.ceil(lyMax) + 14 + (o.ylab ? 14 : 0), droite = 12, haut = 10, bas = 30 + (o.xlab ? 14 : 0);
      const pw = W - gauche - droite, ph = H - haut - bas;
      const X = x => gauche + (x - xa) / (xb - xa) * pw, Y = y => haut + (1 - (y - ya) / (yb - ya)) * ph;
      g.clearRect(0, 0, W, H);
      // grille
      g.lineWidth = 1; g.strokeStyle = grille; g.fillStyle = encre2; g.textAlign = "right"; g.textBaseline = "middle";
      for (const v of ty) { const y = Math.round(Y(v)) + .5; if (y < haut - 1 || y > haut + ph + 1) continue;
        g.beginPath(); g.moveTo(gauche, y); g.lineTo(gauche + pw, y); g.stroke(); g.fillText(o.yfmt(v), gauche - 6, y); }
      g.textAlign = "center"; g.textBaseline = "top";
      for (const v of tx) { const x = Math.round(X(v)) + .5; if (x < gauche - 1 || x > gauche + pw + 1) continue;
        g.beginPath(); g.moveTo(x, haut); g.lineTo(x, haut + ph); g.stroke(); g.fillText(o.xfmt(v), x, haut + ph + 6); }
      g.strokeStyle = encre2; g.beginPath(); g.moveTo(gauche, haut + ph + .5); g.lineTo(gauche + pw, haut + ph + .5); g.stroke();
      if (ya < 0 && yb > 0) { g.strokeStyle = encre2; g.globalAlpha = .5; g.beginPath(); g.moveTo(gauche, Math.round(Y(0)) + .5); g.lineTo(gauche + pw, Math.round(Y(0)) + .5); g.stroke(); g.globalAlpha = 1; }
      if (o.xlab) { g.fillStyle = encre2; g.textAlign = "center"; g.fillText(o.xlab, gauche + pw / 2, H - 15); }
      if (o.ylab) { g.save(); g.translate(11, haut + ph / 2); g.rotate(-Math.PI / 2); g.textAlign = "center"; g.textBaseline = "middle"; g.fillStyle = encre2; g.fillText(o.ylab, 0, 0); g.restore(); }
      // séries
      g.save(); g.beginPath(); g.rect(gauche - 4, haut - 4, pw + 8, ph + 8); g.clip();
      const nb = o.series.filter(s => s.type === "barres").length;
      let ib = 0;
      for (const s of o.series) {
        const c = couleur(s.couleur);
        if (s.bande) {
          g.fillStyle = c; g.globalAlpha = .14; g.beginPath();
          let deb = true;
          for (let i = 0; i < s.x.length; i++) { const lo = s.bande[0][i]; if (!isFinite(lo)) continue; deb ? g.moveTo(X(s.x[i]), Y(lo)) : g.lineTo(X(s.x[i]), Y(lo)); deb = false; }
          for (let i = s.x.length - 1; i >= 0; i--) { const hi = s.bande[1][i]; if (!isFinite(hi)) continue; g.lineTo(X(s.x[i]), Y(hi)); }
          g.closePath(); g.fill(); g.globalAlpha = 1;
        }
        if (s.type === "barres") {
          const n = s.x.length, pasPx = n > 1 ? Math.abs(X(s.x[1]) - X(s.x[0])) : pw / 4, larg = Math.max(2, pasPx * 0.8 / nb - 2);
          g.fillStyle = c;
          for (let i = 0; i < n; i++) { if (!isFinite(s.y[i])) continue; const x = X(s.x[i]) - pasPx * 0.4 + ib * (larg + 2), y0 = Y(0), y1 = Y(s.y[i]);
            g.fillRect(x, Math.min(y0, y1), larg, Math.max(1, Math.abs(y1 - y0))); }
          ib++; continue;
        }
        if (s.type !== "points") {
          g.strokeStyle = c; g.lineWidth = s.epaisseur || 2; g.setLineDash(s.tirets ? [5, 4] : []); g.lineJoin = "round"; g.beginPath();
          let deb = true, py = 0;
          for (let i = 0; i < s.x.length; i++) {
            if (!isFinite(s.y[i]) || !isFinite(s.x[i])) { deb = true; continue; }
            const x = X(s.x[i]), y = Y(s.y[i]);
            if (deb) g.moveTo(x, y); else if (s.type === "marches") { g.lineTo(x, py); g.lineTo(x, y); } else g.lineTo(x, y);
            deb = false; py = y;
          }
          g.stroke(); g.setLineDash([]);
        }
        if (s.type === "points" || s.marqueurs) {
          const rayon = s.rayon || (s.type === "points" ? 3.5 : 3);
          for (let i = 0; i < s.x.length; i++) { if (!isFinite(s.y[i])) continue;
            g.beginPath(); g.arc(X(s.x[i]), Y(s.y[i]), rayon, 0, 2 * Math.PI);
            g.fillStyle = s.creux ? fond : c; g.fill(); g.lineWidth = 1.5; g.strokeStyle = s.creux ? c : fond; g.stroke(); }
        }
      }
      g.restore();
      // survol : point le plus proche
      el.onpointermove = ev => {
        const rc = cv.getBoundingClientRect(), ech = rc.width ? W / rc.width : 1;
        const mx = (ev.clientX - rc.left) * ech, my = (ev.clientY - rc.top) * ech;
        let best = null, bd = 28 * 28;
        for (const s of o.series) for (let i = 0; i < s.x.length; i++) {
          if (!isFinite(s.y[i]) || !isFinite(s.x[i])) continue;
          const dx = X(s.x[i]) - mx, dy = Y(s.y[i]) - my, d = dx * dx + dy * dy;
          if (d < bd) { bd = d; best = {s, i}; }
        }
        if (!best) { bulle.style.display = "none"; return; }
        const {s, i} = best;
        const ex = s.etiquettes ? s.etiquettes[i] : (o.xlab ? o.xlab + " = " : "x = ") + o.xfmt(s.x[i]);
        bulle.innerHTML = (s.nom ? `<b>${s.nom}</b><br>` : "") + ex + "<br>" + (o.ylab ? o.ylab + " = " : "y = ") + (o.yfmtBulle || o.yfmt)(s.y[i]);
        bulle.style.display = "block";
        const bx = X(s.x[i]) / ech, by = Y(s.y[i]) / ech, bw = bulle.offsetWidth;
        bulle.style.left = Math.min(Math.max(0, bx + 12), el.clientWidth - bw - 2) + "px";
        bulle.style.top = Math.max(0, by - bulle.offsetHeight - 10) + "px";
      };
      el.onpointerleave = () => { bulle.style.display = "none"; };
    });
  }

  // carte 2D (zones d'exercice, etc.)
  function carte(el, o) {
    legende(el, []);
    surveiller(el, () => {
      const {g, W, H, bulle, cv} = preparer(el);
      const encre2 = css("--ink2"), mono = css("--mono") || "monospace", fond = css("--panel");
      g.clearRect(0, 0, W, H); g.font = "11px " + mono;
      const gauche = 54, bas = 40, haut = 8, droite = 10;
      const cote = Math.min(W - gauche - droite, H - haut - bas), pw = cote, ph = cote;
      const X = x => gauche + (x - o.x0) / (o.x1 - o.x0) * pw, Y = y => haut + (1 - (y - o.y0) / (o.y1 - o.y0)) * ph;
      const [ca, cb] = (o.couleurs || ["c2", "c5"]).map(couleur);
      const dx = pw / o.nx, dy = ph / o.ny;
      for (let i = 0; i < o.nx; i++) for (let j = 0; j < o.ny; j++) {
        const v = o.valeur(i, j); if (!isFinite(v)) continue;
        g.fillStyle = v >= 0.5 ? cb : ca; g.globalAlpha = 0.25 + 0.55 * Math.abs(v - 0.5) * 2;
        g.fillRect(gauche + i * dx, haut + ph - (j + 1) * dy, Math.ceil(dx), Math.ceil(dy));
      }
      g.globalAlpha = 1;
      for (const p of (o.points || [])) { g.fillStyle = couleur(p.couleur || "--ink"); g.fillRect(X(p.x) - 1, Y(p.y) - 1, 2, 2); }
      g.fillStyle = encre2; g.strokeStyle = encre2; g.textAlign = "center"; g.textBaseline = "top";
      for (const v of graduations(o.x0, o.x1, 5)) g.fillText(+v.toPrecision(4) + "", X(v), haut + ph + 6);
      g.textAlign = "right"; g.textBaseline = "middle";
      for (const v of graduations(o.y0, o.y1, 5)) g.fillText(+v.toPrecision(4) + "", gauche - 6, Y(v));
      g.strokeRect(gauche + .5, haut + .5, pw, ph);
      g.textAlign = "center"; g.textBaseline = "top";
      if (o.xlab) g.fillText(o.xlab, gauche + pw / 2, H - 15);
      if (o.ylab) { g.save(); g.translate(12, haut + ph / 2); g.rotate(-Math.PI / 2); g.textBaseline = "middle"; g.fillText(o.ylab, 0, 0); g.restore(); }
      el.onpointermove = ev => {
        const rc = cv.getBoundingClientRect(), ech = rc.width ? W / rc.width : 1;
        const mx = (ev.clientX - rc.left) * ech, my = (ev.clientY - rc.top) * ech;
        if (mx < gauche || mx > gauche + pw || my < haut || my > haut + ph) { bulle.style.display = "none"; return; }
        const x = o.x0 + (mx - gauche) / pw * (o.x1 - o.x0), y = o.y0 + (1 - (my - haut) / ph) * (o.y1 - o.y0);
        const i = Math.min(o.nx - 1, Math.floor((mx - gauche) / dx)), j = Math.min(o.ny - 1, Math.floor((haut + ph - my) / dy));
        bulle.innerHTML = `${o.xlab || "x"} = ${x.toFixed(1)}<br>${o.ylab || "y"} = ${y.toFixed(1)}<br>${o.etiquette ? o.etiquette(o.valeur(i, j)) : o.valeur(i, j)}`;
        bulle.style.display = "block";
        bulle.style.left = Math.min(mx / ech + 12, el.clientWidth - bulle.offsetWidth - 2) + "px";
        bulle.style.top = Math.max(0, my / ech - bulle.offsetHeight - 10) + "px";
      };
      el.onpointerleave = () => { bulle.style.display = "none"; };
    });
  }
  window.Graphe = {tracer, carte, fmt: fmtDefaut};
})();
