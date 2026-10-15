/* ui.js — onglets (avec l'adresse #onglet), lecteur du code source, utilitaires de formulaire. */
(function () {
  "use strict";
  // ---------- onglets ----------
  function onglets() {
    const barre = document.querySelector(".onglets");
    if (!barre) return;
    const boutons = [...barre.querySelectorAll("button[data-vue]")];
    const montrer = (id, pousser) => {
      boutons.forEach(b => b.setAttribute("aria-selected", String(b.dataset.vue === id)));
      document.querySelectorAll(".vue").forEach(v => { v.hidden = v.id !== id; });
      if (pousser) history.replaceState(null, "", "#" + id);
      document.dispatchEvent(new CustomEvent("vue", {detail: id}));
    };
    boutons.forEach(b => b.addEventListener("click", () => montrer(b.dataset.vue, true)));
    const depuisAdresse = () => {
      const id = decodeURIComponent(location.hash.slice(1));
      montrer(boutons.some(b => b.dataset.vue === id) ? id : boutons[0].dataset.vue, false);
    };
    window.addEventListener("hashchange", depuisAdresse);
    depuisAdresse();
  }

  // ---------- lecteur de code ----------
  const MOTS = {
    py: "False None True and as assert break class continue def del elif else except finally for from global if import in is lambda not or pass raise return try while with yield self",
    cpp: "auto bool break case catch class const constexpr continue default delete do double else enum explicit extern false float for if inline int long namespace new nullptr operator private protected public return short signed sizeof static std struct switch template this throw true try typedef typename unsigned using virtual void while include define pragma once",
    js: "async await break case catch class const continue default delete do else export false for function if import in instanceof let new null return switch this throw true try typeof undefined var while",
    vb: "As Boolean ByVal ByRef Call Dim Do Double Else ElseIf End Exit For Function GoTo If In Long Loop Next Not Nothing On Private Public Set String Sub Then To Variant While With",
  };
  const echap = s => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  function colorer(src, lang) {
    const mots = new Set((MOTS[lang] || "").split(" "));
    const com = lang === "py" ? /#.*/y : lang === "vb" ? /'.*/y : /\/\/.*|\/\*[\s\S]*?\*\//y;
    let out = "", i = 0;
    while (i < src.length) {
      com.lastIndex = i;
      let m = com.exec(src);
      if (m) { out += `<span class="c">${echap(m[0])}</span>`; i += m[0].length; continue; }
      const ch = src[i];
      if (lang === "py" && (src.startsWith('"""', i) || src.startsWith("'''", i))) {
        const q = src.slice(i, i + 3), j = src.indexOf(q, i + 3), fin = j < 0 ? src.length : j + 3;
        out += `<span class="s">${echap(src.slice(i, fin))}</span>`; i = fin; continue;
      }
      if (ch === '"' || (ch === "'" && lang !== "vb") || (ch === "`" && lang === "js")) {
        let j = i + 1;
        while (j < src.length && src[j] !== ch && src[j] !== "\n") j += src[j] === "\\" ? 2 : 1;
        out += `<span class="s">${echap(src.slice(i, j + 1))}</span>`; i = j + 1; continue;
      }
      const w = /[A-Za-z_]\w*/y; w.lastIndex = i; m = w.exec(src);
      if (m) { out += mots.has(m[0]) ? `<span class="k">${m[0]}</span>` : m[0]; i += m[0].length; continue; }
      const n = /\d[\d_.eE]*/y; n.lastIndex = i; m = n.exec(src);
      if (m) { out += `<span class="n">${m[0]}</span>`; i += m[0].length; continue; }
      out += echap(ch); i++;
    }
    return out.split("\n").map(l => `<span class="l">${l || " "}</span>`).join("");
  }
  function lecteurCode(conteneur, fichiers, depot) {
    // fichiers : [{chemin, groupe, lang, texte, lien?}]
    conteneur.innerHTML = `<div class="code-espace"><nav class="code-liste" aria-label="Fichiers"></nav>
      <div class="code-vue"><div class="code-tete"><span class="nom"></span><span class="badge lg"></span>
      <span class="act"><button type="button" class="copier">Copier</button><a class="gh" target="_blank" rel="noopener">Ouvrir sur GitHub ↗</a></span></div>
      <pre class="code" tabindex="0"></pre></div></div>`;
    const liste = conteneur.querySelector(".code-liste"), pre = conteneur.querySelector("pre"), nom = conteneur.querySelector(".nom");
    const gh = conteneur.querySelector(".gh"), lg = conteneur.querySelector(".lg");
    let grp = null, courant = null;
    fichiers.forEach((f, k) => {
      if (f.groupe !== grp) { grp = f.groupe; const d = document.createElement("div"); d.className = "grp"; d.textContent = grp; liste.append(d); }
      const b = document.createElement("button"); b.type = "button"; b.textContent = f.chemin.split("/").pop();
      b.title = f.chemin; b.addEventListener("click", () => ouvrir(k)); f.btn = b; liste.append(b);
    });
    function ouvrir(k) {
      const f = fichiers[k]; courant = f;
      fichiers.forEach(x => x.btn.setAttribute("aria-current", String(x === f)));
      nom.textContent = f.chemin;
      lg.textContent = ({py: "Python", cpp: "C++", js: "JavaScript", vb: "VBA", txt: "texte"})[f.lang] || f.lang;
      const lien = f.lien || (depot ? depot + "/" + f.chemin.split("/").map(encodeURIComponent).join("/") : "");
      gh.hidden = !lien; gh.href = lien;
      pre.innerHTML = colorer(f.texte, f.lang);
      pre.scrollTop = 0;
    }
    conteneur.querySelector(".copier").addEventListener("click", ev => {
      if (!courant) return;
      const b = ev.currentTarget;
      const ok = () => { b.textContent = "Copié ✓"; setTimeout(() => (b.textContent = "Copier"), 1400); };
      if (navigator.clipboard) navigator.clipboard.writeText(courant.texte).then(ok, () => {});
    });
    if (fichiers.length) ouvrir(0);
  }

  // ---------- formulaires ----------
  const lireNombre = (el, defaut) => { const v = parseFloat(String(el.value).replace(",", ".")); return isFinite(v) ? v : defaut; };
  const listeNombres = s => String(s).split(/[\s;,]+/).map(x => x.trim()).filter(Boolean).map(Number);
  function fmt(v, d = 6) {
    if (v == null || !isFinite(v)) return "—";
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-4 || a >= 1e7)) return v.toExponential(Math.max(2, d - 2));
    return v.toFixed(d);
  }
  const pct = (v, d = 3) => isFinite(v) ? (v * 100).toFixed(d) + " %" : "—";
  const bp = (v, d = 2) => isFinite(v) ? (v * 1e4).toFixed(d) + " bp" : "—";
  // différer un calcul lourd après la saisie (évite de tout relancer à chaque frappe)
  function differe(f, ms = 250) { let t = 0; return (...a) => { clearTimeout(t); t = setTimeout(() => f(...a), ms); }; }
  // travail lourd hors du fil principal : un Worker construit depuis le texte d'un <script type="text/plain">
  function travailleur(...ids) {
    const src = ids.map(id => document.getElementById(id).textContent).join("\n;\n");
    const url = URL.createObjectURL(new Blob([src], {type: "text/javascript"}));
    return new Worker(url);
  }
  document.addEventListener("DOMContentLoaded", onglets);
  window.UI = {lecteurCode, lireNombre, listeNombres, fmt, pct, bp, differe, travailleur};
})();
