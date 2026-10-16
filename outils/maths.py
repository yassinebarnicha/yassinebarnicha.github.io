"""LaTeX formulas in the pages, rendered once at build time with KaTeX: no math script runs in the browser,
the page only loads kit/katex/katex.min.css and the few woff2 fonts the formulas use.

Markup in a page source, outside <script>:
    <m>\\sigma\\sqrt{T}</m>      inline formula
    <md>V = S\\,f(K/S)</md>      displayed formula (put it between paragraphs, not inside a <p>)
Inside a <script>, a double-quoted literal tex("...") is rendered as well; the results go into window.__TEX,
which the run-time function tex() reads (kit/ui.js). Use it for formulas written by JavaScript.

Rendered formulas are cached in outils/maths_cache.json, so a rebuild with unchanged formulas needs only Python.
A new formula needs Node.js (it runs outils/katex_rendu.js); a TeX error stops the build with the formula.
"""
import html as _html
import json
import re
import shutil
import subprocess
from pathlib import Path

ICI = Path(__file__).resolve().parent
CACHE = ICI / "maths_cache.json"
SCRIPT = re.compile(r"(<script\b[^>]*>)(.*?)(</script>)", re.S | re.I)
BALISE = re.compile(r"<(md|m)>(.*?)</\1>", re.S)
APPEL = re.compile(r'\btex\("((?:[^"\\\n]|\\.)*)"\)')


def _cle(tex, display):
    return ("D|" if display else "I|") + tex


def _rendre(demandes):
    """{key: (tex, display)} -> {key: html}, from the cache, Node only for what is missing."""
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    manquants = [k for k in demandes if k not in cache]
    if manquants:
        node = shutil.which("node")
        if not node:
            raise SystemExit("Node.js is needed to render these new formulas:\n  " +
                             "\n  ".join(demandes[k][0] for k in manquants[:10]))
        entree = json.dumps([[demandes[k][0], demandes[k][1]] for k in manquants])
        p = subprocess.run([node, str(ICI / "katex_rendu.js")], input=entree, capture_output=True,
                           text=True, encoding="utf-8")
        if p.returncode:
            raise SystemExit("KaTeX: " + p.stderr.strip().splitlines()[0] if p.stderr.strip() else "KaTeX failed")
        for k, h in zip(manquants, json.loads(p.stdout)):
            cache[k] = h
        CACHE.write_text(json.dumps(dict(sorted(cache.items())), ensure_ascii=False, indent=0) + "\n",
                         encoding="utf-8", newline="\n")
    return {k: cache[k] for k in demandes}


def _enrober(tex, display, rendu):
    label = _html.escape(tex, quote=True)
    if display:
        return f'<div class="tex-bloc" role="math" aria-label="{label}">{rendu}</div>'
    return f'<span class="tex" role="math" aria-label="{label}">{rendu}</span>'


def rendre_page(page, feuille):
    """Renders the formulas of an HTML page. `feuille` is the href of katex.min.css as seen from the page.
    Returns the page; adds the stylesheet link and window.__TEX only when the page has formulas."""
    morceaux = SCRIPT.split(page)
    # SCRIPT.split gives: text, <script>, body, </script>, text, ...
    demandes, appels = {}, []
    for i, m in enumerate(morceaux):
        if i % 4 == 0:
            for b in BALISE.finditer(m):
                tex = _html.unescape(b.group(2)).strip()
                demandes[_cle(tex, b.group(1) == "md")] = (tex, b.group(1) == "md")
        elif i % 4 == 2:
            for a in APPEL.finditer(m):
                tex = json.loads('"' + a.group(1) + '"')
                demandes[_cle(tex, False)] = (tex, False)
                appels.append(tex)
    if not demandes:
        return page
    rendus = _rendre(demandes)
    for i in range(0, len(morceaux), 4):
        morceaux[i] = BALISE.sub(lambda b: _enrober(_html.unescape(b.group(2)).strip(), b.group(1) == "md",
                                                    rendus[_cle(_html.unescape(b.group(2)).strip(), b.group(1) == "md")]),
                                 morceaux[i])
    page = "".join(morceaux)
    tete = f'<link rel="stylesheet" href="{feuille}">\n'
    if appels:
        table = {t: _enrober(t, False, rendus[_cle(t, False)]) for t in appels}
        tete += "<script>window.__TEX=" + json.dumps(table, ensure_ascii=False).replace("</", "<\\/") + ";</script>\n"
    if "</head>" not in page:
        raise SystemExit("no </head> in the page")
    return page.replace("</head>", tete + "</head>", 1)
