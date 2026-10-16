"""Builds the project pages of the site: <project>/index.html from <project>/fabrique/.

Each project folder contains:
  fabrique/page.html   the page template
  fabrique/config.py   MOTEURS (marker -> JavaScript file of fabrique/) and FICHIERS (code shown in the "Source code" tab)
  fabrique/*.js        the computation engine
  code/                the original code of the project

Markers replaced in the template: /*__STYLE__*/, /*__GRAPHE__*/, /*__UI__*/ (files of kit/), one marker per engine,
/*__CODE__*/ (JSON list of the source files), __DEPOT__ (address of the project folder on GitHub).
LaTeX written as <m>…</m>, <md>…</md> or tex("…") is rendered at build time (outils/maths.py).

Usage: python construire.py                 (all pages)
       python construire.py swaptions-hjm   (one page)
Python standard library only; Node.js only when a formula is new (rendered formulas are cached).
"""
import importlib.util
import json
import re
import sys
from pathlib import Path

RACINE = Path(__file__).resolve().parent
KIT = RACINE / "kit"
GITHUB = "https://github.com/yassinebarnicha/yassinebarnicha.github.io"
LANG = {".py": "py", ".ipynb": "py", ".cpp": "cpp", ".h": "cpp", ".hpp": "cpp", ".js": "js", ".bas": "vb", ".tex": "tex"}
sys.path.insert(0, str(RACINE / "outils"))
import maths  # noqa: E402


def lire_source(chemin: Path) -> str:
    """Text of a file; for a notebook, its cells (code and text) one after the other, without the outputs."""
    if chemin.suffix == ".ipynb":
        nb = json.loads(chemin.read_text(encoding="utf-8"))
        morceaux = []
        for i, cel in enumerate(nb.get("cells", [])):
            src = "".join(cel.get("source", []))
            if not src.strip():
                continue
            if cel.get("cell_type") == "code":
                morceaux.append(f"# ---- cell {i + 1} ----\n{src}")
            else:
                morceaux.append(f"# ---- cell {i + 1} (text) ----\n" + "\n".join("# " + l for l in src.splitlines()))
        return "\n\n".join(morceaux) + "\n"
    return chemin.read_bytes().decode("utf-8", errors="replace").replace("\r\n", "\n")


def charger_config(dossier: Path):
    spec = importlib.util.spec_from_file_location(f"config_{dossier.name}", dossier / "fabrique" / "config.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def code_json(dossier: Path, fichiers) -> str:
    liste = []
    for groupe, motifs in fichiers:
        for motif in motifs:
            trouves = sorted(p for p in dossier.glob(motif) if p.is_file())
            if not trouves:
                raise SystemExit(f"{dossier.name}: file not found: {motif}")
            for f in trouves:
                liste.append({"chemin": f.relative_to(dossier).as_posix(), "groupe": groupe,
                              "lang": LANG.get(f.suffix.lower(), "txt"), "texte": lire_source(f)})
    # "</" escaped so that the JSON never closes the <script> element holding it
    return json.dumps(liste, ensure_ascii=False).replace("</", "<\\/")


def construire(dossier: Path):
    cfg = charger_config(dossier)
    html = (dossier / "fabrique" / "page.html").read_text(encoding="utf-8")
    # "this is the template" banner: only for someone opening page.html by mistake
    html = re.sub(r"<!--GABARIT-->.*?<!--/GABARIT-->\s*", "", html, flags=re.S)
    # formulas first, on the page's own text only (never on the inserted kit, engines or source code)
    html = maths.rendre_page(html, "../kit/katex/katex.min.css")
    contenus = {"STYLE": KIT / "style.css", "GRAPHE": KIT / "graphe.js", "UI": KIT / "ui.js"}
    contenus.update({cle: dossier / "fabrique" / nom for cle, nom in cfg.MOTEURS.items()})
    textes = {}
    for cle, f in contenus.items():
        t = f.read_text(encoding="utf-8")
        if re.search(r"</script", t, re.I):
            raise SystemExit(f"{f} contains \"</script\": it cannot be inlined as is")
        textes[cle] = t
    reperes = re.findall(r"/\*__([A-Z_]+)__\*/", html)
    inconnus = [r for r in reperes if r not in textes and r != "CODE"]
    absents = [c for c in textes if c not in reperes]
    if inconnus or absents:
        raise SystemExit(f"{dossier.name}: markers without content {inconnus}, content without marker {absents}")
    depot = f"{GITHUB}/blob/main/{dossier.name}"
    html = html.replace("__DEPOT__", depot)
    # a single pass: inserted text is never read again (it may itself quote these markers)
    html = re.sub(r"/\*__([A-Z_]+)__\*/",
                  lambda m: code_json(dossier, cfg.FICHIERS) if m.group(1) == "CODE" else textes[m.group(1)], html)
    sortie = dossier / "index.html"
    sortie.write_text(html, encoding="utf-8", newline="\n")
    print(f"{sortie.relative_to(RACINE).as_posix()} ({len(html) // 1024} KB)")


def main():
    projets = sorted(p.parent.parent for p in RACINE.glob("*/fabrique/config.py"))
    if len(sys.argv) > 1:
        projets = [p for p in projets if p.name in sys.argv[1:]]
        if not projets:
            raise SystemExit(f"unknown project: {sys.argv[1:]}")
    for p in projets:
        construire(p)


if __name__ == "__main__":
    main()
