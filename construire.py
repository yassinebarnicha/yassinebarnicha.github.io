"""Construit les pages de projets du site : <projet>/index.html à partir de <projet>/fabrique/.

Chaque dossier de projet contient :
  fabrique/page.html   le gabarit de la page
  fabrique/config.py   MOTEURS (repère -> fichier JavaScript de fabrique/) et FICHIERS (code montré dans l'onglet « Code »)
  fabrique/*.js        le moteur de calcul
  code/                le code d'origine du projet

Repères remplacés dans le gabarit : /*__STYLE__*/, /*__GRAPHE__*/, /*__UI__*/ (fichiers de kit/), un repère par moteur,
/*__CODE__*/ (liste JSON des fichiers source), __DEPOT__ (adresse du dossier du projet sur GitHub).

Usage : python construire.py            (toutes les pages)
        python construire.py swaptions-hjm   (une seule)
Aucune dépendance hors bibliothèque standard.
"""
import importlib.util
import json
import re
import sys
from pathlib import Path

RACINE = Path(__file__).resolve().parent
KIT = RACINE / "kit"
GITHUB = "https://github.com/yassinebarnicha/yassinebarnicha.github.io"
LANG = {".py": "py", ".ipynb": "py", ".cpp": "cpp", ".h": "cpp", ".hpp": "cpp", ".js": "js", ".bas": "vb"}


def lire_source(chemin: Path) -> str:
    """Texte d'un fichier ; pour un notebook, ses cellules (code et texte) mises bout à bout, sans les sorties."""
    if chemin.suffix == ".ipynb":
        nb = json.loads(chemin.read_text(encoding="utf-8"))
        morceaux = []
        for i, cel in enumerate(nb.get("cells", [])):
            src = "".join(cel.get("source", []))
            if not src.strip():
                continue
            if cel.get("cell_type") == "code":
                morceaux.append(f"# ---- cellule {i + 1} ----\n{src}")
            else:
                morceaux.append(f"# ---- cellule {i + 1} (texte) ----\n" + "\n".join("# " + l for l in src.splitlines()))
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
                raise SystemExit(f"{dossier.name} : fichier introuvable : {motif}")
            for f in trouves:
                liste.append({"chemin": f.relative_to(dossier).as_posix(), "groupe": groupe,
                              "lang": LANG.get(f.suffix.lower(), "txt"), "texte": lire_source(f)})
    # « </ » échappé pour ne jamais fermer la balise <script> qui contient le JSON
    return json.dumps(liste, ensure_ascii=False).replace("</", "<\\/")


def construire(dossier: Path):
    cfg = charger_config(dossier)
    html = (dossier / "fabrique" / "page.html").read_text(encoding="utf-8")
    # bandeau « ceci est le gabarit » : seulement pour qui ouvrirait page.html par erreur
    html = re.sub(r"<!--GABARIT-->.*?<!--/GABARIT-->\s*", "", html, flags=re.S)
    contenus = {"STYLE": KIT / "style.css", "GRAPHE": KIT / "graphe.js", "UI": KIT / "ui.js"}
    contenus.update({cle: dossier / "fabrique" / nom for cle, nom in cfg.MOTEURS.items()})
    textes = {}
    for cle, f in contenus.items():
        t = f.read_text(encoding="utf-8")
        if re.search(r"</script", t, re.I):
            raise SystemExit(f"{f} contient « </script » : impossible à intégrer tel quel")
        textes[cle] = t
    reperes = re.findall(r"/\*__([A-Z_]+)__\*/", html)
    inconnus = [r for r in reperes if r not in textes and r != "CODE"]
    absents = [c for c in textes if c not in reperes]
    if inconnus or absents:
        raise SystemExit(f"{dossier.name} : repères sans contenu {inconnus}, contenus sans repère {absents}")
    depot = f"{GITHUB}/blob/main/{dossier.name}"
    html = html.replace("__DEPOT__", depot)
    # un seul passage : le texte inséré n'est jamais relu (il peut lui-même citer ces repères)
    html = re.sub(r"/\*__([A-Z_]+)__\*/",
                  lambda m: code_json(dossier, cfg.FICHIERS) if m.group(1) == "CODE" else textes[m.group(1)], html)
    sortie = dossier / "index.html"
    sortie.write_text(html, encoding="utf-8", newline="\n")
    print(f"{sortie.relative_to(RACINE).as_posix()} ({len(html) // 1024} Ko)")


def main():
    projets = sorted(p.parent.parent for p in RACINE.glob("*/fabrique/config.py"))
    if len(sys.argv) > 1:
        projets = [p for p in projets if p.name in sys.argv[1:]]
        if not projets:
            raise SystemExit(f"projet inconnu : {sys.argv[1:]}")
    for p in projets:
        construire(p)


if __name__ == "__main__":
    main()
