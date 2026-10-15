"""Page « Swaptions bermudéennes sous HJM » (lue par construire.py)."""

# repère du gabarit -> moteur JavaScript de fabrique/
MOTEURS = {"HJM": "hjm.js"}

# code montré dans l'onglet « Code source » : (groupe, [motifs relatifs au dossier du projet])
FICHIERS = [
    ("Paquet Python hjm1f", ["code/hjm1f/__init__.py", "code/hjm1f/curve.py", "code/hjm1f/instruments.py",
                             "code/hjm1f/market_quotes.py", "code/hjm1f/hjm_1f.py", "code/hjm1f/solvers.py",
                             "code/hjm1f/calibration_1f.py", "code/hjm1f/pde_1f.py", "code/hjm1f/mc_1f.py"]),
    ("Exemples", ["code/run_example.py", "code/Main.ipynb"]),
    ("Portage JavaScript (cette page)", ["fabrique/hjm.js"]),
]
