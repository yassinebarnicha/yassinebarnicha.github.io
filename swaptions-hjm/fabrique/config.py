"""Page "Bermudan swaptions under HJM" (read by construire.py)."""

# template marker -> JavaScript engine of fabrique/
MOTEURS = {"HJM": "hjm.js"}

# code shown in the "Source code" tab: (group, [patterns relative to the project folder])
FICHIERS = [
    ("Python package hjm1f", ["code/hjm1f/__init__.py", "code/hjm1f/curve.py", "code/hjm1f/instruments.py",
                              "code/hjm1f/market_quotes.py", "code/hjm1f/hjm_1f.py", "code/hjm1f/solvers.py",
                              "code/hjm1f/calibration_1f.py", "code/hjm1f/pde_1f.py", "code/hjm1f/mc_1f.py"]),
    ("Examples", ["code/run_example.py", "code/Main.ipynb"]),
    ("JavaScript port (this page)", ["fabrique/hjm.js"]),
    ("Report (LaTeX source)", ["rapport/report.tex"]),
]
