"""Page "Bermudan swaptions under HJM" (read by construire.py)."""

# template marker -> JavaScript file of fabrique/ (HJM, HJM2: computation engines run in a worker; DONNEES: market data)
MOTEURS = {"HJM": "hjm.js", "HJMDEUX": "hjm2.js", "DONNEES": "donnees.js"}

# code shown in the "Source code" tab: (group, [patterns relative to the project folder])
FICHIERS = [
    ("Python package hjm1f (1 factor, corrected)", ["code/hjm1f/__init__.py", "code/hjm1f/curve.py", "code/hjm1f/instruments.py",
                              "code/hjm1f/market_quotes.py", "code/hjm1f/hjm_1f.py", "code/hjm1f/solvers.py",
                              "code/hjm1f/calibration_1f.py", "code/hjm1f/pde_1f.py", "code/hjm1f/mc_1f.py"]),
    ("Python package hjm2f (2 factors)", ["code/hjm2f/__init__.py", "code/hjm2f/hjm_2f.py", "code/hjm2f/calibration_2f.py",
                              "code/hjm2f/calibration_global.py", "code/hjm2f/optimize.py", "code/hjm2f/quad_2f.py",
                              "code/hjm2f/mc_2f.py", "code/hjm2f/caps.py", "code/hjm2f/market_data.py"]),
    ("Examples", ["code/run_example.py", "code/run_example_2f.py", "code/run_market_2f.py", "code/check_1f_consistency.py", "code/validate_2f.py",
                  "code/Main.ipynb"]),
    ("Original code, before the review", ["code/original/hjm_1f.py", "code/original/mc_1f.py"]),
    ("JavaScript ports (this page)", ["fabrique/hjm.js", "fabrique/hjm2.js"]),
    ("Report (LaTeX source)", ["rapport/report.tex"]),
]
