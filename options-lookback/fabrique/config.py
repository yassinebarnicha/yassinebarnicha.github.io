"""Page "Lookback options" (read by construire.py)."""

# template marker -> JavaScript engine of fabrique/
MOTEURS = {"LB": "lookback.js"}

# code shown in the "Source code" tab: (group, [patterns relative to the project folder])
FICHIERS = [
    ("Barnicha engine: exact Brownian bridge", ["code/Barnicha/include/lookback/*.h", "code/Barnicha/src/*.cpp",
                                                "code/Barnicha/app/main.cpp"]),
    ("Chivet engine: time stepping + BGK", ["code/Chivet/include/*.hpp", "code/Chivet/src/*.cpp", "code/Chivet/app/main.cpp"]),
    ("Unified DLL and Excel", ["code/unified/unified_dll.cpp", "code/CMakeLists.txt", "code/excel/modLookback_ttm.bas"]),
    ("JavaScript port (this page)", ["fabrique/lookback.js"]),
    ("Report (LaTeX source)", ["rapport/report.tex"]),
]
