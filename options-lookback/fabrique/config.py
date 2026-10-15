"""Page « Options lookback » (lue par construire.py)."""

# repère du gabarit -> moteur JavaScript de fabrique/
MOTEURS = {"LB": "lookback.js"}

# code montré dans l'onglet « Code source » : (groupe, [motifs relatifs au dossier du projet])
FICHIERS = [
    ("Moteur Barnicha : pont brownien exact", ["code/Barnicha/include/lookback/*.h", "code/Barnicha/src/*.cpp",
                                               "code/Barnicha/app/main.cpp"]),
    ("Moteur Chivet : schéma en temps + BGK", ["code/Chivet/include/*.hpp", "code/Chivet/src/*.cpp", "code/Chivet/app/main.cpp"]),
    ("DLL unifiée et Excel", ["code/unified/unified_dll.cpp", "code/CMakeLists.txt", "code/excel/modLookback_ttm.bas"]),
    ("Portage JavaScript (cette page)", ["fabrique/lookback.js"]),
]
