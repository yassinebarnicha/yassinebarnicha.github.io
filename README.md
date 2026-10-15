# yassinebarnicha.github.io

Page d'accueil de mes projets et leurs applications en ligne : **https://yassinebarnicha.github.io/**

| Page | Adresse |
|---|---|
| Accueil : tous les projets | https://yassinebarnicha.github.io/ |
| Swaptions bermudéennes sous HJM | https://yassinebarnicha.github.io/swaptions-hjm/ |
| Options lookback | https://yassinebarnicha.github.io/options-lookback/ |
| Calibration d'options Bitcoin (dépôt à part) | https://yassinebarnicha.github.io/calibration-btc-svjd/ |

Chaque page tourne entièrement dans le navigateur, sans serveur ni installation, et montre le code d'origine du projet
(onglet « Code source »).

## Organisation

```
index.html               page d'accueil
img/                     captures de la page d'accueil
kit/                     style, graphiques et lecteur de code communs à toutes les pages
construire.py            construit <projet>/index.html pour chaque projet
<projet>/index.html      la page (générée : ne pas modifier à la main)
<projet>/fabrique/       page.html (gabarit), config.py, moteur de calcul en JavaScript
<projet>/code/           le code d'origine du projet
```

Après une modification d'un gabarit, d'un moteur ou du kit : `python construire.py` (toutes les pages) ou
`python construire.py swaptions-hjm` (une seule). Aucune dépendance hors bibliothèque standard.

Publiée par GitHub Pages depuis la branche `main`, dossier racine.

Mohammed-Yassine Barnicha — licence [MIT](LICENSE).
