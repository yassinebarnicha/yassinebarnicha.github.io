# yassinebarnicha.github.io

Home page of my projects and their online applications: **https://yassinebarnicha.github.io/**

| Page | Address | Report |
|---|---|---|
| Home: all projects | https://yassinebarnicha.github.io/ | |
| Bermudan swaptions under HJM (1 factor, 2 factors, 2 factors + stochastic volatility) | https://yassinebarnicha.github.io/swaptions-hjm/ | [PDF](swaptions-hjm/rapport/report.pdf) |
| Lookback options | https://yassinebarnicha.github.io/options-lookback/ | [PDF](options-lookback/rapport/report.pdf) |
| Bitcoin option calibration (separate repository) | https://yassinebarnicha.github.io/calibration-btc-svjd/ | in its repository |

Each page runs entirely in the browser, with no server and nothing to install, and shows the original code of the project
("Source code" tab). Formulas are written in LaTeX and rendered once, at build time, with KaTeX: the pages load no math script.

## Layout

```
index.html               home page
img/                     screenshots used by the home page
kit/                     style, charts, code viewer and KaTeX stylesheet/fonts shared by every page
outils/                  build-time LaTeX rendering (maths.py, katex_rendu.js, KaTeX, cache of rendered formulas)
construire.py            builds <project>/index.html for each project
<project>/index.html     the page (generated: do not edit by hand)
<project>/fabrique/      page.html (template), config.py, JavaScript computation engine
<project>/code/          the original code of the project
<project>/rapport/       the report: report.tex, report.pdf and the data of its figures
interview-prep/          password-protected hub: index.html, verrou.js (lock), data/ (encrypted PDFs)
outils/chiffrer_hub.mjs  encrypts the interview-prep PDFs (Node.js)
```

The Bermudan swaptions page has a model switch (1 factor / 2 factors / 2 factors + stochastic volatility) with its own set of tabs for each model, a shared *Source code* tab,
and three JavaScript engines (`hjm.js`, `hjm2.js`, `hjm3.js`) plus the market data of its real-data example (`donnees.js`).

After changing a template, an engine or the kit: `python construire.py` (all pages) or `python construire.py swaptions-hjm`
(one page). Python standard library only; Node.js is needed only when a formula is new (rendered formulas are cached in
`outils/maths_cache.json`). Reports: `pdflatex report.tex` twice in `<project>/rapport/`.

## Interview prep (password-protected)

`interview-prep/` is a private hub of interview notes: one PDF per topic. The lock button at the top right of the home
page asks for the password; once unlocked, an *Interview prep* link appears and the hub lists the topics, each opening
it inside the page (pages shown as images, so phones without a PDF viewer work too; *PDF ↓* downloads the PDF; the
phone's back button closes it). Pressing the button again, or *Lock* on the hub, locks it.

The repository is public, so the PDFs and their page images are published **encrypted** (`interview-prep/data/*.bin`, AES-256-GCM, key derived
from the password with PBKDF2-SHA256, 310 000 iterations) and decrypted in the browser (WebCrypto, no library). The
LaTeX sources and the plain PDFs are not in this repository. The key stays on the device (localStorage) until locked.

To change the password or update the documents:

```
node outils/chiffrer_hub.mjs <folder with catalogue.json and the PDFs>
```

The script asks the password twice, replaces `interview-prep/data/` and needs Node.js only (standard library).
Changing the password locks every device that had unlocked with the previous one.

Published by GitHub Pages from the `main` branch, root folder.

Mohammed-Yassine Barnicha — [MIT](LICENSE) licence. KaTeX (`kit/katex/`, `outils/katex.min.js`) is under its own MIT licence.
