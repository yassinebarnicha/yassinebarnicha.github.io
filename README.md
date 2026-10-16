# yassinebarnicha.github.io

Home page of my projects and their online applications: **https://yassinebarnicha.github.io/**

| Page | Address | Report |
|---|---|---|
| Home: all projects | https://yassinebarnicha.github.io/ | |
| Bermudan swaptions under HJM | https://yassinebarnicha.github.io/swaptions-hjm/ | [PDF](swaptions-hjm/rapport/report.pdf) |
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
```

After changing a template, an engine or the kit: `python construire.py` (all pages) or `python construire.py swaptions-hjm`
(one page). Python standard library only; Node.js is needed only when a formula is new (rendered formulas are cached in
`outils/maths_cache.json`). Reports: `pdflatex report.tex` twice in `<project>/rapport/`.

Published by GitHub Pages from the `main` branch, root folder.

Mohammed-Yassine Barnicha — [MIT](LICENSE) licence. KaTeX (`kit/katex/`, `outils/katex.min.js`) is under its own MIT licence.
