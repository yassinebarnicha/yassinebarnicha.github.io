// Build-time LaTeX rendering (called by outils/maths.py, never loaded by a page).
// stdin : JSON [[tex, display], ...]   stdout : JSON [html, ...]
// KaTeX 0.19.0 (MIT licence, see kit/katex/LICENSE) is vendored next to this file.
"use strict";
const katex = require("./katex.min.js");
let entree = "";
process.stdin.setEncoding("utf8");
process.stdin.on("data", d => (entree += d));
process.stdin.on("end", () => {
  const sortie = JSON.parse(entree).map(([tex, display]) =>
    katex.renderToString(tex, {displayMode: !!display, throwOnError: true, output: "html", strict: "ignore"}));
  process.stdout.write(JSON.stringify(sortie));
});
