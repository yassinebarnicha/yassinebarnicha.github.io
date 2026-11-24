// Encrypts the interview-prep PDFs for the password-protected hub (interview-prep/).
//
// The repository is public: only encrypted files are published. Each PDF is encrypted with AES-256-GCM,
// with a key derived from the password by PBKDF2-SHA256; the page derives the same key in the browser (WebCrypto)
// and decrypts the PDF locally. Without the password the files are unreadable.
//
// Usage:  node outils/chiffrer_hub.mjs <folder>
//   <folder> holds catalogue.json (list of {id, part, title, summary, file, pages, images}) with the PDFs and the
//   page images (WebP) it names
//   (produced by construire_matieres.py, kept outside this repository).
//   The password is asked twice (hidden). For a non-interactive run: environment variable HUB_PASSWORD.
// Output: interview-prep/data/  key.json (salt, iterations), catalogue.bin, one <random>.bin per PDF.
// Previous .bin files are removed. Node.js standard library only.

import { createCipheriv, pbkdf2Sync, randomBytes } from "node:crypto";
import { mkdirSync, readFileSync, readdirSync, rmSync, writeFileSync } from "node:fs";
import { dirname, join, resolve } from "node:path";
import { createInterface } from "node:readline";
import { fileURLToPath } from "node:url";

const ITERATIONS = 310000;
const RACINE = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const SORTIE = join(RACINE, "interview-prep", "data");

function demander(question) {
  return new Promise((ok) => {
    const rl = createInterface({ input: process.stdin, output: process.stdout, terminal: true });
    rl._writeToOutput = (s) => { if (s.startsWith(question)) process.stdout.write(question); };
    rl.question(question, (rep) => { rl.close(); process.stdout.write("\n"); ok(rep); });
  });
}

async function motDePasse() {
  if (process.env.HUB_PASSWORD) return process.env.HUB_PASSWORD;
  const a = await demander("Password: ");
  const b = await demander("Password again: ");
  if (a !== b) throw new Error("the two passwords differ");
  return a;
}

function chiffrer(cle, donnees) {
  const iv = randomBytes(12);
  const c = createCipheriv("aes-256-gcm", cle, iv);
  const corps = Buffer.concat([c.update(donnees), c.final()]);
  return Buffer.concat([iv, corps, c.getAuthTag()]);   // iv | ciphertext | tag (layout expected by WebCrypto)
}

const dossier = process.argv[2];
if (!dossier) {
  console.error("usage: node outils/chiffrer_hub.mjs <folder with catalogue.json and the PDFs>");
  process.exit(1);
}
const catalogue = JSON.parse(readFileSync(join(dossier, "catalogue.json"), "utf8"));
const mdp = await motDePasse();
if (!mdp) throw new Error("empty password");
if (mdp.length < 8) console.warn("warning: a short password can be guessed in seconds; the encryption is only as strong as the password");

const sel = randomBytes(16);
const cle = pbkdf2Sync(mdp.normalize("NFC"), sel, ITERATIONS, 32, "sha256");

mkdirSync(SORTIE, { recursive: true });
for (const f of readdirSync(SORTIE)) if (f.endsWith(".bin")) rmSync(join(SORTIE, f));

const liste = [];
for (const m of catalogue) {
  const pdf = readFileSync(join(dossier, m.file));
  const nom = randomBytes(8).toString("hex") + ".bin";
  writeFileSync(join(SORTIE, nom), chiffrer(cle, pdf));
  // page images (WebP) shown inside the hub page; the PDF itself stays available for download
  const images = (m.images || []).map((im) => {
    const b = randomBytes(8).toString("hex") + ".bin";
    writeFileSync(join(SORTIE, b), chiffrer(cle, readFileSync(join(dossier, im.file))));
    return { blob: b, w: im.w, h: im.h };
  });
  liste.push({ id: m.id, part: m.part, title: m.title, summary: m.summary, pages: m.pages, bytes: pdf.length, blob: nom, images });
  console.log(`${m.id.padEnd(22)} ${String(pdf.length).padStart(9)} bytes, ${images.length} page images -> ${nom}`);
}
writeFileSync(join(SORTIE, "catalogue.bin"), chiffrer(cle, Buffer.from(JSON.stringify(liste), "utf8")));
writeFileSync(join(SORTIE, "key.json"),
  JSON.stringify({ version: 1, kdf: "PBKDF2-SHA256", iterations: ITERATIONS, salt: sel.toString("base64"), cipher: "AES-256-GCM" }, null, 1) + "\n");
console.log(`${liste.length} documents encrypted into interview-prep/data/`);
