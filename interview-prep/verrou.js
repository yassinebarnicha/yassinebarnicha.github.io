/* Lock of the interview-prep hub, shared by the home page and interview-prep/index.html.
   The PDFs are published encrypted (AES-256-GCM, key = PBKDF2-SHA256 of the password, see outils/chiffrer_hub.mjs).
   Unlocking derives the key in the browser and keeps it on this device (localStorage) until "Lock" is pressed.
   No external library: WebCrypto only (needs https or localhost). */
(function () {
  "use strict";
  var DATA = new URL("data/", document.currentScript.src);
  var CLE_STOCKAGE = "yb.interview-prep.key";

  function b64(buf) { var s = "", o = new Uint8Array(buf); for (var i = 0; i < o.length; i++) s += String.fromCharCode(o[i]); return btoa(s); }
  function deB64(s) { var b = atob(s), o = new Uint8Array(b.length); for (var i = 0; i < b.length; i++) o[i] = b.charCodeAt(i); return o; }
  function lire() { try { return JSON.parse(localStorage.getItem(CLE_STOCKAGE) || "null"); } catch (e) { return null; } }
  function ecrire(v) { try { v ? localStorage.setItem(CLE_STOCKAGE, JSON.stringify(v)) : localStorage.removeItem(CLE_STOCKAGE); } catch (e) { /* private mode: memory only */ } memoire = v; }
  var memoire = lire();

  function disponible() { return !!(window.crypto && window.crypto.subtle); }

  function chercher(nom, type) {
    return fetch(new URL(nom, DATA), { cache: "no-cache" }).then(function (r) {
      if (!r.ok) throw new Error("missing file " + nom + " (" + r.status + ")");
      return type === "json" ? r.json() : r.arrayBuffer();
    });
  }

  function dechiffrer(cle, buf) {
    var o = new Uint8Array(buf);
    return crypto.subtle.decrypt({ name: "AES-GCM", iv: o.slice(0, 12) }, cle, o.slice(12));
  }

  function importer(brut) { return crypto.subtle.importKey("raw", brut, { name: "AES-GCM" }, false, ["decrypt"]); }

  /* password -> key, checked by decrypting the catalogue; resolves to the catalogue (list of documents) */
  function deverrouiller(motDePasse) {
    if (!disponible()) return Promise.reject(new Error("this browser has no WebCrypto (https needed)"));
    var meta;
    return chercher("key.json", "json").then(function (m) {
      meta = m;
      return crypto.subtle.importKey("raw", new TextEncoder().encode(motDePasse.normalize("NFC")), "PBKDF2", false, ["deriveBits"]);
    }).then(function (base) {
      return crypto.subtle.deriveBits({ name: "PBKDF2", hash: "SHA-256", salt: deB64(meta.salt), iterations: meta.iterations }, base, 256);
    }).then(function (bits) {
      return importer(bits).then(function (cle) {
        return chercher("catalogue.bin").then(function (c) { return dechiffrer(cle, c); }).then(function (clair) {
          ecrire({ salt: meta.salt, key: b64(bits) });
          return JSON.parse(new TextDecoder().decode(clair));
        }, function () { throw new Error("wrong password"); });
      });
    });
  }

  /* stored key -> catalogue; null when locked or when the files were re-encrypted with another password */
  function catalogue() {
    if (!memoire || !disponible()) return Promise.resolve(null);
    var stocke = memoire;
    return chercher("key.json", "json").then(function (m) {
      if (m.salt !== stocke.salt) { ecrire(null); return null; }
      return importer(deB64(stocke.key)).then(function (cle) {
        return chercher("catalogue.bin").then(function (c) { return dechiffrer(cle, c); });
      }).then(function (clair) { return JSON.parse(new TextDecoder().decode(clair)); });
    }).catch(function (e) { if (/missing file/.test(e.message)) throw e; ecrire(null); return null; });
  }

  /* decrypts one encrypted file; resolves to a blob: URL of the given type */
  var cleCache = null;
  function dechiffrerURL(blob, type) {
    if (!memoire) return Promise.reject(new Error("locked"));
    if (!cleCache || cleCache.k !== memoire.key) cleCache = { k: memoire.key, p: importer(deB64(memoire.key)) };
    return cleCache.p.then(function (cle) {
      return chercher(blob).then(function (c) { return dechiffrer(cle, c); });
    }).then(function (clair) { return URL.createObjectURL(new Blob([clair], { type: type })); });
  }
  function ouvrir(doc) { return dechiffrerURL(doc.blob, "application/pdf"); }           // the PDF (download)
  function page(image) { return dechiffrerURL(image.blob, "image/webp"); }              // one page image (viewer)

  window.VerrouHub = {
    estOuvert: function () { return !!memoire; },
    deverrouiller: deverrouiller,
    verrouiller: function () { ecrire(null); },
    catalogue: catalogue,
    ouvrir: ouvrir,
    page: page,
    disponible: disponible
  };
})();
