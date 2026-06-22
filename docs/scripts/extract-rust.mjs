#!/usr/bin/env node
// Extract rustdoc JSON for selected crates and reshape into a compact API model.
// Run from the docs/ directory; invokes `cargo +nightly rustdoc` at the workspace root.

import { execFileSync } from "node:child_process";
import { readFileSync, writeFileSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = dirname(fileURLToPath(import.meta.url));
const REPO = resolve(__dirname, "..", "..");
const CRATES = ["ppvm-traits", "ppvm-pauli-word", "ppvm-pauli-sum", "ppvm-tableau", "ppvm-sym"];

const CRATE_BLURB = {
  "ppvm-traits": "Trait system, the `Config` bundle, the `Pauli` alphabet, and map impls.",
  "ppvm-pauli-word": "Packed Pauli strings: `PauliWord`, phased, lossy, and pattern variants.",
  "ppvm-pauli-sum": "The `PauliSum` engine, truncation strategies, and concrete config bundles.",
  "ppvm-tableau": "Generalized stabilizer tableau simulator (Clifford + non-Clifford).",
  "ppvm-sym": "Symbolic, parametric Pauli propagation.",
};

const KIND_ORDER = ["module", "struct", "enum", "trait", "type_alias", "function", "macro"];

function runRustdoc(crate) {
  console.error(`[rust] building rustdoc JSON for ${crate}`);
  execFileSync(
    "cargo",
    ["+nightly", "rustdoc", "-p", crate, "--lib", "--", "-Z", "unstable-options", "--output-format", "json"],
    { cwd: REPO, stdio: ["ignore", "inherit", "inherit"] },
  );
  const path = resolve(REPO, "target", "doc", crate.replace(/-/g, "_") + ".json");
  return JSON.parse(readFileSync(path, "utf8"));
}

function itemKind(item) {
  if (!item.inner || typeof item.inner !== "object") return null;
  return Object.keys(item.inner)[0];
}

function pathOf(doc, id) {
  const p = doc.paths?.[id];
  if (!p) return null;
  return p.path.join("::");
}

// Render a minimal Rust signature line for a function/method.
function fnSignature(name, fn) {
  if (!fn?.sig) return name;
  const inputs = (fn.sig.inputs || [])
    .map(([n, _t]) => n)
    .join(", ");
  const ret = fn.sig.output ? " -> …" : "";
  return `fn ${name}(${inputs})${ret}`;
}

function shortType(generics) {
  if (!generics?.params?.length) return "";
  const names = generics.params.map((p) => p.name);
  return "<" + names.join(", ") + ">";
}

function summarize(docs) {
  if (!docs) return "";
  const firstPara = docs.split(/\n{2,}/)[0].trim();
  return firstPara.length > 240 ? firstPara.slice(0, 237) + "…" : firstPara;
}

// Convert a rustdoc `links` map (text → item_id) into a stable map
// (text → resolved path), so the front-end can rewrite `[`Foo`]` and
// `[`Foo`](crate::path::Foo)` into anchor links without knowing the
// rustdoc id space.
function resolveLinks(doc, links) {
  if (!links) return {};
  const out = {};
  for (const [text, id] of Object.entries(links)) {
    const target = doc.paths?.[id];
    if (!target || !target.path) continue;
    out[text] = target.path.join("::");
  }
  return out;
}

function processCrate(crateName, doc) {
  const items = [];
  for (const item of Object.values(doc.index)) {
    if (item.visibility !== "public") continue;
    if (!item.docs) continue; // only items with docstrings
    const kind = itemKind(item);
    if (!kind || !KIND_ORDER.includes(kind)) continue;
    const path = pathOf(doc, item.id);
    if (!path) continue;
    if (path.split("::")[0] !== crateName.replace(/-/g, "_")) continue;

    let signature = item.name;
    if (kind === "function") signature = fnSignature(item.name, item.inner.function);
    else if (kind === "struct") signature = "struct " + item.name + shortType(item.inner.struct.generics);
    else if (kind === "enum") signature = "enum " + item.name + shortType(item.inner.enum.generics);
    else if (kind === "trait") signature = "trait " + item.name + shortType(item.inner.trait.generics);
    else if (kind === "type_alias") signature = "type " + item.name;
    else if (kind === "macro") signature = item.name + "!";
    else if (kind === "module") signature = "mod " + item.name;

    items.push({
      name: item.name,
      kind,
      path,
      signature,
      summary: summarize(item.docs),
      docs: item.docs,
      span: item.span ? `${item.span.filename}:${item.span.begin[0]}` : null,
      links: resolveLinks(doc, item.links),
    });
  }
  items.sort((a, b) => (KIND_ORDER.indexOf(a.kind) - KIND_ORDER.indexOf(b.kind)) || a.path.localeCompare(b.path));
  return { crate: crateName, blurb: CRATE_BLURB[crateName] || "", items };
}

const result = { generatedAt: new Date().toISOString(), crates: [] };
for (const c of CRATES) {
  const doc = runRustdoc(c);
  result.crates.push(processCrate(c, doc));
}

const out = resolve(__dirname, "..", "src", "data", "rust-api.json");
mkdirSync(dirname(out), { recursive: true });
writeFileSync(out, JSON.stringify(result, null, 2));
console.error(`[rust] wrote ${out} (${result.crates.reduce((n, c) => n + c.items.length, 0)} items)`);
