#!/usr/bin/env node
// Extract Python API + docstrings for the `ppvm` package via `griffe dump`.
// griffe is invoked through `uv run --with griffe` so no global install is needed.

import { execFileSync } from "node:child_process";
import { writeFileSync, mkdirSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = dirname(fileURLToPath(import.meta.url));
const REPO = resolve(__dirname, "..", "..");
const PY_PROJECT = resolve(REPO, "ppvm-python");

function runGriffe() {
  console.error("[python] running griffe dump ppvm (google-style)");
  const buf = execFileSync(
    "uv",
    [
      "run", "--project", PY_PROJECT, "--with", "griffe",
      "griffe", "dump", "ppvm",
      "-s", resolve(PY_PROJECT, "src"),
      "-d", "google",
      "-f",
    ],
    { stdio: ["ignore", "pipe", "inherit"], maxBuffer: 200 * 1024 * 1024 },
  );
  const all = JSON.parse(buf.toString("utf8"));
  return all.ppvm;
}

function exprName(e) {
  if (e == null) return "";
  if (typeof e === "string") return e;
  if (e.name) return e.name;
  if (e.canonical_path) return e.canonical_path;
  return "";
}

function fmtParam(p) {
  let s = p.name;
  const ann = exprName(p.annotation);
  if (ann) s += ": " + ann;
  if (p.default !== null && p.default !== undefined) s += " = " + p.default;
  return s;
}

function fnSignature(name, obj, isMethod) {
  let params = (obj.parameters || []).map(fmtParam);
  if (isMethod && params[0] && params[0].startsWith("self")) params = params.slice(1);
  const ret = exprName(obj.returns);
  return `${name}(${params.join(", ")})${ret ? " -> " + ret : ""}`;
}

function summarize(docs) {
  if (!docs) return "";
  const para = docs.split(/\n{2,}/)[0].trim().replace(/\s+/g, " ");
  return para.length > 240 ? para.slice(0, 237) + "…" : para;
}

// Normalize griffe's parsed sections into a stable shape for the front-end.
// Consecutive text sections get merged; multi-record return descriptions
// (a quirk of the Google parser) get coalesced.
function getParsedSections(obj) {
  const parsed = obj?.docstring?.parsed;
  if (!Array.isArray(parsed)) return null;
  const sections = [];
  for (const sec of parsed) {
    if (!sec || !sec.kind) continue;
    if (sec.kind === "text") {
      const value = typeof sec.value === "string" ? sec.value : "";
      const last = sections[sections.length - 1];
      if (last && last.kind === "text") {
        last.value += "\n\n" + value;
      } else {
        sections.push({ kind: "text", value });
      }
    } else if (sec.kind === "parameters") {
      const items = (sec.value || []).map((p) => ({
        name: p.name || "",
        annotation: exprName(p.annotation),
        default: p.default == null ? null : String(p.default),
        description: typeof p.description === "string" ? p.description : "",
      }));
      sections.push({ kind: "parameters", items });
    } else if (sec.kind === "returns" || sec.kind === "yields") {
      const merged = (sec.value || []).map((r) => ({
        name: r.name || "",
        annotation: exprName(r.annotation),
        description: typeof r.description === "string" ? r.description : "",
      }));
      const grouped = [];
      for (const entry of merged) {
        const last = grouped[grouped.length - 1];
        if (last && last.name === entry.name && last.annotation === entry.annotation) {
          last.description = last.description
            ? last.description + " " + entry.description
            : entry.description;
        } else {
          grouped.push({ ...entry });
        }
      }
      sections.push({ kind: sec.kind, items: grouped });
    } else if (sec.kind === "raises") {
      const items = (sec.value || []).map((r) => ({
        annotation: exprName(r.annotation),
        description: typeof r.description === "string" ? r.description : "",
      }));
      sections.push({ kind: "raises", items });
    } else if (sec.kind === "examples" || sec.kind === "admonition") {
      sections.push({
        kind: sec.kind,
        value: typeof sec.value === "string" ? sec.value : JSON.stringify(sec.value),
      });
    } else {
      sections.push({ kind: sec.kind, value: JSON.stringify(sec.value) });
    }
  }
  return sections;
}

const items = [];

function walk(obj, parentKind) {
  if (!obj || typeof obj !== "object") return;
  const kind = obj.kind;
  if (kind === "alias") return; // skip re-exports; we cover them in their original location
  if (obj.is_private) return;
  if (obj.name && obj.name.startsWith("_") && !obj.name.startsWith("__")) return;
  if (kind === "function" && obj.name && /^__\w+__$/.test(obj.name)) return; // dunder methods

  const docs = obj.docstring?.value || null;
  const path = obj.path;

  const displayKind =
    kind === "function" && parentKind === "class" ? "method" : kind;

  if (docs && path && (kind === "module" || kind === "class" || kind === "function" || kind === "attribute")) {
    let signature = obj.name;
    if (kind === "function") signature = (displayKind === "method" ? "" : "def ") + fnSignature(obj.name, obj, displayKind === "method");
    else if (kind === "class") signature = "class " + obj.name;
    else if (kind === "module") signature = "module " + path;
    else if (kind === "attribute") {
      const ann = exprName(obj.annotation);
      signature = obj.name + (ann ? ": " + ann : "");
    }
    items.push({
      name: obj.name,
      kind: displayKind,
      path,
      signature,
      summary: summarize(docs),
      docs,
      sections: getParsedSections(obj),
    });
  }

  if (obj.members) for (const child of Object.values(obj.members)) walk(child, kind);
}

const tree = runGriffe();
walk(tree, null);

const KIND_ORDER = ["module", "class", "function", "method", "attribute"];
items.sort((a, b) => {
  const ka = KIND_ORDER.indexOf(a.kind);
  const kb = KIND_ORDER.indexOf(b.kind);
  if (ka !== kb) return ka - kb;
  return a.path.localeCompare(b.path);
});

const out = resolve(__dirname, "..", "src", "data", "python-api.json");
mkdirSync(dirname(out), { recursive: true });
writeFileSync(out, JSON.stringify({ generatedAt: new Date().toISOString(), items }, null, 2));
console.error(`[python] wrote ${out} (${items.length} items)`);
