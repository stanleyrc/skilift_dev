#!/usr/bin/env node
// Site-evidence verdicts of one patient, computed by the gOS front end's own code (siteEvidence.js rankFindings),
// so the exporter's "rescue" decision and the app show the same verdict.
//   node hotspot_verdicts.js <gOS_dev src dir or worktree root> <patient data dir>
// Reads <patient>/precompute/{region_calls,site_evidence}.json; prints JSON { siteId: { level, label, summary, checks, stats } }.
// Needs @babel/core and @babel/plugin-transform-modules-commonjs + pirates from the gOS_dev node_modules (siteEvidence.js is ESM).
const fs = require("fs");
const path = require("path");
const root = path.resolve(process.argv[2] || "");
const dir = path.resolve(process.argv[3] || "");
const src = fs.existsSync(path.join(root, "src")) ? path.join(root, "src") : root;
const nm = path.join(src, "..", "node_modules");
const req = (m) => require(require.resolve(m, { paths: [nm] }));
const babel = req("@babel/core");
const { addHook } = req("pirates");
addHook((code, filename) => babel.transformSync(code, { filename, babelrc: false, configFile: false, plugins: [req("@babel/plugin-transform-modules-commonjs")] }).code, {
  exts: [".js"],
  matcher: (f) => f.startsWith(path.join(src, "helpers")),
});
const { rankFindings } = require(path.join(src, "helpers", "singleCell", "siteEvidence.js"));
const calls = JSON.parse(fs.readFileSync(path.join(dir, "precompute", "region_calls.json"), "utf8"));
const ev = JSON.parse(fs.readFileSync(path.join(dir, "precompute", "site_evidence.json"), "utf8"));
const out = {};
rankFindings(calls, ev, new Map()).all.forEach((f) => {
  out[f.id] = {
    level: f.verdict.level,
    label: f.verdict.label,
    summary: f.verdict.summary,
    checks: f.verdict.checks.map((c) => ({ key: c.key, status: c.status, text: c.text })),
    carrierCells: f.stats.carrierCells,
    coveredTumour: f.stats.coveredTumour,
    nTumourCells: f.stats.nTumourCells,
  };
});
process.stdout.write(JSON.stringify(out));
