// Regenerates figure_exports_fixture.json (synthetic inputs) and figure_exports_golden.json (what the gOS browser code computes on them),
// the reference for tests/test_gos_sc_figure_exports.py.
// usage (node 18): GOS_SRC=<gOS_dev>/src/helpers/singleCell SETTINGS=<gOS_dev>/public/settings.json \
//        node --experimental-specifier-resolution=node make_figure_golden.mjs <out_dir>
// (run from a directory whose package.json has {"type":"module"} and which can import $GOS_SRC; see the test README comment)
import fs from "fs";
import path from "path";
import { pathToFileURL } from "url";

const SRC = process.env.GOS_SRC;
const imp = (f) => import(pathToFileURL(path.join(SRC, f)).href);
const { layoutTree, parseNewick } = await imp("newick.js");
const { phyloSignal, walkGroups, columnCorrelations, isNormalClone } = await imp("figures.js");
const { binAt, buildBinIndex } = await imp("matrix.js");
const out = process.argv[2];

// --- synthetic fixture (own PRNG, deterministic) ---
let seed = 12345;
const rnd = () => ((seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0) / 4294967296);
const nCells = 16;
const ids = Array.from({ length: nCells }, (_, i) => `FX_cell${String(i + 1).padStart(2, "0")}`);
const clone = ids.map((_, i) => (i === 3 || i === 11 ? "Normal" : i < 8 ? "Clone1" : "Clone2"));
const newick = `((((${ids[0]}:0.1,${ids[1]}:0.2):0.3,(${ids[2]}:0.15,${ids[3]}:0.25):0.2):0.1,((${ids[4]}:0.2,${ids[5]}:0.1,${ids[6]}:0.3):0.25,(${ids[7]}:0.1,${ids[8]}:0.2):0.4):0.2,((${ids[9]}:0.3,${ids[10]}:0.2):0.1,(${ids[11]}:0.2,(${ids[12]}:0.1,${ids[13]}:0.1,${ids[14]}:0.2):0.3):0.15,${ids[15]}:0.5):0.15);`.replace(");", "));");
const datafiles = ids.map((id, i) => ({
  pair: id, entry_type: "cell", patient_id: "FXP", clone_id: clone[i],
  ploidy: Math.round((1.8 + rnd() * 0.6) * 100) / 100, qc_depth: Math.round((6 + rnd() * 5) * 1000) / 1000, snv_count: Math.round(900 + rnd() * 400),
}));
const copies = (lo, hi, p) => Object.fromEntries(ids.map((id, i) => [id, rnd() < p ? Math.round(lo + rnd() * (hi - lo)) : 0]));
const walks = {
  patient: "FXP", cells: ids,
  walks: [
    { id: "1", label: "EGFR", driver_genes: "EGFR", curated: true, nodes: [{ chromosome: "7", start: 54500000, end: 55500000 }, { chromosome: "7", start: 55800000, end: 56200000 }], cells: copies(3, 30, 0.6) },
    { id: "2", label: "MDM2", driver_genes: "MDM2", curated: true, nodes: [{ chromosome: "12", start: 68000000, end: 70000000 }], cells: copies(2, 20, 0.5) },
    { id: "3", label: "EGFR2", driver_genes: "EGFR", curated: true, nodes: [{ chromosome: "7", start: 55000000, end: 55400000 }], cells: copies(1, 10, 0.4) },
  ],
};
const segs = (id) => {
  // total CN segments on chr7 and chr12 (irregular breakpoints), amplified inside the walk loci in "carrier" cells
  const w = walks.walks;
  const out = [];
  const add = (chr, bps, vals) => bps.slice(0, -1).forEach((b, k) => out.push({ chromosome: chr, startPoint: b + (k ? 1 : 0), endPoint: bps[k + 1], iid: out.length + 1, y: vals[k], type: "interval" }));
  const c1 = w[0].cells[id] + w[2].cells[id];
  add("7", [1, 30000000, 54500000, 55500000, 55800000, 56200000, 80000000, 159345973], [2, 3, 2 + c1, 2 + Math.round(c1 / 2), 2 + c1, 3, 2]);
  add("12", [1, 40000000, 68000000, 70000000, 133275309], [2, 2, 2 + w[1].cells[id], 2]);
  return { settings: {}, intervals: out, connections: [] };
};
const complex = Object.fromEntries(ids.map((id) => [id, segs(id)]));
const fixture = { newick, datafiles, walks, complex };
fs.writeFileSync(path.join(out, "figure_exports_fixture.json"), JSON.stringify(fixture));

// --- gOS computations ---
// chromoBins as updateChromoBins
const settings = JSON.parse(fs.readFileSync(process.env.SETTINGS, "utf8"));
const set = settings.coordinates.sets.hg38;
let boundary = 0;
const chromoBins = {};
set.forEach((c) => {
  chromoBins[c.chromosome] = { ...c, length: c.endPoint, startPlace: boundary + c.startPoint, endPlace: boundary + c.endPoint };
  boundary += c.endPoint;
});
// phylo signal as phyloSignalPanel
const layout = layoutTree(parseNewick(newick));
const rec = new Map(datafiles.map((c) => [c.pair, { ...c, cell_id: c.pair }]));
const cellIds = walks.cells;
const tumor = cellIds.filter((id) => rec.has(id) && !isNormalClone(rec.get(id).clone_id));
const groups = walkGroups(walks.walks, cellIds, {});
const CONTROLS = [["Seq. depth", "qc_depth"], ["Ploidy", "ploidy"], ["SNVs / cell", "snv_count"]];
const vars = [
  ...CONTROLS.map(([key, field]) => ({ key, value: (id) => Number(rec.get(id)?.[field]) })),
  ...groups.flatMap((g) => g.walks.map((w) => ({ key: w.label || w.name || `${w.id}`, value: (id) => Math.log1p(Number(w.cells?.[id]) || 0) }))),
];
const phylo = phyloSignal(layout, tumor, vars, { nPerm: 199 });

// triangle correlation, as cohortFiguresPanel.lociDomains + triangleCorrelation.js
const PAD = 1.5e6;
const globalPos = (chr, pos) => chromoBins[`${chr}`.replace(/^chr/, "")].startPlace + Number(pos);
function lociDomains(loci, pad = PAD) {
  const by = new Map();
  loci.forEach(([chr, s, e]) => {
    const key = `${chr}`.replace(/^chr/, "");
    if (!chromoBins?.[key] || !Number.isFinite(s)) return;
    const cur = by.get(key);
    by.set(key, cur ? [Math.min(cur[0], s), Math.max(cur[1], e)] : [s, e]);
  });
  return [...by.entries()]
    .map(([chr, [s, e]]) => {
      const bin = chromoBins[chr];
      const p = Math.max(pad, 0.25 * (e - s));
      const a = globalPos(chr, Math.max(bin.startPoint ?? 1, s - p));
      const b = globalPos(chr, Math.min(bin.endPoint ?? e + p, e + p));
      return [a, b, chr];
    })
    .filter((d) => d[1] > d[0])
    .sort((x, y) => x[0] - y[0]);
}
const walkLoci = (ws) => ws.flatMap((w) => (w.nodes || []).map((n) => [n.chromosome, Number(n.start), Number(n.end)]));
const K = 96;
function triangle(ws) {
  const domains = lociDomains(walkLoci(ws));
  const total = domains.reduce((s, [a, b]) => s + (b - a), 0);
  const positions = [];
  domains.forEach(([a, b], k) => {
    const n = Math.max(2, Math.round((K * (b - a)) / Math.max(1, total)));
    for (let i = 0; i < n; i += 1) positions.push({ g: a + ((i + 0.5) / n) * (b - a), d: k });
  });
  const cn = (id) => ws.reduce((s, w) => s + (Number(w.cells?.[id]) || 0), 0);
  const carriers = tumor.filter((id) => cn(id) >= 1);
  const others = tumor.filter((id) => cn(id) < 1);
  const rowOf = (id) => {
    const iv = complex[id].intervals;
    const binIndex = buildBinIndex({ chromosome: iv.map((d) => d.chromosome), start: iv.map((d) => d.startPoint), end: iv.map((d) => d.endPoint) }, chromoBins);
    return { binIndex, values: Float32Array.from(iv, (d) => Number(d.y)) };
  };
  const cache = new Map(tumor.map((id) => [id, rowOf(id)]));
  const corr = (idl) => {
    const rows = idl.map((id) => cache.get(id));
    const k = positions.length;
    const M = new Float32Array(rows.length * k);
    rows.forEach((row, i) => {
      let sum = 0;
      positions.forEach(({ g }, j) => {
        const b = binAt(row.binIndex, g);
        M[i * k + j] = b >= 0 ? Math.log1p(Math.max(0, row.values[b])) : 0;
        sum += M[i * k + j];
      });
      const mean = sum / Math.max(1, k);
      for (let j = 0; j < k; j += 1) M[i * k + j] -= mean;
    });
    return Array.from(columnCorrelations(M, rows.length, k));
  };
  return {
    k: positions.length, windowOf: positions.map((p) => p.d), carriers, others,
    positionCn: tumor.map((id) => positions.map(({ g }) => { const b = binAt(cache.get(id).binIndex, g); return b >= 0 ? cache.get(id).values[b] : null; })),
    corrCarriers: corr(carriers), corrOthers: corr(others),
  };
}
const tri = { "w:1": triangle([walks.walks[0]]), "w:2": triangle([walks.walks[1]]), "g:EGFR": triangle(groups.find((g) => g.key === "EGFR").walks) };
fs.writeFileSync(path.join(out, "figure_exports_golden.json"), JSON.stringify({ tumor, phylo, triangle: tri, groupKeys: groups.map((g) => g.key) }));
console.log("golden written", phylo.map((p) => `${p.key}: I=${p.I.toFixed(4)} z=${p.z.toFixed(2)} p=${p.p}`).join("\n"));
