# gOS single-cell Python exporters

Moved from srctools `inst/gos_single_cell/` on 2026-10-18 (copy commits name the source srctools commit). Each script is
run as `python3 <name>.py ...` with the gOS data directory (`gOS_GBM_sc/data`) given as argument or env var; see each header
for arguments. From R: `sc_run_python_exporter("gos_sc_<name>", args)`. The pipeline needs pysam/pandas/numpy/scipy where noted
For numerical exporters/tests, the cluster's `fishhook_scrna` Python provides numpy/pandas/scipy. The `mskilab_ne1`
Python does not currently provide scipy. Pileup exporters additionally require pysam; choose a Python containing the
needed dependencies explicitly via `GOS_PYTHON` or the wrapper's `python=` argument.

Tests (synthetic fixtures only): `python3 -m unittest discover inst/python/gos_single_cell/tests`.
Verified on 2026-10-19: `/gpfs/commons/groups/imielinski_lab/Software/miniforge3/envs/fishhook_scrna/bin/python -m unittest discover inst/python/gos_single_cell/tests` — 72 tests pass (2026-10-20: 109 with test_snv_columnar).

`gos_sc_allsnv_pileup.py` and `gos_sc_hotspot_rescue.py` still read `precompute/steps/site_evidence.py` and
`precompute/config/regions/gbm.tsv` from srctools; set `SRCTOOLS_GOS_SC` to override that location.

| script | what it does (header) | writes |
|---|---|---|
| `gos_sc_allsnv.py` | gos_sc_allsnv: driver events from ALL called SNVs in driver genes, not only the tree sites. | data/<P>/driver_all_snvs.json, filtered.events.json (staged MODEs; WORK dir) |
| `gos_sc_allsnv_pileup.py` | Per-cell read pileups at the candidate driver-gene SNVs of gos_sc_allsnv (read evidence for sites that are not tree sites). | WORK pileup tables (region_calls-like), not a gOS file directly |
| `gos_sc_ase.py` | Allele-specific expression (ASE) against allelic copy number in the same single cells (G&T-seq: DNA + RNA). | data/<P>/rna/ase.json (+ cohort allelic.json) |
| `gos_sc_clinical.py` | Parse the clinical assay reports (MSK-style Oncopanel and Genexus free text) of the single-cell patients. | data/_cohort/clinical_assays.json |
| `gos_sc_clinical_annot.py` | Local ClinVar (GRCh38 VCF) and gnomAD r2.1 genomes (GRCh37) lookups for clinical SNV / indel sites. No web queries: patient | ClinVar/gnomAD annotation columns for clinical_sites (local lookups) |
| `gos_sc_clinical_events.py` | Clinical-assay alterations as filtered.events-style rows, so the gOS Drivers views can show every clinical finding. | data/<P>/filtered.events.json rows, data/_cohort/clinical_events.json |
| `gos_sc_clinical_recovery.py` | Recovery of clinical-assay alterations (Oncopanel / Genexus) in the single-cell data, real counts only. | data/_cohort/clinical_recovery.json |
| `gos_sc_clinical_sites.py` | HGVS cDNA (c.) -> hg38 genomic coordinates for the clinical-assay variants, to look them up in the single-cell data. | data/_cohort/clinical_sites.json |
| `gos_sc_cohort_event_assoc.py` | Cross-patient associations of driver / copy-number events with cell state and RNA programs. | data/_cohort/event_associations.json |
| `gos_sc_cohort_evolution.py` | Cohort evolution exports for the gOS single-cell cohort views (#19 aligned trees, #21 associations). | data/_cohort/aligned_trees.json, associations.json |
| `gos_sc_complex_twins.py` | Mark complex SV events that the gGnome classifier reported twice under different names. | data/<P>/filtered.events.json (twin marks) |
| `gos_sc_driver_evidence.py` | RNA evidence for the DNA fusions of single-cell patients, for the gOS driver deep-dive. | data/<P>/rna/evidence.json (fusion RNA evidence) |
| `gos_sc_ecdna_plasticity.py` | ecDNA copy number as a source of plasticity. | data/_cohort ecDNA plasticity json (see script header) |
| `gos_sc_event_families.py` | Junction identity, per-cell read support and junction-keyed event families (stage 1) for the gOS single-cell export. | data/<P>/ junction identity / event families json (see header) |
| `gos_sc_homdel_qc.py` | Size of the copy-number-0 segment behind each HOMDEL event, per patient, written onto the event. | data/<P>/filtered.events.json (HOMDEL segment size), broad_cn0.json |
| `gos_sc_hotspot_rescue.py` | Hotspot rescue: configured hotspot drivers (TERT promoter C228T / C250T, IDH1 R132, H3-3A K27M / G34, BRAF V600E) | data/<P>/filtered.events.json, hotspot_rescue.json |
| `gos_sc_pooled_drivers.py` | Import high-confidence pooled-scan driver calls (srctools feature/driver-sensitivity) into data/<P>/filtered.events.json. | data/<P>/filtered.events.json |
| `gos_sc_known_cn_rescue.py` | Tier K ("known biology") rescue of driver SNVs in copy-number-altered regions, with dropout- and CN-aware clonality. | data/<P>/filtered.events.json, <P>/known_cn_rescue.json, _cohort/known_cn_rescue.json (to a staging --out) |
| `gos_sc_portrait.py` | Multi-modal patient portrait + cohort "canonical GBM biology vs novel" export for the gOS single-cell cohort. | data/<P>/ portrait json and _cohort canonical-vs-novel (see header) |
| `gos_sc_rna_covariate_de.py` | Cross-patient RNA differential on the per-cell covariates (rna/covariates.json) -> _cohort/rna_covariate_de.json. | data/_cohort/rna_covariate_de.json, rna_covariate_de_cells.json |
| `gos_sc_rna_covariates.py` | Extra per-cell covariates for the cross-patient RNA differential (data/<P>/rna/covariates.json). | data/<P>/rna/covariates.json |
| `gos_sc_site_qc.py` | Export per-site QC scores to <outdir>/<P>/snv_qc.json, keyed by snv_matrix variant id. | <outdir>/<P>/snv_qc.json |
| `gos_sc_snv_columnar.py` | Columnar SNV site sets ("gos-sc-snv-cols/1"): site table + CSR per-cell reads (byte range per branch, gzipped copy); `tree` converts snv_matrix.json, `cn_mapped` exports the CN-aware tree mapping (srctools cn_tree_mapping) and its cancer-gene driver rows (evidence "CN-aware mapped"). | data/<P>/snv_matrix.sites.json + .cells.bin(.gz), cn_mapped_snvs.sites.json + .cells.bin(.gz), filtered.events.json, _cohort/cn_mapped_summary.json |
| `gos_sc_splice_exon_annot.py` | Adds gene exon context to splicing.json files (patient + _cohort) so single-exon / intergenic calls can be filtered. | splicing.json (patient + _cohort), in place |
| `gos_sc_splice_exons.py` | Exon models of the genes in the single-cell splicing exports, so gOS can draw exons under | data/_cohort/splice_exons.json |
| `gos_sc_splice_support.py` | Where each splice alteration comes from: DNA (a splice-site SNV in filtered.events.json), RNA (aberrant splicing in | data/_cohort/splice_support.json |
| `gos_sc_subclonal_drivers.py` | Cohort subclonal drivers and subclonal signature shifts -> <data>/_cohort/subclonal_drivers.json | data/_cohort/subclonal_drivers.json |
| `gos_sc_tert_sv.py` | Mark TERT promoter rearrangements in gOS single-cell filtered.events.json. | data/<P>/filtered.events.json (TERT promoter SV mark) |
| `gos_sc_tmb_assoc.py` | Cohort TMB associations -> analysis/tmb_assoc/out/tmb_associations.json (copied to _cohort/tmb_associations.json). | data/_cohort/tmb_associations.json |
| `gos_sc_tmb_clade.py` | Clade-level TMB and its associations -> _cohort/tmb_clade_assoc.json  (unit = tree clade / branch, not cell). | data/_cohort/tmb_clade_assoc.json |
| `gos_sc_tree_support.py` | Bootstrap support (Felsenstein bootstrap proportion, FBP, and optionally | tree bootstrap support json for Compare Trees (see header) |

`hotspot_verdicts.js` is the node helper used by allsnv and hotspot_rescue (runs gOS `siteEvidence.js`).
