# Decision log and open review items

Everything of importance decided, discovered or deferred while building deliverable 1
(the data-readiness package) is recorded here, together with the evidence. Generated
reports (`reports/qc/*.md`) carry the numbers. This file carries the reasoning.

- **Open review items** (section 1) need a decision from Arman, Caroline or the team.
  When one is settled, record the outcome and date, and update the matching
  `status:` in `config/analysis.yaml` if it has one.
- **Decisions** (section 2) are numbered `D-NNN` and never renumbered. A reversed decision
  is marked *superseded*, with a pointer to the new one.
- **Known limitations** (section 3) feed the unresolved-issues list (deliverable item 15),
  which Stage 11 generates.

Protocol source: Caroline Duncan's Slack canvas "Ancestry calibration audit"
(F0C3D5MUBBL, #4-scientific-outputs, created 2026-09-21). Background: the
"Preprint Intro" canvas (F0C3RGDFBC4).

---

## 1. Open review items

| ID | Item | Proposed resolution | Where | Status |
|---|---|---|---|---|
| R-01 | **Which sample QC to use.** gnomAD's `high_quality` flag removes every Biaka, Mbuti, San, Bougainville, PapuanHighlands and PapuanSepik sample, which deletes the OCE group. | Use the Koenig et al. (2024) post-QC release: 4,094 samples in 80 populations (D-013). | `reports/qc/sample_count_report.md` | open — needs sign-off |
| R-02 | **Unrelated count: 3,400 or 3,376.** Koenig's released unrelated/related PCA sets are 3,400 / 694, and our derivation matches them exactly. The `n_unrelated` column of `post_qc_summary.tsv` sums to 3,376, with 25 populations differing in both directions. | Use 3,400 (D-015). Ask the Koenig et al. authors, or check the paper's supplement, to find which relatedness pass that column comes from. | `tables/population_counts.tsv` (`diff_unrelated_vs_koenig_summary`) | open |
| R-03 | **List of recently admixed populations.** ACB, ASW, CLM, MXL, PEL and PUR are flagged. None are excluded (Step 5). | Keep this list. Decide whether any HGDP populations (for example Hazara or Uygur) should also be flagged. | `tables/population_mapping.tsv` | open |
| R-04 | **Two unexplained exclusions.** HGDP01371 (Basque) and LP6005441-DNA_A09 (Naxi) pass gnomAD's hard filters and are not PCA outliers, but are absent from Koenig's post-QC release. The reason is not documented. | Exclude them, following the release. Record the reason as `not_in_koenig_post_qc_release`. | `tables/sample_metadata.tsv` | open |
| R-05 | **OCE may miss the inferential criteria.** OCE has 30 post-QC samples, 27 unrelated, across 3 populations. It may not reach ≥200 evaluable variants and ≥5 genes (Step 18). | Stage 11 will estimate this at `dms_genes` scope. If OCE fails, it stays descriptive only, as the protocol already specifies. | Stage 11 | watch |
| R-06 | **Proposed analysis settings awaiting confirmation** (all `status: proposed` in `config/analysis.yaml`). | Chromosomes chr1–22 + chrX (no Y/M) · CADD **v1.7** as the comparator · GENCODE **v50** MANE Select as the coding definition · bootstrap seed 20260926 · ClinVar = latest weekly release, date recorded at fetch · DMS replicates = the assay-reported score (mean where only replicates exist) · missing AVI score = no record for the exact chrom:pos:ref:alt, or a null/NaN score · population mapping = `config/population_groups.yaml`. | `reports/qc/config_status.md` | open |
| R-07 | **Items still pending on data.** These are settled by later stages: AVI score field and direction (Stage 8, needs Atlas access), DMS assay list and per-assay direction and quality (Stage 5, review stop), and the DMS gene universe (Stage 9, review stop). | — | `config/analysis.yaml` | pending |

---

## 2. Decisions

### Infrastructure (Stages 0–2)

**D-001 — Staged pipeline with gates.** Each stage is a Snakemake target that writes
`reports/status/stageN.json` (`pass` / `fail` / `blocked`) and exits non-zero on `fail`.
`blocked` (missing access or budget) is not a failure; it feeds item 15. Snakemake runs
with `--keep-incomplete` so that a failed stage's status file is not deleted.

**D-002 — 1 GB local data budget** (Arman, 2026-09-26). The budget covers `data_store/` plus the
generated output directories. `ancal budget-check` enforces it before and after every stage and
while streaming. Large sources are read by remote tabix range queries or streamed through a
filter; they are never stored whole. Scope (`smoke` → `dms_genes` → `exome`) is a config
switch. See `SCALING_UP.md` for moving to more storage.

**D-003 — Conda (bioconda), not uv.** bcftools, htslib and samtools are only packaged there.
The `ancal` environment lives in `~/miniforge3/envs/ancal` and is **not** counted against the
data budget (about 0.8 GB). There is no Hail or Spark. `zstandard` was added in Stage 3 (D-016).

**D-004 — Manifests are tracked in git** (`reports/manifests/`), not in the git-ignored data
store, so the provenance record (Step 1) is version-controlled. Small generated outputs
(`tables/`, `reports/`, `deliverables/`) are also tracked. Data, logs and region BEDs are not.

**D-005 — Provenance for every file.** Each manifest row records: source URL, ETag,
Last-Modified, size, SHA-256 of the source, SHA-256 of the output, command, software and
version, code commit, date, operator, stage and scope. Streamed sources are hashed as they
pass, so the full-file checksum is recorded even though the file is never stored.
Remote-tabix outputs record the command line instead of a source hash (`not_hashed:remote_region_query`).
`fetch()` and `record()` are idempotent. Code is committed **before** each stage runs, so
manifest rows reference a clean commit.

**D-006 — Source locations verified 2026-09-26** (HTTP HEAD):
- **gnomAD v3.1.2 HGDP+1KG dense VCFs:** per chromosome, **about 100–280 GB each** (chr17 is 108 GB),
  tabix-indexed. They are always read by region.
- **Sample metadata:** the correct path is `hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv`,
  dated 2024-05-01. My first guess had the `hgdp_1kg` and `hgdp_1kg_v2` folders swapped.
- **CADD v1.7 GRCh38 SNVs:** 87 GB. **AlphaMissense hg38:** 643 MB. **ClinVar GRCh38 VCF:** 198 MB
  (weekly; last modified 2026-09-24).
- **MaveDB API:** version 2026.2.7.3. **GENCODE:** v50 is the latest release (v51 returns 404).

**D-007 — Reference FASTA from UCSC per-chromosome files** (`hg38/chromosomes/{contig}.fa.gz`),
so that only in-scope contigs are downloaded (chr17 is about 25 MB). The files are uppercased
and bgzipped locally (Stage 4).

**D-008 — AlphaMissense is used only as a lookup from protein changes to genomic
coordinates.** AlphaMissense is a component of AVI, so this use is logged as a
leakage consideration for Step 35.

**D-009 — Analysis configuration.** Values the protocol fixes are written in directly and
checked by the schema: allele-frequency bins (checked for no gaps or overlaps), the 200 / 5 /
0.30 eligibility criteria, 2,000 hierarchical bootstrap iterations, joins on
chrom/pos/ref/alt only (never rsID), `exclude_admixed: false`, HPRC excluded from the primary
analysis, and the seven fixed group labels. Anything open is a `pending` or `proposed` block.
`ancal freeze-config` refuses to run until every block is `confirmed`, and never overwrites an
existing frozen file.

**D-010 — chrX is included.** Its allele numbers must come from genotypes (hemizygous males),
not from 2 × N.

**D-011 — CADD licence.** CADD is free for non-commercial use only. This is listed for item 15.
CADD is not fetched in deliverable 1.

**D-012 — LF line endings in all written TSVs.** Python's `csv` module defaults to CRLF; this
was fixed in `c6ab6c8`.

### Samples and populations (Stage 3)

**D-013 — Sample set = Koenig et al. (2024) post-QC release.** This is defined as the samples
in `pca/pc_scores_without_outliers/GLOBAL_scores_without_outliers.txt.bgz`: **4,094 samples
in 80 populations**, with per-population counts equal to `post_qc_summary.tsv` for all 80.
Evidence:
- **Why not gnomAD's flag:** gnomAD's `high_quality` flag marks 3,942 and drops six
  divergent populations (R-01).
- **Filters versus the release:** "not hard-filtered and not a PCA outlier" gives 4,096. The
  release differs by 2 samples (R-04).
- **Callset:** the VCF header has 4,151 samples, 4,150 of them in the metadata. The extra one is
  `CHMI_CHMI3_WGS2`, the CHM13 synthetic-diploid control, which is excluded.
- **Exclusions:** 57 in total. 31 gnomAD hard filters (24 sex aneuploidy, 6 ambiguous sex,
  1 bad QC metrics), 23 PCA outliers, 2 not in the release, and 1 missing from the metadata.

**D-014 — Population label = the `population` column (80 labels).** gnomAD's
`hgdp_tgp_meta.Population` has 78 labels. It merges NorthernHan into Han and PapuanHighlands
and PapuanSepik into Papuan, and uses older names for Biaka, Mbuti, Mongolian, BergamoItalian
and Bougainville. Both labels are kept in `tables/sample_metadata.tsv`. The seven-group label
is `hgdp_tgp_meta.Genetic.region`, which matches the protocol's AFR, AMR, CSA, EAS, EUR, MID
and OCE exactly.

**D-015 — Unrelated set = post-QC samples not in Koenig's `related_sample_ids` table.** This
gives 3,400 unrelated (primary) and 694 related (sensitivity only), equal to the sample
counts of their released `unrelateds_without_outliers` and `relateds_without_outliers` PCA
MatrixTables. The published per-population `n_unrelated` disagrees (R-02). gnomAD's own
relatedness flags are **not** used: they were computed across all of gnomAD and disagree
with Koenig in 41 populations.

**D-016 — The related-sample Hail Table is read without Hail.** A strict decoder
(`src/ancal/sources/hail_table.py`) reads exactly the one layout it has been verified
against (string key, LEB128 + zstd blocks, Hail 0.2.130). It refuses anything else and checks
every partition's row count against the table's metadata. This avoids adding about 700 MB of
Hail, Spark and Java. It decoded 698 IDs, of which 694 are post-QC.

**D-017 — Principal components = Koenig's 20 global PCs** (post-QC, without outliers). They
are kept for every post-QC sample for Step 28. Koenig also publishes per-group PC files
(AFR, AMR and so on), which are available if Step 28 needs within-group structure.

**D-018 — PC discordance is a review aid only.** Each post-QC sample's nearest group centroid
(median of PCs 1–10) is compared with its assigned group. 173 samples are discordant, mostly
AMR (108; PUR alone accounts for 81) and CSA (43). Nobody is excluded on this basis (Step 5).

**D-019 — Sex = `sex_imputation.sex_karyotype`.** It is available for every sample
(XX 1,904, XY 2,216, plus X, XXY, XYY and ambiguous, which are mostly hard-filtered). The
self-reported `sex` column exists for only 921 samples and is kept, but not used.

**D-020 — Genome build confirmed.** The VCF header's chr17 length (83,257,441) equals GRCh38.

---

## 3. Known limitations (inputs to item 15)

| ID | Limitation | Effect | Mitigation |
|---|---|---|---|
| L-01 | **Storage access not yet arranged.** Only 1 GB of local data is allowed (D-002). | Items 7, 8, 10, 12, 13 and 14 can only be produced at `smoke` scope for now. | `SCALING_UP.md` gives the config-only switch to `dms_genes` and `exome`. |
| L-02 | **AVI (AlphaGenome Atlas) access and format unknown.** | Stage 8 (item 8) is `blocked` until the source is set. | `SCALING_UP.md` §3. |
| L-03 | **No per-population callability mask** is distributed with the HGDP+1KG callset. | "Callable" cannot be tested directly (Step 7). | AN > 0 is used as a proxy, and this is documented. |
| L-04 | **CADD is licensed for non-commercial use only.** | This may restrict use of the comparator in a commercial context. | Confirm the licence position with the team. |
| L-05 | **Per-dataset licences of MaveDB score sets vary.** | Some DMS sets may not be redistributable. | Stage 5 records each score set's licence. |
| L-06 | **Hail Table decoder covers one layout only** (D-016). | A future re-release in a different format would stop Stage 3. | The decoder fails loudly; the fallback is to read the table with Hail. |
