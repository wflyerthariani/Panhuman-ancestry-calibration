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
| R-06 | **Proposed analysis settings awaiting confirmation** (the smoke run used ClinVar **2026-09-23**; see D-024) (all `status: proposed` in `config/analysis.yaml`). | Chromosomes chr1–22 + chrX (no Y/M) · CADD **v1.7** as the comparator · GENCODE **v50** MANE Select as the coding definition · bootstrap seed 20260926 · ClinVar = latest weekly release, date recorded at fetch · DMS replicates = the assay-reported score (mean where only replicates exist) · missing AVI score = no record for the exact chrom:pos:ref:alt, or a null/NaN score · population mapping = `config/population_groups.yaml`. | `reports/qc/config_status.md` | open |
| R-08 | **DMS score direction and quality review (smoke scope).** 23 BRCA1/TP53 assays need a biologist's decision. Proposals: 14 lower-is-damaging, 3 higher-is-damaging, **2 conflicts**, **4 unknown**. The conflicts are the TP53 assays in p53 wild-type cells under nutlin (`00000068-a-1`, `00001235-a-1`): dominant-negative variants enrich there, while truncating variants cannot act dominant-negatively. So "damaging" needs defining for those assays. | A reviewer follows `reports/review/dms_direction_review.md` and fills `config/dms_curation.tsv`. Stage 6 will not run until every in-scope assay has a decision. | `config/dms_curation.tsv` | **open — review stop** |
| R-09 | **Review workload at `dms_genes` scope.** 1,033 candidate assays across 527 genes. At 5–10 minutes each, that is about 85–170 hours. | Options for the team: (a) first drop genes that cannot reach ≥20 population-observed variants (this does not depend on direction, so it could run before review); (b) accept proposals where calibration and controls agree (`high` confidence), with a spot-check sample, and review the rest in full; (c) split the review across several people by gene family. | — | open — needs a decision before scaling up |
| R-07 | **Items still pending on data.** These are settled by later stages: AVI score field and direction (Stage 8, needs Atlas access), DMS assay list and per-assay direction and quality (Stage 5 review, R-08), and the DMS gene universe (Stage 9, review stop). | — | `config/analysis.yaml` | pending |

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

### Reference and annotation (Stage 4)

**D-021 — Smoke scope = BRCA1 and TP53 (chr17).** Both have large, well-known DMS datasets.
Their MANE Select CDS covers 32 merged intervals and 6,768 bp. The smoke run stored 38.6 MB of
data in total. Stage 5 will confirm that both have usable MaveDB score sets.

**D-022 — GENCODE v50 basic, MANE Select only.** 19,256 MANE Select transcripts, exactly one per
gene. chrY and chrM rows (including chrY PAR copies) are excluded, matching the chromosome setting
(D-010). One transcript (TMEM247, ENST00000434431.2) has a CDS length that is not a multiple of 3.
It is not a DMS gene and is left as is.

**D-023 — Reference FASTA is stored per contig** (`GRCh38.<contig>.fa.gz`). It is uppercased
(UCSC soft-masks repeats in lower case) and bgzipped, with `.fai` and `.gzi` indexes. Each contig's
length is checked against GRCh38, and its alphabet is checked to be ACGTN only.

**D-024 — ClinVar is fixed to a dated release.** "latest" is read only to learn its
`##fileDate`. Data is then read from the dated file
(`clinvar_20260923.vcf.gz`, checked unchanged before and after reading), so a weekly
replacement cannot corrupt a run. The permanent location is `archive_2.0/<year>/`. The release
used in the smoke run is **2026-09-23**. It must be written into `config/analysis.yaml →
clinvar.release` when the config is frozen (R-06). ClinVar names chromosomes `1..22, X`; they are
renamed to `chr1..` on extraction. At smoke scope: 15,315 records, 11,046 of them SNVs.

**D-025 — AlphaMissense subset for protein-to-genome lookup.** 15,032 rows at smoke scope.
Every REF allele matches GRCh38. Its transcripts are older Ensembl versions (for example
`ENST00000335137.4`), not necessarily MANE Select, so Stage 6 must match on protein change *and*
check the transcript. Its licence is CC BY-NC-SA 4.0 (L-07).

**D-026 — Remote tabix indexes are fetched into the data store** and passed with htslib's
`URL##idx##local.tbi` syntax. Otherwise htslib silently downloads `.tbi` files into the current
working directory. No stray index files were found after Stage 4.

**D-027 — Outputs are deterministic.** Re-running a stage must give byte-identical tables and
reports. Ties in sort order are broken explicitly (a tie in the Stage 3 exclusion table was found
and fixed), and reports do not contain run dates; run times are in the manifest.

**D-028 — Stage 4 depends on Stage 5 at `dms_genes` scope.** That scope takes its genes from
`tables/dms_candidate_inventory.tsv`. The inventory is scope-independent (all human MaveDB
metadata), so the one made during the smoke run is reused. If it is missing, Stage 4 reports
`blocked` with instructions rather than failing.

### DMS candidate inventory (Stage 5)

**D-029 — MaveDB snapshot.** There are 1,201 published human score sets (MaveDB API
2026.2.7.3, taken 2026-09-27). Full records, the mapped-gene table and each downloaded score
CSV are saved under `data_store/dms/mavedb/` and recorded in the manifest. Re-runs reuse the
snapshot, so the inventory does not change silently as MaveDB grows; `ancal stage5 --refresh`
takes a new one. Snapshots are gzip-written with a fixed timestamp, so identical content gives
identical checksums.

**D-030 — Each score set is resolved to one MANE Select gene**, in this order: MaveDB's own
GRCh38 mapping (527), an Ensembl gene ID (42), UniProt accession to primary gene name (434,
UniProt release 2026_03), then the first word of the target name (165). 31 are unresolved,
and 2 have multi-gene targets. When methods disagree, the record is flagged
`gene_resolution_conflict`; there were none in this snapshot. UniProt was added as a source
because many targets are named informally (e.g. "p53", "alpha-synuclein").

**D-031 — The automatic pre-screen applies only exclusions the protocol requires outright.**
168 are excluded:
- 98 meta-analyses: combined or classifier scores built from other score sets, not direct
  measurements. The underlying assays stay in.
- 41 noncoding targets: regulatory elements, which Step 32 keeps as a separate benchmark.
- 17 unresolved genes, 7 with fewer than 20 variants, 3 genes outside MANE Select, and
  2 multi-gene panels.
- None superseded: MaveDB records supersession only on the newer record, and the code handles
  that.

**1,033 candidates across 527 genes remain.** Every record, including excluded ones, is in
`tables/dms_candidate_inventory.tsv` with its reason.

**D-032 — Direction evidence, gathered for in-scope candidates only.**
- **A. Calibrations.** Normal and abnormal ranges from MaveDB score calibrations. Only 82 of
  1,201 score sets have one.
- **B. Built-in controls.** Truncating variants (nonsense *and* frameshift) against synonymous
  variants, measured as AUC. |AUC − 0.5| ≥ 0.25 counts as a clear direction, and at least 10 of
  each class are needed.
  - If there are no synonymous variants, truncating variants are compared with missense
    instead. This is marked *weaker* and lowers confidence.
  - Protein HGVS is read from `hgvs_pro`, or from a submitter column such as
    `HGVS(protein)`, with any `NP_…:` prefix removed.
- **C. Author text.** Excerpts are shown verbatim and never parsed into a decision.
- **Confidence levels.** `high` means A and B agree with B against synonymous. `medium` or
  `low` means one source, or the weaker B. `conflict` or `unknown` means a person must decide.
- **Development check:** the first version counted only nonsense variants and missed
  frameshifts and custom protein columns. 13 of 23 assays then had no control evidence. After
  the fix, 19 of 23 have a proposal.

**D-033 — The review sheet is a TSV (`config/dms_curation.tsv`), not the YAML in the plan.**
Biologists can edit it in Excel or Google Sheets. Evidence columns are regenerated on every run.
The five reviewer columns (`reviewer_decision`, `reviewer_direction`, `reviewer_notes`,
`reviewer_name`, `reviewed_on`) are always preserved. MaveDB free text is stripped of line breaks
and tabs, because they were breaking TSV rows. The reviewer's guide is
`reports/review/dms_direction_review.md`: numbered steps, background, one evidence card per
assay, and sibling score sets from the same experiment (to avoid counting one measurement
twice).

**D-034 — Scores are downloaded only for in-scope candidates.** At smoke scope that is 23
files, 9.9 MB. Downloading all 1,033 candidates at `dms_genes` scope is estimated at about
**0.25 GB** (3.6 M variants at about 69 bytes each).

**D-035 — Stage 6 gate.** DMS harmonization will refuse to run until every in-scope candidate has
a `reviewer_decision`, and every included one has a `reviewer_direction`. Scores whose raw
direction is `lower_is_more_damaging` are reversed there. The original score is always kept
(Step 13).

---

## 3. Known limitations (inputs to item 15)

| ID | Limitation | Effect | Mitigation |
|---|---|---|---|
| L-01 | **Storage access not yet arranged.** Only 1 GB of local data is allowed (D-002). | Items 7, 8, 10, 12, 13 and 14 can only be produced at `smoke` scope for now. | `SCALING_UP.md` gives the config-only switch to `dms_genes` and `exome`. |
| L-02 | **AVI (AlphaGenome Atlas) access and format unknown.** | Stage 8 (item 8) is `blocked` until the source is set. | `SCALING_UP.md` §3. |
| L-03 | **No per-population callability mask** is distributed with the HGDP+1KG callset. | "Callable" cannot be tested directly (Step 7). | AN > 0 is used as a proxy, and this is documented. |
| L-04 | **CADD is licensed for non-commercial use only.** | This may restrict use of the comparator in a commercial context. | Confirm the licence position with the team. |
| L-05 | **Per-dataset licences of MaveDB score sets vary.** Of 1,201: 1,170 CC0, 14 CC BY 4.0, 11 CC BY-SA 4.0 and **6 CC BY-NC-SA 4.0**. | The NC-SA sets may not be usable commercially or redistributable. | Each licence is recorded in the inventory; NC sets are flagged `licence_CC_BY-NC-SA_4.0` on review cards. |
| L-06 | **Hail Table decoder covers one layout only** (D-016). | A future re-release in a different format would stop Stage 3. | The decoder fails loudly; the fallback is to read the table with Hail. |
| L-07 | **AlphaMissense is licensed CC BY-NC-SA 4.0 (non-commercial, share-alike).** | This may restrict commercial use of anything derived from the protein-to-genome lookup. AVI itself includes AlphaMissense, so the AVI licence needs checking too (L-02). | Confirm the licence position with the team, alongside CADD (L-04). |
