# Scaling up: from the 1 GB smoke run to the full analysis

Right now the pipeline runs at `scope: smoke` with a **1 GB local data cap**. This
proves every stage end-to-end on a few genes (`make d1`). When more storage is available,
follow these steps in order. **Only config changes are needed, not code changes.**

## 0. Using files you already have

Every `uri:` in `config/data_sources.yaml` accepts a local path in place of a URL.

- **gnomAD HGDP+1KG VCFs**: keep `{contig}` in the path; the `.tbi` must sit beside each file.
  Local files are read with `tabix`, which is much faster than the remote reader.
- **Reference**: either one file per chromosome (`.../{contig}.fa.gz`) or a single genome FASTA
  with a `.fai` index (no `{contig}` in the path).
- **ClinVar**: a local VCF is used as it is. Its date comes from `release:` (YYYYMMDD), or from
  the file's `##fileDate` header when `release: latest`.
- **GENCODE, AlphaMissense**: the path to the downloaded file.
- **Sample metadata**: local copies work, but the files are small (6 MB) and the URLs are fine.
- **MaveDB and UniProt** are live APIs and always need internet access.

If a local copy is a different version from the one named in the config, record that in
`docs/DECISIONS.md`: the Stage 3 and Stage 4 gates check counts and reference alleles against
the named versions and will fail on a mismatch.

> Sizes and runtimes below are estimates. The budget guard (`make budget`) checks
> actual usage before and after every stage.

---

## 1. Decide where the data lives

Edit `config/data_sources.yaml`:

```yaml
storage:
  root: /Volumes/BigDisk/ancal_data   # or a mounted cloud disk; any local path works
  local_budget_gb: 50                 # set to what you can actually spare
```

- **Existing smoke data**: to reuse it, move the contents of `./data_store/` into the new
  root. Otherwise it is simply refetched.
- **Check**: run `make budget` and confirm it shows the new root and budget.

## 2. Move to `dms_genes` scope (needed to finish deliverable 1)

**Needs:** about 5–15 GB. Most of this is the full GRCh38 FASTA (about 3 GB) plus gnomAD
HGDP+1KG sites across the coding exons of every candidate DMS gene. Genotypes are
never stored.

```sh
make clean-scope SCOPE=smoke     # optional: archives smoke outputs to reports/archive/smoke_<date>/
# set `scope: dms_genes` in config/data_sources.yaml
make stage4     # regions BED (genes from the Stage 5 inventory made during the smoke run),
                # reference for every contig in scope, AlphaMissense + ClinVar subsets
make stage5     # DMS inventory + scores for all candidate genes
                # STOP: biologist review of config/dms_curation.tsv
                # (guide: reports/review/dms_direction_review.md; see R-09 on workload)
make stage6     # DMS harmonization (blocked until the review is done, unless
                # review_gates.dms_review is `provisional`)
make stage7     # gnomAD sites; the longest step, resumable per chunk of 20 intervals
make stage8     # AVI join and coverage (needs AVI access, see section 3)
make stage9     # proposed DMS gene universe
                # STOP: team approves tables/proposed_dms_gene_universe.tsv, then:
ancal freeze-config
make stage10 stage11
```

`make d1` runs the same sequence in one command; it does not stop for the two reviews, so use
it with `review_gates.dms_review: provisional` for a first pass and re-run it after each review.

- **Check**: `make status` shows every stage `pass`, or `blocked` for a documented reason.
  `deliverables/d1/README.md` lists every item with `scope=dms_genes`.

## 3. When AVI access arrives

Fill in the `avi:` block in `config/data_sources.yaml`:

```yaml
avi:
  uri: <path or URL to the AlphaGenome Atlas AVI file>
  format: <tsv | parquet | vcf | ...>
  score_field: <column name>
  direction: higher_is_more_damaging   # or lower_is_more_damaging; the stage verifies and documents it
  # optional:
  version: <Atlas release>
  columns: {chrom: ..., pos: ..., ref: ..., alt: ..., transcript: ...}   # if not CHROM/POS/REF/ALT
  damaging_threshold: <number>         # only if the Atlas defines one
```

- **Formats**: `parquet` (file, directory or glob), `tsv`/`csv` (plain or gzipped) and `vcf`.
- **Region access**: a `tsv` or `vcf` with a tabix index beside it is read by region. Without an
  index a `tsv`/`csv` is streamed once in full and filtered, which is slow for a genome-wide file.
- **If `direction` is left null**, Stage 8 infers it and asks for confirmation in its report.
- **If a column cannot be found**, Stage 8 is `blocked` and its message lists the file's columns.

Then run:

```sh
make stage8 && make stage9 stage10 stage11
```

- **If the AVI file is genome-wide and too large for the budget**, Stage 8 stays `blocked` with
  the size it needs. It can still run if the file is tabix-indexed (read by region),
  or if Stage 8 is run on the machine where the file lives.

## 4. `exome` scope (full ClinVar benchmark, later)

**Needs:** tens of GB. Same as section 2, but with `scope: exome`, starting from `make stage4`.
Two tables grow large at this scope and have not been run at it: the coding-consequence table
(about 100 million rows) and the population long table (about 1.7 million sites × 176 levels).
Expect to need tens of GB of memory for Stages 6, 7 and 11.

---

## Stage-by-stage reference

The gate checks and common failures for each stage. Runtimes are for `smoke` scope.

| Stage | Gate checks | If it fails |
|---|---|---|
| 0 | bcftools/tabix/bgzip/samtools/snakemake on PATH; Python imports; `pytest`; disk budget | Run `make env`; read `logs/stage0.log` |
| 1 | `config/analysis.yaml` and `config/data_sources.yaml` pass the schema (AF bins contiguous, joins on chrom/pos/ref/alt, seven group labels, every source has a URI except AVI); writes `reports/qc/config_status.md` | Read the validation message in `logs/stage1.log`; fix the YAML, never the generated report |
| 2 | Software manifest written with bcftools/tabix/samtools/snakemake found; one live download (gnomAD `post_qc_summary.tsv`, about 10 KB) recorded with source and output SHA-256; every manifest row has all Step 1 fields | Network: check that `storage.googleapis.com` is reachable. Incomplete row: `logs/stage2.log` names the row and fields |
| 3 | Every post-QC sample has one population label and one group label; 4,094 samples and 80 populations; per-population counts equal Koenig et al.'s published post-QC counts; unrelated/related counts equal their released PCA sets (3,400 / 694); every post-QC sample is in the VCF header; the VCF contig length matches GRCh38 | A count mismatch means the gnomAD release files changed: compare the ETags in `reports/manifests/data_manifest.tsv` and read `reports/qc/sample_count_report.md`. The related-ID decoder refuses unfamiliar Hail formats; if so, read the table with Hail instead. Stage 3 does not depend on scope (sample metadata is always fetched in full, about 7 MB) |
| 4 | 19,000+ MANE Select transcripts, one per gene, with CDS lengths that are multiples of 3; every scope gene found; each reference contig has the GRCh38 length and only ACGTN; every AlphaMissense and ClinVar REF allele equals the GRCh38 base; ClinVar comes from a dated release that did not change during the read; disk budget | `blocked` at `dms_genes` scope if `tables/dms_candidate_inventory.tsv` is missing: run Stage 5 first. A REF mismatch means a build or contig-naming problem: stop and investigate, do not filter it away. About 3 minutes at smoke scope, most of it streaming AlphaMissense (643 MB) |
| 5 | MaveDB snapshot taken or reused; every published human score set classified (candidate or excluded, with a reason); every scope gene has at least one candidate assay; disk budget. The stage **passes** before review, and its status note shows how many in-scope assays are reviewed | A MaveDB outage gives HTTP errors after 5 retries, so run it again later. To pick up new MaveDB records, run `ancal stage5 --refresh` (this changes the inventory, so note it in `docs/DECISIONS.md`). About 1 minute for a fresh snapshot |
| 6 | Review gate (see `review_gates.dms_review`); coding-consequence table built; its protein changes agree with AlphaMissense for at least 98% of shared SNVs; every included assay has a score file; no duplicate (variant, assay) rows; percentiles in (0, 1] | `blocked` in `strict` mode until `config/dms_curation.tsv` is fully reviewed. An assay listed as `assay_numbering_mismatch` does not use MANE Select residue numbering and no single offset fits: the reviewer should check its construct. Seconds |
| 7 | All post-QC samples found in the VCF; long table = variants × 176 levels; no allele number above 2 × samples; population counts sum to group and overall counts; REF equals GRCh38 for at least 99.9% of variants | A chunk that fails 4 times stops the stage with the bcftools or network error; re-running resumes from the finished chunks in `data_store/variants/chunks/`. About 1 minute for remote files |
| 8 | AVI records found in scope; AVI REF equals GRCh38 for at least 99% of coding SNV records; score direction agrees with the two sign checks; every population variant has one join class; coverage within [0, 1] | `blocked` while `sources.avi.uri` is null or a column is missing. A direction failure means `avi.direction` contradicts the evidence in `reports/qc/avi_join_report.md` |
| 9 | Every candidate gene in scope has a status; every proposed exclusion has a reason | `blocked` if Stage 6 has not run. The result is a proposal, never a freeze |
| 10 | Every ClinVar record classified; no VUS in the eligible set; MANE Select missense agrees with ClinVar's consequence for at least 95% of coding SNVs; no eligible variant listed twice | `blocked` if the Stage 4 ClinVar subset or the Stage 6 consequence table is missing |
| 11 | Counts never increase along the attrition funnel; every lost variant has a reason; seven groups in the forecast; every earlier stage is `pass` or `blocked`; every passed item has its files | A `fail` here usually means an earlier stage failed or was never run: check `make status` |
