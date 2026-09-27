# Scaling up: from the 1 GB smoke run to the full analysis

Right now the pipeline runs at `scope: smoke` with a **1 GB local data cap**. This
proves every stage end-to-end on a few genes. When more storage is available,
follow these steps in order. **Only config changes are needed, not code changes.**

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
make stage6     # DMS harmonization
make stage7     # gnomAD sites; the longest step, resumable per chromosome
make stage8     # AVI join and coverage (needs AVI access, see section 3)
make stage9     # proposed DMS gene universe
                # STOP: team approves tables/proposed_dms_gene_universe.tsv, then:
ancal freeze-config
make stage10 stage11
```

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
```

Then run:

```sh
make stage8 && make stage9 stage10 stage11
```

- **If the AVI file is genome-wide and too large for the budget**, Stage 8 stays `blocked` with
  the size it needs. It can still run if the file is tabix-indexed (read by region),
  or if Stage 8 is run on the machine where the file lives.

## 4. `exome` scope (full ClinVar benchmark, later)

**Needs:** tens of GB. Same as section 2, but with `scope: exome`, starting from `make stage4`.

---

## Stage-by-stage reference

The gate checks, common failures and expected runtimes for each stage are added
below as each stage is implemented.

| Stage | Gate checks | If it fails |
|---|---|---|
| 0 | bcftools/tabix/bgzip/samtools/snakemake on PATH; Python imports; `pytest`; disk budget | Run `make env`; read `logs/stage0.log` |
| 1 | `config/analysis.yaml` and `config/data_sources.yaml` pass the schema (AF bins contiguous, joins on chrom/pos/ref/alt, seven group labels, every source has a URI except AVI); writes `reports/qc/config_status.md` | Read the validation message in `logs/stage1.log`; fix the YAML, never the generated report |
| 2 | Software manifest written with bcftools/tabix/samtools/snakemake found; one live download (gnomAD `post_qc_summary.tsv`, about 10 KB) recorded with source and output SHA-256; every manifest row has all Step 1 fields | Network: check that `storage.googleapis.com` is reachable. Incomplete row: `logs/stage2.log` names the row and fields |
| 3 | Every post-QC sample has one population label and one group label; 4,094 samples and 80 populations; per-population counts equal Koenig et al.'s published post-QC counts; unrelated/related counts equal their released PCA sets (3,400 / 694); every post-QC sample is in the VCF header; the VCF contig length matches GRCh38 | A count mismatch means the gnomAD release files changed: compare the ETags in `reports/manifests/data_manifest.tsv` and read `reports/qc/sample_count_report.md`. The related-ID decoder refuses unfamiliar Hail formats; if so, read the table with Hail instead. Stage 3 does not depend on scope (sample metadata is always fetched in full, about 7 MB) |
| 4 | 19,000+ MANE Select transcripts, one per gene, with CDS lengths that are multiples of 3; every scope gene found; each reference contig has the GRCh38 length and only ACGTN; every AlphaMissense and ClinVar REF allele equals the GRCh38 base; ClinVar comes from a dated release that did not change during the read; disk budget | `blocked` at `dms_genes` scope if `tables/dms_candidate_inventory.tsv` is missing: run Stage 5 first. A REF mismatch means a build or contig-naming problem: stop and investigate, do not filter it away. About 3 minutes at smoke scope, most of it streaming AlphaMissense (643 MB) |
| 5 | MaveDB snapshot taken or reused; every published human score set classified (candidate or excluded, with a reason); every scope gene has at least one candidate assay; disk budget. The stage **passes** before review, and its status note shows how many in-scope assays are reviewed | A MaveDB outage gives HTTP errors after 5 retries, so run it again later. To pick up new MaveDB records, run `ancal stage5 --refresh` (this changes the inventory, so note it in `docs/DECISIONS.md`). About 1 minute for a fresh snapshot |
