# PanHuman: AVI ancestry-calibration audit

Pipeline for **Outcome 1** of the PanHuman scientific outputs: does the AlphaGenome Atlas
**Aggregate Variant Impact (AVI)** score perform consistently across the 80 populations
of the gnomAD-harmonized HGDP + 1000 Genomes callset?

The analysis protocol is Caroline Duncan's "Ancestry calibration audit" canvas
(#4-scientific-outputs). This repository builds **deliverable 1, the data-readiness package**
(protocol Steps 1–16, 18, 29, 35 and 39, up to the Step 40 review). The performance analysis of
the second deliverable (Steps 17 and 19–38) is not implemented here. Deliverable 1 calculates **no performance metrics**
(no AVI–DMS correlations, CADD comparisons or ClinVar AUROCs). Those wait until
after the Step 40 review and config freeze.

## Decisions and open review items

**[docs/DECISIONS.md](docs/DECISIONS.md)** records every important decision, discovery and deferred
question, with its evidence. It also holds the open review checklist for the team. It is
updated in every stage. Start there before reading any output.

## Rules

- **Never edit files by hand** in `data_store/`, `tables/`, `reports/` or `deliverables/`.
  Every file there is produced by a script in `src/ancal/` run through `workflow/Snakefile`.
- **Every download** is recorded in `reports/manifests/data_manifest.tsv` (protocol Step 1), with source
  URL, ETag, size and SHA-256, the output checksum, the command and the code commit.
  Large sources are hashed as they stream past but never stored in full.
- **Local disk is capped** by `storage.local_budget_gb` in `config/data_sources.yaml`
  (1 GB for now). Stages refuse downloads that would exceed it. See [SCALING_UP.md](SCALING_UP.md).

## Setup

```sh
make env          # create/update the `ancal` conda env (bcftools, samtools, snakemake, polars, ...)
make stage0       # Stage 0 gate: tools, imports, unit tests, disk budget
make status       # one line per stage: pass / fail / blocked
```

## Stages

`make d1` runs every stage that is out of date, in order, and prints the status of each.
`make stageN` runs up to one stage. Each stage writes `reports/status/stageN.json` and stops the
pipeline if its gate fails. A `blocked` stage (missing access, a pending review, too little disk)
is not a failure: later stages still run and report what they can, and the block is listed
automatically in the unresolved-issues report (item 15).

Two settings in `config/data_sources.yaml` decide how far an end-to-end run gets:

- **`review_gates.dms_review`**: `provisional` lets Stage 6 run on the pipeline's proposed score
  directions for assays nobody has reviewed yet, and marks everything downstream PROVISIONAL.
  `strict` blocks Stage 6 until the review in `config/dms_curation.tsv` is complete, and is
  required before the config can be frozen.
- **`sources.avi`**: while `uri` is null, Stage 8 is blocked and Stages 9–11 report AVI-dependent
  counts as not available (evaluable-variant counts become upper bounds).

| Stage | Make target | Deliverable items | Protocol steps | Status |
|---|---|---|---|---|
| 0 Environment and scaffold | `make stage0` | 1 | 1 | implemented |
| 1 Analysis configuration | `make stage1` | 2 | 2 | implemented — open items in `reports/qc/config_status.md` |
| 2 Provenance and manifests | `make stage2` | 3 | 1 | implemented — `reports/manifests/` |
| 3 Samples and populations | `make stage3` | 4, 5, 6 | 3–5 | implemented — **review stop**: `reports/qc/sample_count_report.md` |
| 4 Reference and annotation (scope-aware) | `make stage4` | (inputs) | 6, 29 | implemented — `reports/qc/reference_report.md` |
| 5 DMS candidate inventory | `make stage5` | 9 | 12 | implemented — **review stop**: biologist guide `reports/review/dms_direction_review.md` → `config/dms_curation.tsv` |
| 6 DMS harmonization | `make stage6` | 11 | 13 | implemented — `reports/qc/dms_harmonization_report.md` |
| 7 Population variants | `make stage7` | 7 | 6–8 | implemented — `reports/qc/variant_normalization_report.md` |
| 8 AVI join and coverage | `make stage8` | 8 | 9–11 | implemented — **blocked** until the `avi:` source is set in `config/data_sources.yaml` |
| 9 Proposed DMS gene universe | `make stage9` | 10 | 12, 35 | implemented — **review stop**: `reports/qc/gene_universe_report.md` |
| 10 ClinVar eligibility | `make stage10` | 12 | 29 | implemented — `reports/qc/clinvar_eligibility_report.md` |
| 11 Attrition, forecast, issues, package | `make stage11` | 13, 14, 15 | 14, 15, 18, 39 | implemented — `deliverables/d1/README.md` |

## Layout

```sh
config/                 analysis.yaml, data_sources.yaml (scope, storage, URIs), frozen/
src/ancal/              pipeline package (settings, config, budget, status, provenance, consequence,
                        sources/, stages/)
workflow/Snakefile      one rule per stage
tests/                  unit tests
resources/regions/      region BEDs for the current scope (generated)
data_store/             all data (git-ignored): reference, population_metadata,
                        variants, avi, dms, clinvar, intermediate, analysis_datasets
tables/ figures/        generated outputs
reports/                manifests/ (tracked provenance), qc/, exclusions/, status/, archive/
deliverables/d1/        the data-readiness package index (one row per item 1–15)
logs/                   per-stage logs
```

## Scope

`scope` in `config/data_sources.yaml` controls how much data is fetched:

- `smoke` (current): 2–3 DMS-rich genes on chr17. Proves every stage end-to-end within 1 GB.
- `dms_genes`: all candidate DMS genes. Needed for the final gene universe and item counts.
- `exome`: all coding regions, for the full ClinVar benchmark.

Results produced at `smoke` scope are marked as such in `deliverables/d1/README.md`.
