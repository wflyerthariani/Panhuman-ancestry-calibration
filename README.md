# PanHuman: AVI ancestry-calibration audit

Pipeline for **Outcome 1** of the PanHuman scientific outputs: does the AlphaGenome Atlas
**Aggregate Variant Impact (AVI)** score perform consistently across the 80 populations
of the gnomAD-harmonized HGDP + 1000 Genomes callset?

The analysis protocol is Caroline Duncan's "Ancestry calibration audit" canvas
(#4-scientific-outputs). This repository currently builds **deliverable 1, the
data-readiness package**. Deliverable 1 calculates **no performance metrics**
(no AVI–DMS correlations, CADD comparisons or ClinVar AUROCs). Those wait until
after the Step 40 review and config freeze.

## Rules

- **Never edit files by hand** in `data_store/`, `tables/`, `reports/` or `deliverables/`.
  Every file there is produced by a script in `src/ancal/` run through `workflow/Snakefile`.
- **Every download** is recorded in `data_store/manifests/data_manifest.tsv` (protocol Step 1).
- **Local disk is capped** by `storage.local_budget_gb` in `config/data_sources.yaml`
  (1 GB for now). Stages refuse downloads that would exceed it. See [SCALING_UP.md](SCALING_UP.md).

## Setup

```
make env          # create/update the `ancal` conda env (bcftools, samtools, snakemake, polars, ...)
make stage0       # Stage 0 gate: tools, imports, unit tests, disk budget
make status       # one line per stage: pass / fail / blocked
```

## Stages

Run the stages one at a time. Each one writes `reports/status/stageN.json` and stops the
pipeline if its gate fails. A `blocked` stage (missing access or too little disk) is not
a failure. It is listed automatically in the unresolved-issues report (item 15).

| Stage | Make target | Deliverable items | Protocol steps | Status |
|---|---|---|---|---|
| 0 Environment and scaffold | `make stage0` | 1 | 1 | implemented |
| 1 Analysis configuration | `make stage1` | 2 | 2 | planned |
| 2 Provenance and manifests | `make stage2` | 3 | 1 | planned |
| 3 Samples and populations | `make stage3` | 4, 5, 6 | 3–5 | planned |
| 4 Reference and annotation (scope-aware) | `make stage4` | (inputs) | 6, 29 | planned |
| 5 DMS candidate inventory | `make stage5` | 9 | 12 | planned — **review stop** |
| 6 DMS harmonization | `make stage6` | 11 | 13 | planned |
| 7 Population variants | `make stage7` | 7 | 6–8 | planned |
| 8 AVI join and coverage | `make stage8` | 8 | 9–11 | planned (needs AVI access) |
| 9 Proposed DMS gene universe | `make stage9` | 10 | 12, 15 | planned — **review stop** |
| 10 ClinVar eligibility | `make stage10` | 12 | 29 | planned |
| 11 Attrition, forecast, issues, package | `make stage11` | 13, 14, 15 | 14, 18 | planned |

## Layout

```
config/                 analysis.yaml, data_sources.yaml (scope, storage, URIs), frozen/
src/ancal/              pipeline package (settings, budget, status, provenance, sources/, stages/)
workflow/Snakefile      one rule per stage
tests/                  unit tests
resources/regions/      region BEDs for the current scope (generated)
data_store/             all data (git-ignored): manifests, reference, population_metadata,
                        variants, avi, dms, clinvar, intermediate, analysis_datasets
tables/ figures/        generated outputs
reports/                qc/, exclusions/, status/, archive/
deliverables/d1/        the data-readiness package index (one row per item 1–15)
logs/                   per-stage logs
```

## Scope

`scope` in `config/data_sources.yaml` controls how much data is fetched:

- `smoke` (current): 2–3 DMS-rich genes on chr17. Proves every stage end-to-end within 1 GB.
- `dms_genes`: all candidate DMS genes. Needed for the final gene universe and item counts.
- `exome`: all coding regions, for the full ClinVar benchmark.

Results produced at `smoke` scope are marked as such in `deliverables/d1/README.md`.
