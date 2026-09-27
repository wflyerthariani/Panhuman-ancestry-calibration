# Run every target inside the `ancal` conda env.
ENV      ?= ancal
RUN      := conda run --no-capture-output -n $(ENV)
SCOPE    ?=
SMK_ARGS := --snakefile workflow/Snakefile --cores 1 --keep-incomplete $(if $(SCOPE),--config scope=$(SCOPE),)

.PHONY: env check-env test budget status stage0 stage1 stage2 stage3 freeze-config

## Create or update the conda environment and install the ancal package into it.
env:
	@if conda env list | grep -qE '^$(ENV)[[:space:]]'; then \
		mamba env update -n $(ENV) -f environment.yml --prune; \
	else \
		mamba env create -n $(ENV) -f environment.yml; \
	fi
	$(RUN) python -m pip install -e . --no-deps

test:
	$(RUN) pytest -q

## Stage 0 gate: tools, imports, unit tests, disk budget -> reports/status/stage0.json
check-env:
	@if $(RUN) pytest -q; then ok=1; else ok=0; fi; \
	$(RUN) ancal check-env --pytest-ok $$ok

budget:
	$(RUN) ancal budget-check $(if $(SCOPE),--scope $(SCOPE),)

status:
	$(RUN) ancal status

stage0:
	$(RUN) snakemake $(SMK_ARGS) stage0

stage1:
	$(RUN) snakemake $(SMK_ARGS) stage1

## Only succeeds once every deferred item in config/analysis.yaml is confirmed.
freeze-config:
	$(RUN) ancal freeze-config

stage2:
	$(RUN) snakemake $(SMK_ARGS) stage2

stage3:
	$(RUN) snakemake $(SMK_ARGS) stage3
