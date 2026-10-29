SHELL := /usr/bin/env bash

CONFIG ?= conf/demo.yaml
SAMPLES ?=

.PHONY: help lint validate smoke dry-run run clean tree

help:
	@printf "%s\n" \
	"Available targets:" \
	"  make lint   - shell and python syntax checks" \
	"  make validate - check sample sheet and demo input paths" \
	"  make smoke  - run minimal smoke test" \
	"  make dry-run - validate entry arguments without creating outputs" \
	"  make run    - execute workflow entry with config and sample sheet" \
	"  make clean  - remove generated smoke outputs" \
	"  make tree   - print current template layout"

lint:
	bash -n workflow/run.sh
	bash -n tests/run_smoke.sh
	python3 -m py_compile scripts/*.py
	python3 -m py_compile scripts/00_shared/*.py

validate:
	samples_path="${SAMPLES}"; \
	if [[ -z "$$samples_path" ]]; then \
		samples_path="$$(python3 scripts/read_config_value.py --config "$(CONFIG)" --key samples_tsv --required)"; \
	fi; \
	python3 scripts/check_samples.py --samples "$$samples_path" --config "$(CONFIG)"

smoke:
	bash tests/run_smoke.sh

dry-run:
	if [[ -n "$(SAMPLES)" ]]; then \
		bash workflow/run.sh --config "$(CONFIG)" --samples "$(SAMPLES)" --dry-run; \
	else \
		bash workflow/run.sh --config "$(CONFIG)" --dry-run; \
	fi

run:
	if [[ -n "$(SAMPLES)" ]]; then \
		bash workflow/run.sh --config "$(CONFIG)" --samples "$(SAMPLES)"; \
	else \
		bash workflow/run.sh --config "$(CONFIG)"; \
	fi

clean:
	rm -rf tmp/smoke-runtime

tree:
	find . -maxdepth 3 | sort
