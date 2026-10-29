# LR-Celatlas Spatial SLAM Pipeline

This repository contains a reusable, configuration-driven workflow for
LR-Celatlas spatial SLAM data.

The current shape follows the stable conventions that still show up across official Snakemake and nf-core guidance in 2025:

- config lives in `conf/`
- workflow code stays under `workflow/`
- helper code stays under `scripts/`
- runtime artifacts follow the configured `results_root`, `logs_root`, and `tmp_root`
- one top-level entry point remains at `workflow/run.sh`

That keeps the template generic, auditable, and cheap to adapt for future projects without dragging in a full workflow framework too early.

## Default Design Priorities

- reproducible runs with config-driven parameters
- batch-safe multi-sample processing
- resumable execution through stable run directories
- clear separation of logs, results, and scratch space
- delivery-ready outputs plus editable source files
- compact QC summaries and table-first key metrics
- FLAMES/gffcompare full-length transcript structure analysis
- whole-field long-read transcriptome clustering with manually entered parameters
- SLAM overlays on the same whole-transcriptome spatial clusters
- bin50 spatial QC tables with UMI, long-read molecule, detected-gene and sample-level isoform metrics

## Layout

```text
conf/       project config, templates, sample sheets
data/       tiny example inputs only
docs/       design notes and changelog
envs/       reusable environment definitions
logs/       runtime logs only
results/    formal outputs only
scripts/    helper utilities
tests/      smoke tests
tmp/        scratch files only
workflow/   unified entry point and future engine wrappers
```

## Recommended Reuse Pattern

1. Copy `conf/lr_celatlas_spatial_slam.template.yaml` to `conf/<project>.yaml`.
2. Copy `conf/samples.lr_celatlas_spatial_slam.template.tsv` to `conf/<project>.samples.tsv`.
3. Replace every `CHANGE_ME` value and enter the sample paths.
4. Run `make validate CONFIG=conf/<project>.yaml`.
5. Run `make smoke`.
6. Extend `workflow/` or `scripts/` only after config is no longer enough.

## Config Contract

The template keeps the config intentionally small and top-level so it can be read from Bash without extra dependencies.

Key fields:

- `project_name`
- `analysis_type`
- `run_id`
- `samples_tsv`
- `input_root`
- `results_root`
- `logs_root`
- `tmp_root`
- `execution_profile`
- `resume_mode`
- `max_parallel_samples`
- `threads`
- `memory_gb`

Optional feature toggles are also kept as top-level booleans for future modules.

Runtime subdirectories are also configurable through top-level keys such as `tables_dirname`, `figures_dirname`, `reports_dirname`, and `deliverables_dirname`.

For a real dataset, absolute roots may point outside the Git worktree. Runtime
paths are declared in the project config and remain independent of repository
location.

Plot and report defaults are configured through:

- `plot_palette`
- `plot_font_family`
- `plot_base_size`
- `plot_dpi`
- `report_title`
- `report_logo`
- `report_style`

When `enable_report: true`, the template can generate a formal report skeleton with:

- `results/<run_id>/reports/<run_id>_report.md`
- `results/<run_id>/deliverables/<run_id>_report.md`
- `results/<run_id>/editable/<run_id>_report.md`
- `results/<run_id>/editable/report_style.tsv`
- `results/<run_id>/editable/report_sections.tsv`

Use `conf/project.template.yaml` as the single source template for new project configs.

## Commands

```bash
make lint
make validate
make dry-run
make run
make smoke
make clean
```

The reusable LR-Celatlas spatial SLAM workflow is documented in
[`docs/lr-celatlas-spatial-slam.md`](docs/lr-celatlas-spatial-slam.md).

The formal sample workflow ends at `08_report`. The separate PPT and
exploratory isoform package is generated only when requested with
`--promo-figures`; it is not part of the normal stage completion sequence.

`make validate` and `make run` now default to the sample sheet declared inside `CONFIG`, which keeps the project config as the main control surface.

## Notes

- `conf/demo.yaml` is the smoke-test config.
- `conf/project.template.yaml` is the copy-edit starting point for new projects.
- `workflow/run.sh` stays deliberately small; add framework-specific wrappers under `workflow/` instead of replacing the unified entry pattern.
- The default report logo path is `epi-logo.png`, so downstream report steps can include the shared brand asset without project-specific hard coding.
- `scripts/build_lr_spatial_slam_report.py` creates the self-contained Chinese HTML report.
- Final reports are named `<sample_id>_report.html` under the sample `08_report` stage.
- The repository is released under the MIT License; runtime data and machine-specific configuration remain excluded.
