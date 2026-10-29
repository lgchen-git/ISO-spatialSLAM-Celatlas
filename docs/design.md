# Template Design Notes

This file only records the extra design rules that are easy to lose during project reuse.

## Keep Stable

- Keep one top-level entry at `workflow/run.sh`.
- Adapt projects through `conf/` before changing workflow code.
- Keep runtime files under the configured result, log, and scratch roots.
- Prefer sample-first formal results: `{results_root}/{sample_id}/{stage}`.
- Keep large recomputable intermediates under `tmp_root`, not formal results.

## Keep Reusable

- Parameterize sample names, paths, thresholds, and resources.
- Default to multi-sample execution with explicit resource limits.
- Prefer small local edits over layout changes.
- Keep Codex reuse in repo skills and root instructions before adding plugin packaging.

## Keep Deliverable

- Write key metrics as tables before prose summaries.
- Keep plot styling consistent and publication-ready.
- Keep editable report sources separate from final deliverables.
