#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root_dir"

test_cache_root="$root_dir/tmp/test-cache"
export TMPDIR="$test_cache_root/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
export XDG_CACHE_HOME="$test_cache_root/xdg"
export MPLCONFIGDIR="$test_cache_root/matplotlib"
export NUMBA_CACHE_DIR="$test_cache_root/numba"
export PYTHONPYCACHEPREFIX="$test_cache_root/pycache"
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$PYTHONPYCACHEPREFIX"

printf 'Running smoke test in %s\n' "$root_dir"

mkdir -p data/example
smoke_inputs=(
    data/example/demo_A_R1.fastq.gz
    data/example/demo_A_R2.fastq.gz
    data/example/demo_B_R1.fastq.gz
    data/example/demo_B_R2.fastq.gz
)
touch "${smoke_inputs[@]}"
mkdir -p tmp/smoke-fixtures
touch tmp/smoke-fixtures/epi-logo.png
smoke_root="tmp/smoke-runtime"
smoke_results="${smoke_root}/results/demo"
smoke_logs="${smoke_root}/logs/demo"
smoke_tmp="${smoke_root}/tmp/demo"
cleanup_smoke_inputs() {
    rm -f "${smoke_inputs[@]}"
    rm -f tmp/smoke-fixtures/epi-logo.png
    rm -rf "$smoke_root"
    rm -rf "$test_cache_root"
}
trap cleanup_smoke_inputs EXIT

bash -n workflow/run.sh
bash -n tests/run_smoke.sh
python3 -m py_compile scripts/*.py
python3 -m py_compile scripts/00_shared/*.py
python3 -m py_compile tests/test_lr_reusable.py

python3 scripts/check_samples.py --samples "$(python3 scripts/read_config_value.py --config conf/demo.yaml --key samples_tsv --required)"
bash workflow/run.sh --config conf/demo.yaml --dry-run
bash workflow/run.sh --config conf/demo.yaml
python3 scripts/generate_report_skeleton.py \
    --project-name demo_pipeline \
    --analysis-type template \
    --run-id demo \
    --report-title "Demo Report" \
    --report-style clean \
    --report-logo tmp/smoke-fixtures/epi-logo.png \
    --plot-palette okabe_ito \
    --plot-font-family Arial \
    --plot-base-size 10 \
    --plot-dpi 300 \
    --sample-count 2 \
    --config-path conf/demo.yaml \
    --samples-path conf/samples.tsv \
    --summary-path "${smoke_results}/tables/run_summary.tsv" \
    --figures-dir "${smoke_results}/figures" \
    --tables-dir "${smoke_results}/tables" \
    --reports-dir "${smoke_results}/reports" \
    --deliverables-dir "${smoke_results}/deliverables" \
    --editable-dir "${smoke_results}/editable" >/dev/null

test -f "${smoke_results}/run_manifest.tsv"
test -f "${smoke_results}/tables/run_summary.tsv"
test -f "${smoke_logs}/pipeline/run.log"
test -f "${smoke_tmp}/checkpoints/run.state"
test -d "${smoke_tmp}/work"
test -d "${smoke_results}/deliverables"
test -d "${smoke_results}/editable"
test -d "${smoke_results}/figures"
test -d "${smoke_results}/reports"
test -d "${smoke_results}/qc"
test -f "${smoke_results}/reports/demo_report.md"
test -f "${smoke_results}/deliverables/demo_report.md"
test -f "${smoke_results}/editable/demo_report.md"
test -f "${smoke_results}/editable/report_style.tsv"
test -f "${smoke_results}/editable/report_sections.tsv"
test -f "${smoke_results}/editable/epi-logo.png"

grep -q '^run_id[[:space:]]demo$' "${smoke_results}/run_manifest.tsv"
grep -q '^project_name[[:space:]]demo_pipeline$' "${smoke_results}/run_manifest.tsv"
grep -q '^execution_profile[[:space:]]local$' "${smoke_results}/run_manifest.tsv"
grep -q '^report_logo[[:space:]]epi-logo.png$' "${smoke_results}/run_manifest.tsv"
grep -q '^metric[[:space:]]value$' "${smoke_results}/tables/run_summary.tsv"
grep -q '^# Demo Report$' "${smoke_results}/reports/demo_report.md"
grep -q '^| Palette | okabe_ito |$' "${smoke_results}/reports/demo_report.md"
grep -q 'template run completed' "${smoke_logs}/pipeline/run.log"

python3 tests/test_lr_reusable.py

printf 'Eight-stage reusable workflow smoke test passed.\n'
