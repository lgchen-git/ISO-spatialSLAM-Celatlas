#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash workflow/run.sh --config conf/demo.yaml
  bash workflow/run.sh --config conf/demo.yaml --samples conf/samples.tsv
  bash workflow/run.sh --config conf/run.yaml --from 03_align --to 08_report
  bash workflow/run.sh --config conf/run.yaml --from 08_report --to 08_report --report-title "长读长空间转录组与 SLAM 分析报告"
  bash workflow/run.sh --config conf/run.yaml --refresh-visuals
  bash workflow/run.sh --config conf/run.yaml --promo-figures

Purpose:
  Unified template entry point for bioinformatics workflow development.

Notes:
  This script is intentionally minimal.
  Real projects should extend behavior through config and helper scripts.
  --samples overrides samples_tsv from the config file when needed.
  --force-stage may be repeated for explicitly regenerated LR stages.
  --promo-figures consumes completed sample results and writes a PPT figure pack.
EOF
}

CONFIG=""
SAMPLES=""
DRY_RUN="false"
FROM_STAGE="01_read_qc"
TO_STAGE="08_report"
FORCE_STAGES=()
REPORT_TITLE_OVERRIDE=""
REFRESH_VISUALS="false"
PROMO_FIGURES="false"
STAGE_SELECTOR_SET="false"
CONFIG_READER="scripts/read_config_value.py"

read_config_value() {
    local key="$1"
    python3 "$CONFIG_READER" --config "$CONFIG" --key "$key"
}

default_if_empty() {
    local value="$1"
    local fallback="$2"

    if [[ -n "$value" ]]; then
        printf '%s\n' "$value"
    else
        printf '%s\n' "$fallback"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG="${2:-}"
            shift 2
            ;;
        --samples)
            SAMPLES="${2:-}"
            shift 2
            ;;
        --dry-run)
            DRY_RUN="true"
            shift
            ;;
        --from)
            FROM_STAGE="${2:-}"
            STAGE_SELECTOR_SET="true"
            shift 2
            ;;
        --to)
            TO_STAGE="${2:-}"
            STAGE_SELECTOR_SET="true"
            shift 2
            ;;
        --force-stage)
            FORCE_STAGES+=("${2:-}")
            STAGE_SELECTOR_SET="true"
            shift 2
            ;;
        --refresh-visuals)
            REFRESH_VISUALS="true"
            shift
            ;;
        --promo-figures)
            PROMO_FIGURES="true"
            shift
            ;;
        --report-title)
            REPORT_TITLE_OVERRIDE="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'Unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ "$REFRESH_VISUALS" == "true" && "$STAGE_SELECTOR_SET" == "true" ]]; then
    printf '%s\n' '--refresh-visuals cannot be combined with --from, --to, or --force-stage.' >&2
    exit 1
fi

if [[ -z "$CONFIG" ]]; then
    printf '--config is required.\n' >&2
    usage >&2
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    printf 'Config not found: %s\n' "$CONFIG" >&2
    exit 1
fi

if [[ ! -f "$CONFIG_READER" ]]; then
    printf 'Config reader not found: %s\n' "$CONFIG_READER" >&2
    exit 1
fi

if [[ -z "$SAMPLES" ]]; then
    SAMPLES="$(read_config_value samples_tsv)"
fi

if [[ -z "$SAMPLES" ]]; then
    printf 'Sample sheet path is missing. Set --samples or samples_tsv in %s\n' "$CONFIG" >&2
    exit 1
fi

if [[ ! -f "$SAMPLES" ]]; then
    printf 'Sample sheet not found: %s\n' "$SAMPLES" >&2
    exit 1
fi

run_id="$(read_config_value run_id)"
project_name="$(read_config_value project_name)"
data_root="$(read_config_value data_root)"
results_root="$(read_config_value results_root)"
logs_root="$(read_config_value logs_root)"
tmp_root="$(read_config_value tmp_root)"
sample_first_layout="$(read_config_value sample_first_layout)"
git_dirty_policy="$(read_config_value git_dirty_policy)"
execution_profile="$(read_config_value execution_profile)"
resume_mode="$(read_config_value resume_mode)"
sample_parallelism="$(read_config_value sample_parallelism)"
max_parallel_samples="$(read_config_value max_parallel_samples)"
threads="$(read_config_value threads)"
memory_gb="$(read_config_value memory_gb)"
conda_env_path="$(read_config_value conda_env_path)"
tables_dirname="$(read_config_value tables_dirname)"
figures_dirname="$(read_config_value figures_dirname)"
reports_dirname="$(read_config_value reports_dirname)"
deliverables_dirname="$(read_config_value deliverables_dirname)"
editable_dirname="$(read_config_value editable_dirname)"
qc_dirname="$(read_config_value qc_dirname)"
pipeline_log_dirname="$(read_config_value pipeline_log_dirname)"
checkpoints_dirname="$(read_config_value checkpoints_dirname)"
work_dirname="$(read_config_value work_dirname)"
plot_palette="$(read_config_value plot_palette)"
plot_font_family="$(read_config_value plot_font_family)"
plot_base_size="$(read_config_value plot_base_size)"
plot_dpi="$(read_config_value plot_dpi)"
report_title="$(read_config_value report_title)"
report_logo="$(read_config_value report_logo)"
report_style="$(read_config_value report_style)"
analysis_type="$(read_config_value analysis_type)"
enable_report="$(read_config_value enable_report)"
enable_lr_spatial_slam="$(read_config_value enable_lr_spatial_slam)"

if [[ -z "$run_id" ]]; then
    run_id="run"
fi

if [[ -z "$project_name" ]]; then
    project_name="bio_pipeline_project"
fi

if [[ -z "$results_root" ]]; then
    results_root="results"
fi

if [[ -z "$logs_root" ]]; then
    logs_root="logs"
fi

if [[ -z "$tmp_root" ]]; then
    tmp_root="tmp"
fi

data_root="$(default_if_empty "$data_root" "$PWD")"
sample_first_layout="$(default_if_empty "$sample_first_layout" "false")"
git_dirty_policy="$(default_if_empty "$git_dirty_policy" "warn")"
execution_profile="$(default_if_empty "$execution_profile" "local")"
resume_mode="$(default_if_empty "$resume_mode" "resume")"
sample_parallelism="$(default_if_empty "$sample_parallelism" "true")"
max_parallel_samples="$(default_if_empty "$max_parallel_samples" "1")"
threads="$(default_if_empty "$threads" "1")"
memory_gb="$(default_if_empty "$memory_gb" "4")"
conda_env_path="$(default_if_empty "$conda_env_path" "")"
tables_dirname="$(default_if_empty "$tables_dirname" "tables")"
figures_dirname="$(default_if_empty "$figures_dirname" "figures")"
reports_dirname="$(default_if_empty "$reports_dirname" "reports")"
deliverables_dirname="$(default_if_empty "$deliverables_dirname" "deliverables")"
editable_dirname="$(default_if_empty "$editable_dirname" "editable")"
qc_dirname="$(default_if_empty "$qc_dirname" "qc")"
pipeline_log_dirname="$(default_if_empty "$pipeline_log_dirname" "pipeline")"
checkpoints_dirname="$(default_if_empty "$checkpoints_dirname" "checkpoints")"
work_dirname="$(default_if_empty "$work_dirname" "work")"
plot_palette="$(default_if_empty "$plot_palette" "okabe_ito")"
plot_font_family="$(default_if_empty "$plot_font_family" "DejaVu Sans")"
plot_base_size="$(default_if_empty "$plot_base_size" "10.5")"
plot_dpi="$(default_if_empty "$plot_dpi" "300")"
report_title="$(default_if_empty "$report_title" "$project_name report")"
if [[ -n "$REPORT_TITLE_OVERRIDE" ]]; then
    report_title="$REPORT_TITLE_OVERRIDE"
fi
report_logo="$(default_if_empty "$report_logo" "epi-logo.png")"
report_style="$(default_if_empty "$report_style" "clean")"
analysis_type="$(default_if_empty "$analysis_type" "generic")"
enable_report="$(default_if_empty "$enable_report" "false")"
enable_lr_spatial_slam="$(default_if_empty "$enable_lr_spatial_slam" "false")"

if [[ -n "$conda_env_path" && -d "$conda_env_path/bin" ]]; then
    export PATH="${conda_env_path}/bin:${PATH}"
fi

python3 scripts/check_samples.py --samples "$SAMPLES" --config "$CONFIG" >/dev/null

if [[ "$sample_first_layout" == "true" ]]; then
    results_dir="$results_root"
    project_dir="${results_root}/_project"
    logs_dir="${project_dir}/03_logs"
    tmp_dir="$tmp_root"
    tables_dir="${project_dir}/01_tables"
    figures_dir="${project_dir}/02_figures"
    reports_dir="${project_dir}/04_tool_output/reports"
    deliverables_dir="${project_dir}/04_tool_output/deliverables"
    editable_dir="${tmp_dir}/editable"
    qc_dir="${project_dir}/01_tables/qc"
    pipeline_logs_dir="${project_dir}/03_logs"
    checkpoints_dir="${tmp_dir}/checkpoints"
    work_dir="${tmp_dir}/work"
    manifest_path="${project_dir}/05_metadata/run_manifest.tsv"
    summary_path="${project_dir}/01_tables/run_summary.tsv"
    checkpoint_path="${checkpoints_dir}/run.state"
    log_path="${project_dir}/03_logs/run.log"
else
    results_dir="${results_root}/${run_id}"
    logs_dir="${logs_root}/${run_id}"
    tmp_dir="${tmp_root}/${run_id}"
    tables_dir="${results_dir}/${tables_dirname}"
    figures_dir="${results_dir}/${figures_dirname}"
    reports_dir="${results_dir}/${reports_dirname}"
    deliverables_dir="${results_dir}/${deliverables_dirname}"
    editable_dir="${results_dir}/${editable_dirname}"
    qc_dir="${results_dir}/${qc_dirname}"
    pipeline_logs_dir="${logs_dir}/${pipeline_log_dirname}"
    checkpoints_dir="${tmp_dir}/${checkpoints_dirname}"
    work_dir="${tmp_dir}/${work_dirname}"
    manifest_path="${results_dir}/run_manifest.tsv"
    summary_path="${tables_dir}/run_summary.tsv"
    checkpoint_path="${checkpoints_dir}/run.state"
    log_path="${pipeline_logs_dir}/run.log"
fi

if [[ "$DRY_RUN" == "true" ]]; then
    printf 'Dry run OK.\n'
    printf 'project_name: %s\n' "$project_name"
    printf 'run_id: %s\n' "$run_id"
    printf 'config: %s\n' "$CONFIG"
    printf 'samples: %s\n' "$SAMPLES"
    printf 'data_root: %s\n' "$data_root"
    printf 'sample_first_layout: %s\n' "$sample_first_layout"
    printf 'execution_profile: %s\n' "$execution_profile"
    printf 'resume_mode: %s\n' "$resume_mode"
    printf 'sample_parallelism: %s\n' "$sample_parallelism"
    printf 'max_parallel_samples: %s\n' "$max_parallel_samples"
    printf 'threads: %s\n' "$threads"
    printf 'memory_gb: %s\n' "$memory_gb"
    printf 'report_title: %s\n' "$report_title"
    printf 'conda_env_path: %s\n' "$conda_env_path"
    printf 'results_dir: %s\n' "$results_dir"
    printf 'logs_dir: %s\n' "$logs_dir"
    printf 'tmp_dir: %s\n' "$tmp_dir"
    printf 'enable_lr_spatial_slam: %s\n' "$enable_lr_spatial_slam"
    exit 0
fi

mkdir -p \
    "$results_dir" \
    "$logs_dir" \
    "$tmp_dir" \
    "$tables_dir" \
    "$figures_dir" \
    "$reports_dir" \
    "$deliverables_dir" \
    "$editable_dir" \
    "$qc_dir" \
    "$pipeline_logs_dir" \
    "$checkpoints_dir" \
    "$work_dir" \
    "$(dirname "$manifest_path")"

export TMPDIR="$work_dir"
export TMP="$work_dir"
export TEMP="$work_dir"
export XDG_CACHE_HOME="$tmp_dir/xdg-cache"
export MPLCONFIGDIR="$tmp_dir/mplconfig"
export MPLBACKEND="Agg"
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
if [[ -d "$PWD/.git-local" && ! -f "$PWD/.git/HEAD" ]]; then
    export GIT_DIR="$PWD/.git-local"
    export GIT_WORK_TREE="$PWD"
fi

provenance_path="$(dirname "$manifest_path")/git_provenance.tsv"
python3 scripts/capture_git_provenance.py \
    --repo-root "$PWD" \
    --output "$provenance_path" \
    --dirty-policy "$git_dirty_policy"
git_commit="$(awk -F '\t' '$1 == "commit" {print $2}' "$provenance_path")"
git_branch="$(awk -F '\t' '$1 == "branch" {print $2}' "$provenance_path")"
git_dirty="$(awk -F '\t' '$1 == "dirty" {print $2}' "$provenance_path")"
pipeline_release_tag="$(awk -F '\t' '$1 == "pipeline_release_tag" {print $2}' "$provenance_path")"

timestamp="$(date '+%Y-%m-%d %H:%M:%S')"
sample_count="$(( $(wc -l < "$SAMPLES") - 1 ))"

if [[ "$enable_lr_spatial_slam" == "true" || "$analysis_type" == "lr_celatlas_spatial_slam" ]]; then
    lr_args=(
        --config "$CONFIG" \
        --samples "$SAMPLES" \
        --results-dir "$results_dir" \
        --logs-dir "$pipeline_logs_dir" \
        --tmp-dir "$work_dir" \
        --threads "$threads" \
        --from "$FROM_STAGE" \
        --to "$TO_STAGE"
    )
    for stage in "${FORCE_STAGES[@]}"; do
        lr_args+=(--force-stage "$stage")
    done
    if [[ -n "$REPORT_TITLE_OVERRIDE" ]]; then
        lr_args+=(--report-title "$REPORT_TITLE_OVERRIDE")
    fi
    if [[ "$REFRESH_VISUALS" == "true" ]]; then
        lr_args+=(--refresh-visuals)
    fi
    if [[ "$PROMO_FIGURES" == "true" ]]; then
        lr_args+=(--promo-figures)
    fi
    python3 scripts/run_lr_spatial_slam.py "${lr_args[@]}"
fi

{
    printf 'timestamp\t%s\n' "$timestamp"
    printf 'project_name\t%s\n' "$project_name"
    printf 'config\t%s\n' "$CONFIG"
    printf 'samples\t%s\n' "$SAMPLES"
    printf 'data_root\t%s\n' "$data_root"
    printf 'results_root\t%s\n' "$results_root"
    printf 'tmp_root\t%s\n' "$tmp_root"
    printf 'sample_first_layout\t%s\n' "$sample_first_layout"
    printf 'run_id\t%s\n' "$run_id"
    printf 'execution_profile\t%s\n' "$execution_profile"
    printf 'resume_mode\t%s\n' "$resume_mode"
    printf 'sample_parallelism\t%s\n' "$sample_parallelism"
    printf 'max_parallel_samples\t%s\n' "$max_parallel_samples"
    printf 'threads\t%s\n' "$threads"
    printf 'memory_gb\t%s\n' "$memory_gb"
    printf 'conda_env_path\t%s\n' "$conda_env_path"
    printf 'tables_dir\t%s\n' "$tables_dir"
    printf 'figures_dir\t%s\n' "$figures_dir"
    printf 'reports_dir\t%s\n' "$reports_dir"
    printf 'deliverables_dir\t%s\n' "$deliverables_dir"
    printf 'editable_dir\t%s\n' "$editable_dir"
    printf 'qc_dir\t%s\n' "$qc_dir"
    printf 'pipeline_logs_dir\t%s\n' "$pipeline_logs_dir"
    printf 'checkpoints_dir\t%s\n' "$checkpoints_dir"
    printf 'work_dir\t%s\n' "$work_dir"
    printf 'plot_palette\t%s\n' "$plot_palette"
    printf 'plot_font_family\t%s\n' "$plot_font_family"
    printf 'plot_base_size\t%s\n' "$plot_base_size"
    printf 'plot_dpi\t%s\n' "$plot_dpi"
    printf 'report_title\t%s\n' "$report_title"
    printf 'report_logo\t%s\n' "$report_logo"
    printf 'report_style\t%s\n' "$report_style"
    printf 'enable_lr_spatial_slam\t%s\n' "$enable_lr_spatial_slam"
    printf 'sample_count\t%s\n' "$sample_count"
    printf 'git_commit\t%s\n' "$git_commit"
    printf 'git_branch\t%s\n' "$git_branch"
    printf 'git_dirty\t%s\n' "$git_dirty"
    printf 'pipeline_release_tag\t%s\n' "$pipeline_release_tag"
} > "$manifest_path"

{
    printf 'metric\tvalue\n'
    printf 'sample_count\t%s\n' "$sample_count"
    printf 'threads\t%s\n' "$threads"
    printf 'memory_gb\t%s\n' "$memory_gb"
    printf 'max_parallel_samples\t%s\n' "$max_parallel_samples"
    printf 'execution_profile\t%s\n' "$execution_profile"
    printf 'resume_mode\t%s\n' "$resume_mode"
} > "$summary_path"

printf 'status=completed\ntimestamp=%s\n' "$timestamp" > "$checkpoint_path"

if [[ "$enable_report" == "true" ]]; then
    python3 scripts/generate_report_skeleton.py \
        --project-name "$project_name" \
        --analysis-type "$analysis_type" \
        --run-id "$run_id" \
        --report-title "$report_title" \
        --report-style "$report_style" \
        --report-logo "$report_logo" \
        --plot-palette "$plot_palette" \
        --plot-font-family "$plot_font_family" \
        --plot-base-size "$plot_base_size" \
        --plot-dpi "$plot_dpi" \
        --sample-count "$sample_count" \
        --config-path "$CONFIG" \
        --samples-path "$SAMPLES" \
        --summary-path "$summary_path" \
        --figures-dir "$figures_dir" \
        --tables-dir "$tables_dir" \
        --reports-dir "$reports_dir" \
        --deliverables-dir "$deliverables_dir" \
        --editable-dir "$editable_dir" >/dev/null
fi

{
    printf '[%s] workflow started\n' "$timestamp"
    printf '[%s] project_name=%s\n' "$timestamp" "$project_name"
    printf '[%s] config=%s\n' "$timestamp" "$CONFIG"
    printf '[%s] samples=%s\n' "$timestamp" "$SAMPLES"
    printf '[%s] manifest=%s\n' "$timestamp" "$manifest_path"
    printf '[%s] summary=%s\n' "$timestamp" "$summary_path"
    printf '[%s] checkpoint=%s\n' "$timestamp" "$checkpoint_path"
    printf '[%s] enable_report=%s\n' "$timestamp" "$enable_report"
    printf '[%s] enable_lr_spatial_slam=%s\n' "$timestamp" "$enable_lr_spatial_slam"
    printf '[%s] reports_dir=%s\n' "$timestamp" "$reports_dir"
    printf '[%s] deliverables_dir=%s\n' "$timestamp" "$deliverables_dir"
    printf '[%s] tmp_dir=%s\n' "$timestamp" "$tmp_dir"
    printf '[%s] template run completed\n' "$timestamp"
} > "$log_path"

printf 'Template run completed.\n'
printf 'Manifest: %s\n' "$manifest_path"
printf 'Log: %s\n' "$log_path"
