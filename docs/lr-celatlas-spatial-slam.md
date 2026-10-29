# LR-Celatlas spatial SLAM workflow

This workflow is configuration-driven and runs without interactive or AI assistance.

## Prepare a project

1. Copy `conf/lr_celatlas_spatial_slam.template.yaml` to `conf/<project>.yaml`.
2. Copy `conf/samples.lr_celatlas_spatial_slam.template.tsv` to `conf/<project>.samples.tsv`.
3. Replace every `CHANGE_ME` value and set `samples_tsv` to the copied sample sheet.
4. Validate the configuration and input contract:

```bash
make validate CONFIG=conf/<project>.yaml
bash workflow/run.sh --config conf/<project>.yaml --dry-run
```

## Run

```bash
bash workflow/run.sh --config conf/<project>.yaml
```

The eight stages are:

1. `01_read_qc`
2. `02_barcode_cdna`
3. `03_align`
4. `04_full_length_transcriptome`
5. `05_slam`
6. `06_whole_transcriptome`
7. `07_slam_joint`
8. `08_report`

FLAMES and gffcompare provide full-length transcript discovery and reference
comparison. The whole-transcriptome stage clusters all long-read molecules with
the parameters entered in the project config. It does not scan parameter grids.
The joint stage overlays SLAM classes on the same whole-transcriptome clusters.

The normal analysis ends at `08_report`. The optional `09_promo_figures`
directory is a post-report delivery module and is not part of the formal stage
sequence or completion markers. It is enabled only with `--promo-figures`.

Resume from a stage or rebuild one stage with:

```bash
bash workflow/run.sh --config conf/<project>.yaml --from 05_slam --to 08_report
bash workflow/run.sh --config conf/<project>.yaml --from 08_report --to 08_report --force-stage 08_report
bash workflow/run.sh --config conf/<project>.yaml --from 08_report --to 08_report \
  --force-stage 08_report --report-title "长读长空间转录组与 SLAM 分析报告"
```

All scientific parameters are explicit top-level config keys. Edit the config,
then rerun from the affected stage. The stage marker stores the config SHA-256;
changed configurations do not silently reuse a previous completion marker.

To update only figures and the final report from existing stage outputs:

```bash
bash workflow/run.sh --config conf/<project>.yaml --refresh-visuals
```

This mode redraws QC, substitution, SLAM and spatial figures, then rebuilds and
validates the offline report. It does not rerun alignment, molecule calling,
clustering, UMAP or statistical analysis. It requires the compact source tables,
stage metadata and spatial H5AD already present under the sample result tree.
`--refresh-visuals` cannot be combined with `--from`, `--to` or `--force-stage`.
Plot font, base size and resolution are controlled by `plot_font_family`,
`plot_base_size` and `plot_dpi` in the project config. `plot_hist_max_rows`
limits the deterministic cDNA-length display sample while stage summaries retain
the full-read counts.

Every stage records its command log, software versions, Git provenance, output manifest and completion marker. The final report is:

```text
<results_root>/<sample_id>/08_report/<sample_id>_report.html
```

## Sample sheet

Required columns:

```text
sample_id fastq chip_id whitelist_csv species chemistry genome_fa gtf
```

Optional columns:

```text
mmi junc_bed threads
```

Optional experiment-design columns are also accepted:

```text
tissue_type experimental_group label_4su_dose label_4su_duration harvest_time
unlabelled_control biological_replicate section_id tissue_image
```

They generate warnings when absent or empty, but do not block legacy projects.

Paths must point to existing nonempty files. A chip ID may be assigned to only one sample in a run.

## Validation

```bash
make lint
make smoke
python3 scripts/validate_delivery_report.py \
  --report <results_root>/<sample_id>/08_report/<sample_id>_report.html \
  --manifest <results_root>/<sample_id>/08_report/05_metadata/delivery_manifest.tsv
```

To build the separate PPT and exploratory structure package after the normal
report:

```bash
bash workflow/run.sh --config conf/<project>.yaml --promo-figures
```
