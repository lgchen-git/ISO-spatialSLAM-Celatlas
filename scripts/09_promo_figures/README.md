# Optional PPT promotion figure pack

This optional delivery stage reads completed `01_read_qc` through
`07_slam_joint` outputs and writes a reproducible figure pack under the sample
result root. It does not rerun alignment, transcript discovery, clustering, or
SLAM calling.
The delivery directory contains single-panel PNG/PDF figures; composite panels
are not retained there.

The unified entry point is:

```bash
bash workflow/run.sh --config conf/<project>.yaml --promo-figures
```

The pack explicitly records missing H&E, controls, replicate metadata, and
read-to-isoform mappings. Unlabelled single-sample T>C results remain
exploratory.

This is a post-report delivery module, not a formal numbered stage in
`analysis_layout.yaml`. The normal analysis and final report end at
`08_report`; this module is run explicitly with `--promo-figures`.

The optional `promo_isoform_spatial_enabled` stage reconstructs conservative
read-to-isoform assignments from the primary genomic BAM, FLAMES exon models,
and the barcode/UMI per-read table. It writes independent candidate-gene,
isoform-spatial, and structure-event figures plus
`isoform_figure_manifest.tsv`. Ambiguous reads are excluded. If no pair meets
the configured thresholds, no isoform A/B figure is emitted and
`05_metadata/isoform_spatial_qc.json` records the reason.
