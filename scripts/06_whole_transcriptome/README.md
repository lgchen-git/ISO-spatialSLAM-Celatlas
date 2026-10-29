# 06 whole transcriptome

Builds a spatial gene-expression matrix from all long-read molecules and applies
the manually configured Scanpy/Leiden clustering parameters.

After clustering, `summarize_spatial_qc.py` writes bin50 UMI, long-read
molecule, detected-gene and exploratory dynamic-proxy statistics. Isoform
counts are sample-level FLAMES candidate counts; per-bin isoform counts require
the optional read-to-isoform module.
