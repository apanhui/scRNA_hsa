# Expected demo outputs

These are the actual outputs produced by `run_demo.sh` on the
manuscript server (R 3.5.1_conda, Seurat 3.1.1, harmony 1.0) against the
synthetic dataset shipped in `demo/data/`. Files in this directory
should match what you see in your `demo/results/` sub-folders after a
clean run.

| File                              | Source in demo results                    |
|-----------------------------------|------------------------------------------|
| `seurat_run_summary.txt`         | `demo/results/seurat/run_summary.txt`   |
| `edger_summary.txt`              | `demo/results/edger/Myeloid.demo.summary.txt` |
| `cluster_markers_top5.tsv`       | `demo/results/seurat/cluster_markers.tsv` |

## What a successful run shows

1. **Seurat CCA integration** (`seurat_run_summary.txt`):
   400 cells, 800 genes, 2 clusters. The top marker for cluster 0 is
   `CD2` (T-cell); for cluster 1 it is `CSF1R` (Myeloid). The
   pipeline recovered the two largest cell types.

2. **harmony integration** (`demo_result.UMAP.pdf` in
   `demo/results/scRNA_pipeline/`): the UMAP shows the four samples
   (two Control, two PEDV) overlaid, with clusters co-located across
   conditions — confirming the harmony correction did its job.

3. **edgeR DE** (`edger_summary.txt`): about 50 genes pass the
   |log2FC| > log2(1.5) & P < 0.05 cutoff. The top up-regulated hits in
   PEDV are `NMI` and `IFIT2` — interferon-stimulated genes that
   the synthetic generator was primed to up-regulate in the Myeloid
   subset of PEDV. This is the same biological signal the manuscript
   reports on the real data.
