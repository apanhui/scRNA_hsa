# Demo dataset & end-to-end run

This directory contains everything a reviewer needs to verify the
analysis pipeline runs end-to-end on a small synthetic dataset, without
needing access to the manuscript's sequencing files.

## Quick start

```bash
cd demo
bash run_demo.sh
```

That's it. `run_demo.sh` will:

1. Read the merged 10X data (`demo/data/merged/{Control,PEDV}/`) and
   run a Seurat pipeline with CCA integration
   (`demo/code/demo.Seurat.R`). Output goes to `demo/results/seurat/`.
2. Read the per-replicate 10X data
   (`demo/data/raw/{Control_rep1,Control_rep2,PEDV_rep1,PEDV_rep2}/`)
   and run a Seurat pipeline with harmony integration
   (`demo/code/demo.scRNA_pipeline.R`). Output goes to
   `demo/results/scRNA_pipeline/`.
3. Run edgeR DE on the Myeloid count matrix
   (`demo/data/Myeloid.exp.xls`) using
   (`demo/code/demo.edgeR.R`). Output goes to `demo/results/edger/`.

Total wall-clock on a normal desktop: ~3 minutes.

## What the data is

Structure of the synthetic dataset:


| File / folder                              | What it is                              | Used by script      |
|--------------------------------------------|-----------------------------------------|---------------------|
| `demo/data/raw/Control_rep1/` (etc., 4 folders) | per-replicate 10X (100 cells each) | `demo.scRNA_pipeline.R` |
| `demo/data/merged/Control/` (etc., 2 folders)    | merged 10X (200 cells each)         | `demo.Seurat.R`     |
| `demo/data/Myeloid.exp.xls`                     | gene × sample count matrix (Myeloid subset) | `demo.edgeR.R`  |

Ground truth:

* 2 conditions: Control vs PEDV
* 4 cell types per condition: T_cell, Myeloid, Epithelial, Stromal
  (with Myeloid expanded in PEDV to mimic infection-driven
  expansion)
* 30 interferon-stimulated genes (ISG15, MX1, …) up-regulated in
  PEDV-Myeloid cells

## Expected outputs

After `run_demo.sh` completes, see:

| File                                       | What it shows                         |
|--------------------------------------------|----------------------------------------|
| `results/seurat/run_summary.txt`           | cell / cluster / marker counts        |
| `results/seurat/cluster_markers.tsv`       | top 5 markers per cluster             |
| `results/seurat/UMAP.pdf`                  | UMAP coloured by sample + cluster      |
| `results/seurat/TSNE.pdf`                  | t-SNE coloured by sample + cluster     |
| `results/seurat/obj.Rda`                   | Seurat object                          |
| `results/seurat/markers.Rda`               | FindAllMarkers result                  |
| `results/scRNA_pipeline/demo_result.seurat.rds` | Seurat object with harmony reduction |
| `results/scRNA_pipeline/demo_result.metadata.csv` | per-cell metadata                   |
| `results/scRNA_pipeline/demo_result.UMAP.pdf`     | UMAP coloured by sample + cluster |
| `results/edger/Myeloid.demo.summary.txt`   | edgeR DE summary                       |
| `results/edger/Myeloid.demo.all.xls`       | DE table (all expressed genes)         |
| `results/edger/Myeloid.demo.filter.xls`    | DE table filtered at |log2FC| > log2(1.5), P < 0.05 |
