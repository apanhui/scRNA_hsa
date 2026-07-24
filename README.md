# Code & Software Submission

The package layout is:

```
├── README.md                 <- you are reading it
├── code/                     <- production scripts
│   ├── 1.Seurat.R
│   ├── 2.cellchat.R
│   ├── 3.edgeR.R
│   ├── 4.scRNA_pipeline.R
│   └── 5.scRNA_diff_enrich.R
└── demo/                     <- example dataset + runnable end-to-end demo
    ├── README.md             <- demo-specific instructions
    ├── run_demo.sh           <- one-shot demo runner (5 steps)
    ├── code/                 <- demo versions of the scripts
    │   ├── demo.Seurat.R, demo.edgeR.R, demo.scRNA_pipeline.R
    │   └── demo.cellchat_init.R, demo.cellchat_sample.R  (CellChat two-step)
    ├── data/                 <- small synthetic 10X dataset (400 cells x 1000 genes)
    ├── results/              <- populated when demo is run
    └── expected_output/      <- saved summaries for reviewer byte-diff
```

------------------------------------------------------------------------

## 1. System requirements

**Operating system:**

* CentOS 7.9.2009

**Software:**

* **R ≥ 4.1.3** 
* **Java 8+** for clusterProfiler / ReactomePA enrichment plotting.

**R packages:**

All analysis is implemented in **R 4.1.3** with the following R packages
(`code/README.md` lists them in detail):

| Package         | Version (tested)   | Role                                     |
|-----------------|--------------------|------------------------------------------|
| Seurat          | 4.1.0              | QC, integration, clustering, DE          |
| SeuratObject    | 4.0.4              | Data structure backing Seurat            |
| harmony         | 0.1.0              | Multi-sample integration                 |
| edgeR           | 3.36.0             | Bulk-style DE analysis                   |
| clusterProfiler | 4.2.0              | ORA / GSEA enrichment                    |
| org.Hs.eg.db    | 3.15.0             | Human gene annotations                   |
| org.Mm.eg.db    | 3.15.0             | Mouse gene annotations                   |
| ReactomePA      | 1.38.0             | Reactome enrichment                      |
| msigdbr         | 7.4.1              | MSigDB Hallmark gene sets                |
| CellChat        | 2.1.2              | Cell-cell communication                  |
| presto          | 1.0.0              | Fast DE testing                          |
| ggplot2         | 3.3.5              | Plotting                                 |
| dplyr           | 1.0.9              | Data wrangling                           |
| patchwork       | 1.1.1              | Multi-panel plots                        |


------------------------------------------------------------------------


## 2 Installation guide

```bash
# 1. Install R 4.1.3 (skip if already present; conda works too)
#    See https://cran.r-project.org/

# 2. Install the R packages listed in §1.1. We recommend renv:
cd code
R -e 'install.packages("renv"); renv::restore()'

# If you prefer Bioconductor directly:
R -e 'install.packages("BiocManager")'
R -e 'BiocManager::install(c(
        "Seurat", "SeuratObject", "harmony", "edgeR",
        "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db",
        "ReactomePA", "msigdbr", "CellChat", "presto",
        "ggplot2", "dplyr", "patchwork"))'
```

**Typical install time.** On a "normal" desktop (8-core, 16 GB, fresh
conda env):

| Step                                      | Time     |
|-------------------------------------------|----------|
| R + BiocManager install                   | 5–10 min |
| `BiocManager::install(...)`               | 30–45 min|
| **Total**                                 | **≈ 35–55 min** |

## 3 Demo
### 3.1 Dataset of demo

The `demo/data/` directory ships a **simulated** 10X dataset:

* 2 conditions × 2 replicates × 100 cells = **400 cells**
* **800 genes** including real cell-type markers (CD3D, CD14, EPCAM,
  COL1A1, …), interferon-stimulated genes (ISG15, MX1, IFIT1, …), and
  400+ synthetic `GENE_NNNN` placeholders so reviewers can see the
  noise floor.
* 4 ground-truth cell types per condition: **T_cell, Myeloid,
  Epithelial, Stromal** (with proportions skewed toward Myeloid in the
  PEDV condition to mimic infection-driven expansion).
* The dataset is provided in two forms:

  * `demo/data/raw/{Control,PEDV}_rep{1,2}/` — four per-replicate 10X
    folders, used by script 4 (`scRNA_pipeline.R`).
  * `demo/data/merged/{Control,PEDV}/` — replicates concatenated per
    condition, used by script 1 (`Seurat.R`).
  * `demo/data/Myeloid.exp.xls` — gene × sample count matrix restricted
    to the Myeloid lineage, used by script 3 (`edgeR.R`).

### 3.2 Run on data

```bash
cd demo
bash run_demo.sh
```

The script runs five R scripts in sequence and writes their outputs
to `demo/results/`:

| Step                       | Script                            | R env          | Output dir                    |
|----------------------------|-----------------------------------|----------------|-------------------------------|
| Seurat CCA integration     | `demo/code/demo.Seurat.R`         | R 4.1.3        | `demo/results/seurat/`        |
| harmony integration        | `demo/code/demo.scRNA_pipeline.R` | R 4.1.3        | `demo/results/scRNA_pipeline/`|
| edgeR DE on Myeloid subset | `demo/code/demo.edgeR.R`          | R 4.1.3        | `demo/results/edger/`         |
| CellChat: build init.Rds   | `demo/code/demo.cellchat_init.R`  | R 4.1.3        | `demo/results/cellchat/`      |
| CellChat: run inference    | `demo/code/demo.cellchat_sample.R`| R 4.1.3        | `demo/results/cellchat/`      |

### 3.3 Expected output

```
demo/results/
├── seurat/
│   ├── UMAP.pdf              <- UMAP coloured by sample / cluster
│   ├── TSNE.pdf              <- t-SNE coloured by sample / cluster
│   ├── obj.Rda               <- Seurat object (RNA assay + integrated reduction)
│   ├── markers.Rda           <- FindAllMarkers output
│   ├── cluster_markers.tsv   <- top-5 markers per cluster (text)
│   └── run_summary.txt       <- counts of cells / clusters / markers
├── scRNA_pipeline/
│   ├── demo_result.seurat.rds
│   ├── demo_result.metadata.csv
│   └── demo_result.UMAP.pdf
└── edger/
│   ├── Myeloid.demo.all.xls     <- DE table for all expressed genes
│   ├── Myeloid.demo.filter.xls  <- DE table filtered at |log2FC| > log2(1.5) & P < 0.05
│   └── Myeloid.demo.summary.txt
└── cellchat/
    ├── cellchat.init.Rds     <- intermediate CellChat object
    ├── CellChat.LR.all.xls   <- inferred LR table (all)
    ├── CellChat.LR.filtered.xls <- inferred LR table (P < 0.05)
    ├── Pathway.CommunProb.xls <- aggregated pathway x cell-pair table
    ├── pathways.show.list    <- top-5 pathways (one per line)
    ├── Pathway.<name>.hierarchy.pdf / .networks.pdf / .heatmap.pdf  (5 pathways x 3 layouts)
    ├── cellchat.Rds          <- final CellChat object
    └── cellchat_summary.txt
```
### 3.4 Expected run time

| Step                       | Wall-clock |
|----------------------------|-----------:|
| Seurat CCA integration     | ≈ 90 s     |
| harmony integration        | ≈ 60 s     |
| edgeR DE                   | < 5 s      |
| CellChat init + inference  | ≈ 60 s     |
| **Total**                  | **≈ 5 min**|
## 4 Instructions for use

### 4.1 Running the pipeline on your own 10X data

For each `Seurat.R`-style integration run:

1. Place your 10X folders under a single parent directory, e.g.
   `my_data/Control1/`, `my_data/Control2/`, `my_data/Treated1/`, …
   Each subfolder must contain `barcodes.tsv`, `genes.tsv`,
   `matrix.mtx` in 10X v2 format.
2. Edit `code/1.Seurat.R`: change the `data_name` vector to the names
   of your conditions, then run `Rscript code/1.Seurat.R`.
3. For pipeline-style runs that integrate many samples via harmony or
   CCA, edit the `CONFIG$inputs` vector in `code/4.scRNA_pipeline.R`
   to point at each 10X folder and run `Rscript code/4.scRNA_pipeline.R`.

The script produces `obj.Rda` (or `<prefix>.seurat.rds`), a set of
`UMAP.pdf` / `TSNE.pdf` plots, and a `markers.Rda` table of cluster
markers.

### 4.2 Differential expression on a single lineage

`code/3.edgeR.R` is the bulk-style DE script. It expects a tab-separated
count matrix (`Myeloid.exp.xls`) whose first `C_number+T_number` columns
are raw counts (rows = genes). Adjust the four header variables at the
top of the file (`count_file`, `C_number`, `T_number`, `group_colnames`)
and run `Rscript code/3.edgeR.R`.

### 4.3 Cell-cell communication

`code/2.cellchat.R` expects a Seurat object converted to a CellChat
object (`cellchat.init.Rds`). 

### 4.4 Differential expression + enrichment per cluster

`code/5.scRNA_diff_enrich.R` reads `<prefix>.seurat.rds`, runs DE
(presto or FindMarkers), then ORA/GSEA enrichment against GO, KEGG,
Reactome, and Hallmark. 
