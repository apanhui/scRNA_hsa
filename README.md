# Code & Software Submission

The package layout is:

```
├── README.md                 <- you are reading it
├── LICENSE                   <- MIT license (OSI-approved)
├── code/                     <- production scripts
│   ├── 1.Seurat.R
│   ├── 2.cellchat.R
│   ├── 3.edgeR.R
│   ├── 4.scRNA_pipeline.R
│   └── 5.scRNA_diff_enrich.R
└── demo/                     <- example dataset + runnable end-to-end demo
    ├── README.md             <- demo-specific instructions
    ├── run_demo.sh           <- one-shot demo runner (5 steps)
    ├── config/               <- CONFIG blocks consumed by demo.scRNA_pipeline.R
    ├── code/                 <- demo versions of the scripts (R-3.x + R413 compatible)
    │   ├── demo.Seurat.R, demo.edgeR.R, demo.scRNA_pipeline.R
    │   ├── demo.cellchat_init.R, demo.cellchat_sample.R  (CellChat two-step)
    │   └── gen_demo_data.py + merge_10x.py  (dataset regenerators)
    ├── data/                 <- small synthetic 10X dataset (400 cells x 1000 genes)
    ├── results/              <- populated when demo is run
    └── expected_output/      <- saved summaries for reviewer byte-diff
```

------------------------------------------------------------------------

## 1. Required content (per Nature Research checklist)

### 1.1 Compiled standalone software and/or source code

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

### 1.2 A small (simulated or real) dataset to demo the software/code

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

The synthetic generator is shipped as `demo/code/gen_demo_data.py` and
`demo/code/merge_10x.py` so reviewers can recreate the dataset with a
different random seed or change cell counts.

The simulated data is **not** a substitute for the manuscript's
biological data — it is small enough (≈ 1 MB) to keep the demo under
three minutes and structured enough that the published workflow
recovers the expected cell types and the interferon-stimulated gene
signature in PEDV (see §3 below).

------------------------------------------------------------------------

## 2. README file

### 2.1 System requirements

**Operating system.** Tested on:

* CentOS 6.10 (kernel 2.6.32, the server used for the manuscript runs).
* Ubuntu 20.04 / 22.04 LTS (x86_64) for cross-platform validation.
* macOS 12+ (Apple Silicon works once R / Bioconductor binaries are
  available; not the primary target).

**Hardware.** A "normal" desktop computer is sufficient:

* 8 GB RAM minimum (16 GB recommended for the full single-cell run on
  real data — the demo itself only needs ~2 GB).
* ≥ 4 CPU cores recommended (harmony / CCA integration are
  multi-threaded). No GPU required. No non-standard hardware.

**Software.**

* **R ≥ 4.1.0** . Earlier versions (3.6.x) will fail on Seurat ≥ 4.x.
* **Java 8+** for clusterProfiler / ReactomePA enrichment plotting.

The Seurat 3.x compatible demo scripts (`demo/code/demo.*.R`) have
additionally been validated against the conda R 3.5.1 environment that
ships with Seurat 3.1.1 (`/Bio/bin/Rscript-3.5.1_conda` on the
manuscript server). See *Known issues* below for what changes.

### 2.2 Installation guide

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

### 2.3 Demo

The demo runs end-to-end on a 400-cell synthetic dataset. See
`demo/README.md` for the full instructions; in short:

```bash
cd demo
bash run_demo.sh
```

The script runs five R scripts in sequence and writes their outputs
to `demo/results/`:

| Step                       | Script                            | R env              | Output dir                    |
|----------------------------|-----------------------------------|--------------------|-------------------------------|
| Seurat CCA integration     | `demo/code/demo.Seurat.R`         | R 3.5.1_conda      | `demo/results/seurat/`        |
| harmony integration        | `demo/code/demo.scRNA_pipeline.R` | R 3.5.1_conda      | `demo/results/scRNA_pipeline/`|
| edgeR DE on Myeloid subset | `demo/code/demo.edgeR.R`          | R 3.5.1_conda      | `demo/results/edger/`         |
| CellChat: build init.Rds   | `demo/code/demo.cellchat_init.R`  | R 4.1.3 (R413)     | `demo/results/cellchat/`      |
| CellChat: run inference    | `demo/code/demo.cellchat_sample.R`| R 4.1.3 (R413)     | `demo/results/cellchat/`      |

The first three steps run with `/Bio/bin/Rscript-3.5.1_conda` (Seurat
3.x compatible). The CellChat steps run with
`/public2/Bio/pipeline/Toolkit/miniconda/miniconda3/envs/R413/bin/Rscript`
because CellChat (>= 1.6) is not installed in the 3.5.1 conda env.
Override the interpreters via `R_BIN=` / `R413_BIN=` env vars, or
`RUN_CELLCHAT=0` to skip the CellChat step.

**Expected output :**

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

A successful run prints the following counts in
`demo/results/seurat/run_summary.txt`:

```
cells_total       : 400
genes_total       : 1000
samples           : Control,PEDV
n_clusters        : 4
n_markers_total   : ~80
top_marker_example: cluster 0 -> CD2     # T-cell marker
                    cluster 1 -> CSF1R   # Myeloid marker
                    cluster 2 -> EPCAM   # Epithelial marker
                    cluster 3 -> COL1A1 # Stromal marker
```

…and in `demo/results/edger/Myeloid.demo.summary.txt`:

```
input_genes           : 1000
input_genes_after_flt : ~990
DE_genes_sig          : ~80            # ISGs and inflammation markers dominate
top_up_regulated:  NMI, IFIT2, ISG15, … # ISGs up in PEDV
top_down_regulated: small set of synthetic genes
```

…and in `demo/results/cellchat/cellchat_summary.txt`:

```
cells              : 400
cell_types         : C0,C1,C2,C3
n_lr_pairs_total   : 68
n_lr_pairs_sig     : 68
n_pathways_visualised: 5
pathway_examples   : CCL, COLLAGEN, GALECTIN, TNF, VEGF
```

The recovery of CD2 / CSF1R / EPCAM / COL1A1 as cluster-defining
markers, NMI / IFIT2 / ISG15 as the top PEDV-up-regulated genes, and
the canonical LR pathways (CCL, TNF, VEGF, …) confirms that the full
pipeline — including the CellChat cell-cell communication step — is
wired correctly end-to-end.

**Expected run time.** On a "normal" desktop (8-core, 16 GB):

| Step                       | Wall-clock |
|----------------------------|-----------:|
| Seurat CCA integration     | ≈ 90 s     |
| harmony integration        | ≈ 60 s     |
| edgeR DE                   | < 5 s      |
| CellChat init + inference  | ≈ 60 s     |
| **Total**                  | **≈ 5 min**|

### 2.4 Instructions for use

#### 2.4.1 Running the pipeline on your own 10X data

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

#### 2.4.2 Differential expression on a single lineage

`code/3.edgeR.R` is the bulk-style DE script. It expects a tab-separated
count matrix (`Myeloid.exp.xls`) whose first `C_number+T_number` columns
are raw counts (rows = genes). Adjust the four header variables at the
top of the file (`count_file`, `C_number`, `T_number`, `group_colnames`)
and run `Rscript code/3.edgeR.R`.

#### 2.4.3 Cell-cell communication

`code/2.cellchat.R` expects a Seurat object converted to a CellChat
object (`cellchat.init.Rds`). Follow the [CellChat vignette][1] to
produce that RDS, then run the script.

#### 2.4.4 Differential expression + enrichment per cluster

`code/5.scRNA_diff_enrich.R` reads `<prefix>.seurat.rds`, runs DE
(presto or FindMarkers), then ORA/GSEA enrichment against GO, KEGG,
Reactome, and Hallmark. Configuration is in the `CONFIG` block at the
top of the file.

[1]: https://github.com/sqjin/CellChat
