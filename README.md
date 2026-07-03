## Pre-install software and R packages

### Software
+ R-4.1.3

### R packages
+ R version 4.1.3
+ Seurat_4.1.0
+ SeuratObject_4.0.4
+ ggplot2_3.3.5
+ dplyr_1.0.9
+ patchwork_1.1.1
+ edgeR_3.36.0
+ clusterProfiler_4.2.0
+ org.Hs.eg.db_3.15.0
+ org.Mm.eg.db_3.15.0
+ ReactomePA_1.38.0
+ msigdbr_7.4.1
+ CellChat_2.1.2
+ presto_1.0.0
+ harmony_0.1.0

## Scripts

### 1. Single-cell RNA-seq analysis

#### Cell-cluster

Seurat.R -- Cluster analysis with Seurat

#### Cellchat-analysis

cellchat.R -- Cell communication analysis with CellChat

#### Diff-expression

edgeR.R -- Differential expression analysis with edgeR

#### scRNA-pipeline

scRNA_pipeline.R -- Single-cell pipeline: QC, integration (harmony/CCA), clustering, UMAP

#### Diff-enrich

scRNA_diff_enrich.R -- Differential expression and enrichment analysis (ORA/GSEA)

##### shell scripts

Rscript 1.Seurat.R
Rscript 2.cellchat.R
Rscript 3.edgeR.R
Rscript 4.scRNA_pipeline.R
Rscript 5.scRNA_diff_enrich.R
