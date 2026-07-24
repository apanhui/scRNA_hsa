## demo.Seurat.R - Seurat compatible demo version of 1.Seurat.R
##
## This script demonstrates end-to-end usage of the analysis pipeline on
## the small simulated dataset shipped in /demo/data/merged.
##
## Usage:
##   Rscript demo.Seurat.R <input_dir> <output_dir>
## where <input_dir> contains Control/ and PEDV/ subfolders in 10X
## format (matrix.mtx, barcodes.tsv, genes.tsv).
##
## Outputs (written under <output_dir>):
##   - obj.Rda                 Seurat object (RNA assay, integrated reduction)
##   - markers.Rda             FindAllMarkers result
##   - UMAP.pdf                UMAP coloured by sample / cluster
##   - TSNE.pdf                t-SNE coloured by sample / cluster
##   - cluster_markers.tsv     top-5 markers per cluster (text summary)
##   - run_summary.txt         one-line summary of clustering outcome
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript demo.Seurat.R <input_dir> <output_dir>")
}
input_dir  <- args[[1]]
output_dir <- args[[2]]
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

data_name <- c("Control", "PEDV")

## ---- read input ----
message("[demo.Seurat] reading 10X data from ", input_dir)
object.list <- list()
for (i in seq_along(data_name)) {
  mat <- Read10X(data.dir = file.path(input_dir, data_name[i]), gene.column = 1)
  object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i])
}

object <- merge(x = object.list[[1]], y = unlist(object.list[-1]),
                add.cell.ids = data_name)
message("[demo.Seurat] merged: ", ncol(object), " cells, ", nrow(object), " genes")

## ---- normalisation and feature selection ----
object <- NormalizeData(object, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
object <- FindVariableFeatures(object, selection.method = "vst",
                               nfeatures = 2000, verbose = FALSE)

## ---- Seurat 3.x integration via CCA anchors ----
object.list <- SplitObject(object, split.by = "orig.ident")
for (i in seq_along(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = FALSE)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]],
                                           selection.method = "vst",
                                           nfeatures = 2000, verbose = FALSE)
}
anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30)
object <- IntegrateData(anchorset = anchors, dims = 1:30)

## ---- downstream ----
DefaultAssay(object) <- "integrated"
object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object, npcs = 30, verbose = FALSE)
object <- RunUMAP(object, dims = 1:30, verbose = FALSE)
object <- RunTSNE(object, dims = 1:30, verbose = FALSE)
object <- FindNeighbors(object, reduction = "pca", dims = 1:30, verbose = FALSE)
object <- FindClusters(object, resolution = 0.5, verbose = FALSE)

## ---- plots ----
n_clusters <- length(unique(object$seurat_clusters))
cluster_levels <- as.character(sort(unique(object$seurat_clusters)))
object@misc$color.cluster <- rainbow(n_clusters)
names(object@misc$color.cluster) <- cluster_levels

p1 <- DimPlot(object, reduction = "tsne", group.by = "orig.ident",
              cols = object@misc$color.cluster, label = FALSE)
p2 <- DimPlot(object, reduction = "tsne", group.by = "seurat_clusters",
              cols = object@misc$color.cluster, label = TRUE)
ggsave(p1 + p2, file = "TSNE.pdf", width = 12, height = 5)

p1 <- DimPlot(object, reduction = "umap", group.by = "orig.ident",
              cols = object@misc$color.cluster, label = FALSE)
p2 <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters",
              cols = object@misc$color.cluster, label = TRUE)
ggsave(p1 + p2, file = "UMAP.pdf", width = 12, height = 5)

## ---- save object ----
DefaultAssay(object) <- "RNA"
save(object, file = "obj.Rda")

## ---- markers ----
Idents(object) <- "seurat_clusters"
obj.markers <- FindAllMarkers(object = object, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25,
                              return.thresh = 0.01, pseudocount.use = 0)
save(obj.markers, file = "markers.Rda")

## Top-5 per cluster for human-readable summary.
top_markers <- obj.markers %>%
  group_by(cluster) %>%
  top_n(5, avg_logFC)
write.table(top_markers, file = "cluster_markers.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)

## ---- run summary ----
summary_lines <- c(
  sprintf("cells_total       : %d", ncol(object)),
  sprintf("genes_total       : %d", nrow(object)),
  sprintf("samples           : %s", paste(unique(object$orig.ident), collapse = ",")),
  sprintf("n_clusters        : %d", n_clusters),
  sprintf("n_markers_total   : %d", nrow(obj.markers)),
  sprintf("top_marker_example: cluster %s -> %s",
          as.character(top_markers$cluster[1]),
          as.character(top_markers$gene[1]))
)
writeLines(summary_lines, "run_summary.txt")
message("[demo.Seurat] done. outputs:\n  ",
        paste(list.files(output_dir), collapse = ", "))