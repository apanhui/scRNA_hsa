## demo.scRNA_pipeline.R - Seurat 3.x + harmony demo version of 4.scRNA_pipeline.R
##
## Mirrors the design of ../code/4.scRNA_pipeline.R (CONFIG block driving
## the same pipeline), but uses only APIs available in Seurat 3.1.1 +
## harmony 1.0 (no SeuratObject dependency, no Seurat 4-only helpers such
## as SplitObject.Image or FindIntegrationAnchors(... normalization.method)).
##
## Usage:
##   Rscript demo.scRNA_pipeline.R <config.r>
## where <config.r> is a small R script that defines CONFIG (see
## ../config/demo_pipeline_config.R for an example).
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript demo.scRNA_pipeline.R <config.r>")
source(args[[1]], local = FALSE)

stopifnot(exists("CONFIG"), is.list(CONFIG))
inputs <- CONFIG$inputs
outdir <- CONFIG$outdir
prefix <- CONFIG$prefix
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

message("[demo.scRNA_pipeline] inputs: ", paste(inputs, collapse = ", "))

## ---- load ----
object.list <- lapply(inputs, function(p) {
  if (grepl("\\.rds$", p, ignore.case = TRUE)) {
    obj <- readRDS(p)
    if (!inherits(obj, "Seurat")) stop("RDS is not a Seurat object: ", p)
    return(obj)
  }
  mat <- Read10X(data.dir = p, gene.column = 1)
  CreateSeuratObject(counts = mat, project = basename(p))
})

## ---- merge + normalise + HVG ----
object <- merge(x = object.list[[1]], y = unlist(object.list[-1]))
object <- NormalizeData(object, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
object <- FindVariableFeatures(object, selection.method = "vst",
                               nfeatures = 2000, verbose = FALSE)

## ---- harmony integration ----
object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object, npcs = 30, verbose = FALSE)
object <- RunHarmony(object, group.by.vars = "orig.ident",
                      reduction = "pca", dims.use = 1:30)

## ---- cluster + UMAP ----
object <- FindNeighbors(object, reduction = "harmony", dims = 1:30, verbose = FALSE)
object <- FindClusters(object, resolution = 0.5, verbose = FALSE)
object <- RunUMAP(object, reduction = "harmony", dims = 1:30, verbose = FALSE)

## ---- save ----
saveRDS(object, file = file.path(outdir, paste0(prefix, ".seurat.rds")))
write.csv(object[[]], file = file.path(outdir, paste0(prefix, ".metadata.csv")),
          row.names = TRUE)

## ---- plots ----
p1 <- DimPlot(object, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave(p1 + p2, file = file.path(outdir, paste0(prefix, ".UMAP.pdf")),
       width = 12, height = 5)

message("[demo.scRNA_pipeline] wrote ", outdir, " (",
        paste(list.files(outdir), collapse = ", "), ")")