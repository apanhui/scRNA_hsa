## demo.cellchat_init.R - Step 1 of 2 of the CellChat demo
##
## Reads the Seurat object produced by demo.Seurat.R (under
## demo/results/seurat/obj.Rda) and creates the cellchat.init.Rds file
## that cellchat_sample.R expects as input. Mirrors the role of
## ../seurat2cellchat.R in the production pipeline but stays as short
## as possible: human CellChatDB, no species auto-detect, no label
## remapping beyond the "C0/C1" rule CellChat requires.
##
## Usage:
##   Rscript demo.cellchat_init.R <seurat_obj.rda> <out_dir>
##
## Output:
##   <out_dir>/cellchat.init.Rds   - CellChat object ready for inference
##
## Requires CellChat >= 1.5.

suppressMessages({
  library(Seurat)
  library(CellChat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript demo.cellchat_init.R <seurat_obj.rda> <out_dir>")
}
seurat_obj_path <- args[[1]]
out_dir <- args[[2]]
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## demo.Seurat.R saved with save(object, ...); load() assigns it as
## a top-level variable named "object".
obj_name <- load(seurat_obj_path)
obj <- get(obj_name)
rm(list = obj_name)
stopifnot(inherits(obj, "Seurat"))

## CellChat rejects labels that start with a digit ("0", "1" -> error).
obj$cell_type <- factor(paste0("C", obj$seurat_clusters))
Idents(obj) <- "cell_type"

cellchat <- createCellChat(object = obj, group.by = "cell_type")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)

## Minimum options block: the sample step will pick these up.
cellchat@options$select_commun <- list(
  source = NULL, target = NULL, pairs = NULL, drop = FALSE, method = "default"
)

out_path <- file.path(out_dir, "cellchat.init.Rds")
saveRDS(cellchat, file = out_path)

message("[demo.cellchat_init] wrote ", out_path)
message("[demo.cellchat_init] cell types: ",
        paste(levels(cellchat@idents), collapse = ", "),
        " (n=", ncol(cellchat), " cells, ",
        nrow(cellchat@data.signaling), " signaling genes)")