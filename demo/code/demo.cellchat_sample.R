## demo.cellchat_sample.R - Step 2 of 2 of the CellChat demo
##
## Reads <out_dir>/cellchat.init.Rds (produced by demo.cellchat_init.R)
## and runs the same inference block as the production cellchat_sample.R.
##
## Usage:
##   Rscript demo.cellchat_sample.R <work_dir>
## where <work_dir> contains cellchat.init.Rds and is also where outputs
## will be written.
##
## Requires CellChat >= 1.5.

suppressMessages({
  library(CellChat)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript demo.cellchat_sample.R <work_dir>")
work_dir <- args[[1]]
setwd(work_dir)

cellchat <- readRDS("cellchat.init.Rds")

## ---- inference ----
cellchat <- identifyOverExpressedInteractions(
  cellchat, features = rownames(cellchat@data.signaling))
if (nrow(cellchat@LR$LRsig) == 0) {
  stop("No significant LR interactions; check the input CellChat object.")
}
cellchat@idents <- droplevels(cellchat@idents)
population.size <- if (do.call("/", as.list(rev(range(table(cellchat@idents))))) > 20) TRUE else FALSE

## CellChat >= 1.6 replaces SelectCommunication() with filterCommunication()
## (drops low-coverage cell groups) + aggregateNet().
cellchat <- computeCommunProb(cellchat, population.size = population.size)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- aggregateNet(cellchat, thresh = 0.05)

## ---- LR tables ----
## CellChat >= 1.6 dropped GetOnlineData(); subsetCommunication() returns
## the inferred LR table directly. The result already has a `pval`
## column; threshold to P < 0.05 the same way as the production script.
data_all <- as.data.frame(subsetCommunication(cellchat))
write.table(data_all, file = "CellChat.LR.all.xls",
            sep = "\t", quote = FALSE, row.names = FALSE)
if ("pval" %in% colnames(data_all)) {
  data_sig <- subset(data_all, pval <= 0.05)
} else {
  data_sig <- data_all
}
write.table(data_sig, file = "CellChat.LR.filtered.xls",
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- pathway-level inference ----
cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)
## CellChat >= 1.6 only populates netP$prob (per-pathway cell-pair
## probability); netP$pval was dropped. We aggregate the same way the
## production script did: for each pathway, sum probability over all
## cell-pair entries (a pathway is "active" if any cell-pair has
## probability > 0).
pathway_tbl <- data.frame()
if (!is.null(cellchat@netP$prob)) {
  prob_arr <- cellchat@netP$prob      # sources x targets x pathways
  if (length(dim(prob_arr)) >= 3) {
    n_sources <- dim(prob_arr)[1]
    n_targets <- dim(prob_arr)[2]
    pathways  <- cellchat@netP$pathways
    sources   <- levels(cellchat@idents)
    rows <- list()
    for (i in seq_along(pathways)) {
      for (s in seq_len(n_sources)) for (t in seq_len(n_targets)) {
        v <- prob_arr[s, t, i]
        if (!is.na(v)) {
          rows[[length(rows) + 1]] <- data.frame(
            pathway_name = pathways[i],
            source = sources[s],
            target = sources[t],
            prob   = v,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    pathway_tbl <- do.call(rbind, rows)
    write.table(pathway_tbl, file = "Pathway.CommunProb.xls",
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
if (nrow(pathway_tbl) > 0) {
  pathways.show <- (pathway_tbl %>%
                      group_by(pathway_name) %>%
                      summarise(prob = sum(prob), .groups = "drop") %>%
                      arrange(-prob) %>%
                      head(5))[["pathway_name"]]
} else {
  pathways.show <- character()
}
writeLines(pathways.show, "pathways.show.list")

## ---- visualise top-5 pathways ----
for (i in pathways.show) {
  pdf(paste0("Pathway.", i, ".hierarchy.pdf"), 6, 6)
  netVisual_aggregate(cellchat, signaling = i, layout = "hierarchy",
                      vertex.receiver = 1:round(nlevels(cellchat@idents) / 2))
  dev.off()

  pdf(paste0("Pathway.", i, ".networks.pdf"), 6, 6)
  netVisual_aggregate(cellchat, signaling = i, layout = "circle")
  dev.off()

  pdf(paste0("Pathway.", i, ".heatmap.pdf"), 6, 6)
  ht <- netVisual_heatmap(cellchat, signaling = i,
                          color.heatmap = c("white", "darkred"))
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
  dev.off()
}

saveRDS(cellchat, file = "cellchat.Rds")

## ---- summary ----
n_lr_all <- nrow(read.table("CellChat.LR.all.xls", header = TRUE, sep = "\t"))
n_lr_sig <- nrow(read.table("CellChat.LR.filtered.xls", header = TRUE, sep = "\t"))
summary_lines <- c(
  sprintf("cells              : %d", ncol(cellchat)),
  sprintf("cell_types         : %s",
          paste(levels(cellchat@idents), collapse = ",")),
  sprintf("n_lr_pairs_total   : %d", n_lr_all),
  sprintf("n_lr_pairs_sig     : %d", n_lr_sig),
  sprintf("n_pathways_visualised: %d", length(pathways.show)),
  sprintf("pathway_examples   : %s",
          paste(head(pathways.show, 3), collapse = ", "))
)
writeLines(summary_lines, "cellchat_summary.txt")
message("[demo.cellchat_sample] done. outputs:\n  ",
        paste(list.files(), collapse = ", "))