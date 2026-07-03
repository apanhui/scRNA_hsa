library(CellChat)

# Parameters
infile <- "cellchat.init.Rds"
outdir <- "Control"

# Initialize
cellchat <- readRDS(infile)
dir.create(outdir, F, T)
setwd(outdir)

# Cell-cell communication analysis
cellchat <- identifyOverExpressedInteractions(cellchat, features = rownames(cellchat@data.signaling))
cellchat@idents <- droplevels(cellchat@idents)
population.size <- if (do.call(`/`, as.list(rev(range(table(cellchat@idents))))) > 20) TRUE else FALSE

cellchat <- computeCommunProb(cellchat, population.size = population.size)
cellchat <- SelectCommunication(cellchat,
    source = cellchat@options$select_commun$source,
    target = cellchat@options$select_commun$target,
    pairs = cellchat@options$select_commun$pairs, drop = cellchat@options$select_commun$drop,
    method = cellchat@options$select_commun$method)
cellchat <- aggregateNet(cellchat, thresh = 0.05)

# Output LR results
data <- GetOnlineData(cellchat)
WriteTable(data, file = 'CellChat.LR.all.xls')
data <- subset(data, Pvalue <= 0.05)
write.table(data, file = 'CellChat.LR.filtered.xls', row.names = F, sep = "\t", quote = F)

# Pathway analysis
cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)
data <- StatCommunProb(cellchat, outpref = "Pathway.CommunProb", slot.name = "netP", thresh = 0.05)

# Top5 pathway visualization
pathways.show <- (subset(data, prob > 0) %>% group_by(pathway_name) %>% summarise(pval = sum(pval), prob = sum(prob)) %>% arrange(pval, -prob) %>% head(5))[['pathway_name']]
writeLines(pathways.show, 'pathways.show.list')

for (i in pathways.show) {
    pdf(paste0("Pathway.", i, ".hierarchy.pdf"), 6, 6)
    netVisual_aggregate(cellchat, signaling = i, layout = 'hierarchy', vertex.receiver = 1:round(nlevels(cellchat@idents) / 2))
    dev.off()

    pdf(paste0("Pathway.", i, ".networks.pdf"), 6, 6)
    netVisual_aggregate(cellchat, signaling = i, layout = "circle")
    dev.off()

    pdf(paste0("Pathway.", i, ".heatmap.pdf"), 6, 6)
    ht <- netVisual_heatmap(cellchat, signaling = i, color.heatmap = c("white", "darkred"))
    print(ht)
    dev.off()
}

# Save results
saveRDS(cellchat, file = 'cellchat.Rds')
