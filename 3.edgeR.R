library(edgeR)

# Parameters
count_file <- "Myeloid.exp.xls"
C_number <- 2
T_number <- 2
group_colnames <- c("MOCK", "PEDV")
outpfx <- "./Myeloid.MOCK-vs-PEDV"
dispersion <- 0.01
cutp <- 0.05
cutfc <- log2(1.5)
fpkm_file <- "Myeloid.exp.xls"

# Read data
data <- read.table(count_file, header=T, row.names=1, sep='\t', check.names=F, quote="", comment.char="")
DiffMatrix <- data[, 1:(C_number+T_number)]

# Filter low-expressed genes
filter <- 0.001
DiffMatrix <- DiffMatrix[rowSums(DiffMatrix) > filter, ]

# DEG analysis
conditions <- factor(c(rep("control", C_number), rep("treat", T_number)))
Diff_list <- DGEList(counts=DiffMatrix, group=conditions)
Diff_list <- calcNormFactors(Diff_list, method="TMM")
etest <- exactTest(Diff_list, pair=c("control", "treat"), dispersion=dispersion^2)
tTags <- topTags(etest, n=NULL)

# Calculate CPM and mean
data2 <- data
for(i in 1:(C_number+T_number)) {
    data2[, i] <- round(data[, i] * 1000000 / sum(as.numeric(data[, i])), 2)
    colnames(data2)[i] <- paste(colnames(data)[i], "_CPM", sep="")
}

data3 <- data2[, 1:2]
data3[, 1] <- apply(data.frame(data2[, 1:C_number]), 1, mean)
data3[, 2] <- apply(data.frame(data2[, (C_number+1):(C_number+T_number)]), 1, mean)
log2fc <- log2((data3[, 2] + 1e-10) / (data3[, 1] + 1e-10))
colnames(data3) <- paste(group_colnames, "_mean", sep="")

# Merge and output
data4 <- data.frame(id=rownames(data), data, data2, data3, log2fc=log2fc, check.names=FALSE)
datat <- tTags@.Data[[1]][, c("logFC", "logCPM", "PValue", "FDR")]
datat$id <- rownames(datat)
data5 <- merge(data4, datat, by="id")

# Output
write.table(data5, file=paste0(outpfx, ".all.xls"), sep='\t', quote=F, row.names=F)
data6 <- data5[which(abs(data5$log2fc) > cutfc & data5$PValue < cutp), ]
write.table(data6, file=paste0(outpfx, ".filter.xls"), sep='\t', quote=F, row.names=F)
