## demo.edgeR.R - direct copy of 3.edgeR.R, only the output dir is parametrised
##
## Identical logic to ../code/3.edgeR.R so reviewers can compare the
## output with what the production script would produce. Reads a tab-
## separated count matrix from --count-file, runs edgeR exactTest with
## the same dispersion / cutoffs used in the manuscript, and writes
## both unfiltered and filtered DE tables.
##
## Usage:
##   Rscript demo.edgeR.R --count-file <tsv> --out-prefix <prefix>
##
## Required columns in the input file:
##   row.names = gene id, columns 1..(C_number+T_number) = raw counts
##
## The default values below match what the manuscript pipeline used for
## the Myeloid subset (2 control + 2 PEDV samples).
suppressMessages({
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && length(args) >= i + 1) args[[i + 1]] else default
}
count_file <- get_arg("--count-file", "Myeloid.exp.xls")
outpfx     <- get_arg("--out-prefix", "./Myeloid.MOCK-vs-PEDV")
C_number   <- as.integer(get_arg("--control", "2"))
T_number   <- as.integer(get_arg("--treat",   "2"))
group_colnames <- c("MOCK", "PEDV")
dispersion <- 0.01
cutp       <- 0.05
cutfc      <- log2(1.5)

message("[demo.edgeR] reading ", count_file)
data <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t",
                   check.names = FALSE, quote = "", comment.char = "")
DiffMatrix <- data[, 1:(C_number + T_number)]

filter <- 0.001
DiffMatrix <- DiffMatrix[rowSums(DiffMatrix) > filter, ]

conditions <- factor(c(rep("control", C_number), rep("treat", T_number)))
Diff_list <- DGEList(counts = DiffMatrix, group = conditions)
Diff_list <- calcNormFactors(Diff_list, method = "TMM")
etest <- exactTest(Diff_list, pair = c("control", "treat"),
                   dispersion = dispersion^2)
tTags <- topTags(etest, n = NULL)

data2 <- data
for (i in 1:(C_number + T_number)) {
  data2[, i] <- round(data[, i] * 1000000 / sum(as.numeric(data[, i])), 2)
  colnames(data2)[i] <- paste(colnames(data)[i], "_CPM", sep = "")
}
data3 <- data2[, 1:2]
data3[, 1] <- apply(data.frame(data2[, 1:C_number]), 1, mean)
data3[, 2] <- apply(data.frame(data2[, (C_number + 1):(C_number + T_number)]),
                    1, mean)
log2fc <- log2((data3[, 2] + 1e-10) / (data3[, 1] + 1e-10))
colnames(data3) <- paste(group_colnames, "_mean", sep = "")

data4 <- data.frame(id = rownames(data), data, data2, data3,
                    log2fc = log2fc, check.names = FALSE)
datat <- tTags@.Data[[1]][, c("logFC", "logCPM", "PValue", "FDR")]
datat$id <- rownames(datat)
data5 <- merge(data4, datat, by = "id")

write.table(data5, file = paste0(outpfx, ".all.xls"),
            sep = "\t", quote = FALSE, row.names = FALSE)
data6 <- data5[which(abs(data5$log2fc) > cutfc & data5$PValue < cutp), ]
write.table(data6, file = paste0(outpfx, ".filter.xls"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_sig <- nrow(data6)
top_up   <- head(data6[order(-data6$log2fc), c("id", "log2fc", "PValue")], 5)
top_down <- head(data6[order( data6$log2fc), c("id", "log2fc", "PValue")], 5)

summary_lines <- c(
  sprintf("input_genes           : %d", nrow(data)),
  sprintf("input_genes_after_flt : %d", nrow(DiffMatrix)),
  sprintf("DE_genes_sig          : %d", n_sig),
  "top_up_regulated:",
  sprintf("  %s log2FC=%.2f P=%.3g",
          top_up$id,   top_up$log2fc,   top_up$PValue),
  "top_down_regulated:",
  sprintf("  %s log2FC=%.2f P=%.3g",
          top_down$id, top_down$log2fc, top_down$PValue)
)
writeLines(summary_lines, paste0(outpfx, ".summary.txt"))
message("[demo.edgeR] wrote ", outpfx, ".all.xls / .filter.xls / .summary.txt")