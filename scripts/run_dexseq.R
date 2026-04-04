#!/usr/bin/env Rscript
library(DEXSeq)
library(optparse)

option_list <- list(
    make_option("--counts_dir", type = "character"),
    make_option("--samples",    type = "character"),
    make_option("--gff",        type = "character"),
    make_option("--fdr",        type = "numeric", default = 0.1),
    make_option("--output",     type = "character"),
    make_option("--output_sig", type = "character")
)
opts <- parse_args(OptionParser(option_list = option_list))

# Read sample metadata
sample_info <- read.delim(opts$samples)
count_files <- file.path(opts$counts_dir, paste0(sample_info$sample_id, "_counts.txt"))
stopifnot(all(file.exists(count_files)))

# Create DEXSeqDataSet
dxd <- DEXSeqDataSetFromHTSeq(
    countfiles    = count_files,
    sampleData    = data.frame(condition = factor(sample_info$condition,
                                                levels = c("wildtype", "mutant"))),
    design        = ~ sample + exon + condition:exon,
    flattenedfile = opts$gff
)

# Run DEXSeq pipeline
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd)

# Extract results
res <- as.data.frame(DEXSeqResults(dxd))

# Remove list columns that write.table can't handle
list_cols <- sapply(res, is.list)
res <- res[, !list_cols]

write.table(res, file = opts$output, sep = "\t", quote = FALSE, row.names = TRUE)

sig <- res[which(res$padj < opts$fdr), ]
sig <- sig[order(sig$padj), ]
write.table(sig, file = opts$output_sig, sep = "\t", quote = FALSE, row.names = TRUE)

cat(sprintf("DEXSeq: %d exonic bins tested, %d significant at FDR < %g\n",
            nrow(res), nrow(sig), opts$fdr))
