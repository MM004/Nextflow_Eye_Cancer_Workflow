#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
    make_option("--dexseq", type="character"),
    make_option("--rmats_dir", type="character"),
    make_option("--fdr", type="numeric", default=0.1),
    make_option("--output", type="character")
)
opts <- parse_args(OptionParser(option_list=option_list))

targets <- data.frame(
    ensembl_id = c("ENSG00000114770", "ENSG00000245694", "ENSG00000101019",
                    "ENSG00000228315", "ENSG00000131503", "ENSG00000148543"),
    gene_name  = c("ABCC5", "CRNDE", "UQCC1", "GUSBP11", "ANKHD1", "ADAM12"),
    stringsAsFactors = FALSE
)

# Parse DEXSeq significant results
dexseq <- read.table(opts$dexseq, header=TRUE, row.names=1, sep="\t",
                    stringsAsFactors=FALSE, check.names=FALSE)

dex_results <- lapply(targets$ensembl_id, function(gid) {
    hits <- grepl(gid, rownames(dexseq), fixed=TRUE)
    if (any(hits)) list(count=sum(hits), min_fdr=min(dexseq$padj[hits], na.rm=TRUE))
    else list(count=0, min_fdr=NA)
})

# Parse rMATS results
event_types <- c("SE", "A3SS", "A5SS", "MXE", "RI")
rmats_all <- do.call(rbind, lapply(event_types, function(et) {
    f <- file.path(opts$rmats_dir, paste0(et, ".MATS.JC.txt"))
    if (!file.exists(f)) return(NULL)
    d <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    d$EventType <- et
    d$geneSymbol <- gsub('"', '', d$geneSymbol)
    d$GeneID <- gsub('"', '', d$GeneID)
    d[, c("GeneID", "geneSymbol", "FDR", "EventType")]
}))
rmats_sig <- rmats_all[rmats_all$FDR < opts$fdr, ]

rma_results <- lapply(targets$gene_name, function(gname) {
    hits <- grepl(gname, rmats_sig$geneSymbol)
    if (any(hits)) list(count=sum(hits),
                        types=paste(unique(rmats_sig$EventType[hits]), collapse=","),
                        min_fdr=min(rmats_sig$FDR[hits]))
    else list(count=0, types=NA, min_fdr=NA)
})

# Build summary
summary_df <- data.frame(
    Gene           = targets$gene_name,
    Ensembl_ID     = targets$ensembl_id,
    DEXSeq_bins    = sapply(dex_results, `[[`, "count"),
    DEXSeq_min_FDR = sapply(dex_results, function(x) ifelse(is.na(x$min_fdr), "NA",
formatC(x$min_fdr, format="e", digits=2))),
    rMATS_events   = sapply(rma_results, `[[`, "count"),
    rMATS_types    = sapply(rma_results, function(x) ifelse(is.na(x$types), "NA", x$types)),
    rMATS_min_FDR  = sapply(rma_results, function(x) ifelse(is.na(x$min_fdr), "NA",
formatC(x$min_fdr, format="e", digits=2))),
    Reproduced     = ifelse(sapply(dex_results, `[[`, "count") > 0 | sapply(rma_results, `[[`,
"count") > 0, "YES", "no"),
    stringsAsFactors = FALSE
)

write.table(summary_df, file=opts$output, sep="\t", quote=FALSE, row.names=FALSE)
cat("Summary:\n")
print(summary_df)
