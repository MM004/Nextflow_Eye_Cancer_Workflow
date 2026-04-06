#!/usr/bin/env Rscript
library(ggplot2)
library(optparse)

option_list <- list(
    make_option("--samples", type="character"),
    make_option("--output_dir", type="character", default=".")
)
opts <- parse_args(OptionParser(option_list=option_list))

samples <- read.delim(opts$samples)

genes <- data.frame(
    gene = c("UQCC1", "CRNDE", "ABCC5"),
    file = c("UQCC1_coverage.tsv", "CRNDE_coverage.tsv", "ABCC5_coverage.tsv"),
    stringsAsFactors = FALSE
)

for (i in 1:nrow(genes)) {
    gene <- genes$gene[i]
    cov_file <- genes$file[i]

    if (!file.exists(cov_file) || file.size(cov_file) == 0) {
        cat("Skipping", gene, "- no coverage data\n")
        next
    }

    cov <- read.table(cov_file, header=FALSE, sep="\t",
                    col.names=c("sample_id", "chr", "pos", "depth"))
    cov <- merge(cov, samples[, c("sample_id", "condition")], by="sample_id")

    # Bin to 100bp windows
    cov$bin <- floor(cov$pos / 100) * 100
    agg <- aggregate(depth ~ bin + condition, data=cov, FUN=mean)

    p <- ggplot(agg, aes(x=bin, y=depth, color=condition, fill=condition)) +
        geom_area(alpha=0.3, position="identity") +
        scale_color_manual(values=c(mutant="#E41A1C", wildtype="#377EB8")) +
        scale_fill_manual(values=c(mutant="#E41A1C", wildtype="#377EB8")) +
        labs(title=paste("Read Coverage:", gene),
            subtitle="SF3B1 mutant vs wildtype (mean across samples)",
            x="Genomic Position", y="Mean Read Depth",
            color="SF3B1 Status", fill="SF3B1 Status") +
        theme_minimal() +
        theme(legend.position="bottom")

    outfile <- file.path(opts$output_dir, paste0("coverage_", gene, ".pdf"))
    ggsave(outfile, p, width=10, height=5)
    cat("Saved", outfile, "\n")
}
