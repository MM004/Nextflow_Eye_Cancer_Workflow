#!/usr/bin/env Rscript
library(ggplot2)
library(optparse)

option_list <- list(
    make_option("--samples", type="character"),
    make_option("--output_dir", type="character", default=".")
)
opts <- parse_args(OptionParser(option_list=option_list))

samples <- read.delim(opts$samples)

# --- RPM normalization: get total mapped reads per sample from BAM index ---
library(Rsamtools)

bam_files <- list.files(".", pattern="\\.dedup\\.bam$", full.names=TRUE)
lib_sizes <- data.frame(
    sample_id = sub("\\.dedup\\.bam$", "", basename(bam_files)),
    total_reads = sapply(bam_files, function(f) countBam(f)$records),
    stringsAsFactors = FALSE
)
lib_sizes$rpm_factor <- lib_sizes$total_reads / 1e6
cat("Library sizes:\n")
print(lib_sizes[, c("sample_id", "total_reads")])

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

    # Normalize depth to RPM
    cov <- merge(cov, lib_sizes[, c("sample_id", "rpm_factor")], by="sample_id")
    cov$depth_rpm <- cov$depth / cov$rpm_factor

    # Bin to 100bp windows
    cov$bin <- floor(cov$pos / 100) * 100
    agg <- aggregate(depth_rpm ~ bin + condition, data=cov, FUN=mean)

    p <- ggplot(agg, aes(x=bin, y=depth_rpm, color=condition, fill=condition)) +
        geom_area(alpha=0.3, position="identity") +
        scale_color_manual(values=c(mutant="#E41A1C", wildtype="#377EB8")) +
        scale_fill_manual(values=c(mutant="#E41A1C", wildtype="#377EB8")) +
        labs(title=paste("Read Coverage:", gene),
            subtitle="SF3B1 mutant vs wildtype (RPM-normalized, mean across samples)",
            x="Genomic Position", y="Read Depth (RPM)",
            color="SF3B1 Status", fill="SF3B1 Status") +
        theme_minimal() +
        theme(legend.position="bottom")

    outfile <- file.path(opts$output_dir, paste0("coverage_", gene, ".pdf"))
    ggsave(outfile, p, width=10, height=5)
    cat("Saved", outfile, "\n")
}
