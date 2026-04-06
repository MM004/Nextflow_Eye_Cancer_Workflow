#!/usr/bin/env Rscript
library(pheatmap)
library(optparse)

option_list <- list(
    make_option("--rmats_dir", type="character"),
    make_option("--samples", type="character"),
    make_option("--fdr", type="numeric", default=0.1),
    make_option("--output", type="character")
)
opts <- parse_args(OptionParser(option_list=option_list))

samples <- read.delim(opts$samples)
target_genes <- c("ABCC5", "CRNDE", "UQCC1", "GUSBP11", "ANKHD1", "ADAM12")

# Read all rMATS JC results
all_df <- do.call(rbind, lapply(c("SE", "A3SS", "A5SS", "RI"), function(et) {
    f <- file.path(opts$rmats_dir, paste0(et, ".MATS.JC.txt"))
    if (!file.exists(f)) return(NULL)
    d <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    d$EventType <- et
    d$geneSymbol <- gsub('"', '', d$geneSymbol)
    d[, c("geneSymbol", "FDR","IncLevel1", "IncLevel2", "EventType")]
}))

# Build PSI matrix: genes x samples
psi_matrix <- matrix(NA, nrow=length(target_genes), ncol=nrow(samples))
rownames(psi_matrix) <- target_genes
colnames(psi_matrix) <- samples$sample_id

mutant_ids <- samples$sample_id[samples$condition == "mutant"]
wildtype_ids <- samples$sample_id[samples$condition == "wildtype"]

for (gene in target_genes) {
    hits <- all_df[grepl(gene, all_df$geneSymbol) & all_df$FDR < opts$fdr, ]
    if (nrow(hits) == 0) next

    best <- hits[which.min(hits$FDR), ]
    psi1 <- as.numeric(strsplit(best$IncLevel1, ",")[[1]])
    psi2 <- as.numeric(strsplit(best$IncLevel2, ",")[[1]])

    if (length(psi1) == length(mutant_ids)) psi_matrix[gene, mutant_ids] <- psi1
    if (length(psi2) == length(wildtype_ids)) psi_matrix[gene, wildtype_ids] <- psi2
}

# Remove genes with no data
psi_matrix <- psi_matrix[rowSums(!is.na(psi_matrix)) > 0, , drop=FALSE]

annotation_col <- data.frame(SF3B1=samples$condition, row.names=samples$sample_id)
annotation_colors <- list(SF3B1=c(mutant="#E41A1C", wildtype="#377EB8"))

pdf(opts$output, width=8, height=5)
pheatmap(psi_matrix,
        color=colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        cluster_rows=FALSE, cluster_cols=FALSE,
        annotation_col=annotation_col,
        annotation_colors=annotation_colors,
        main="Splicing Inclusion Levels (PSI) — Target Genes",
        display_numbers=TRUE, number_format="%.2f", fontsize_number=8,
        na_col="grey90")
dev.off()
cat("Heatmap saved to", opts$output, "\n")
