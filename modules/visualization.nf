process SPLICING_SUMMARY {
    publishDir "${params.outdir}/summary", mode: 'symlink'

    input:
    path dexseq_sig
    path rmats_results
    path samples

    output:
    path "splicing_summary.tsv", emit: summary

    script:
    """
    conda run -n splicing-env Rscript ${projectDir}/scripts/summarize_splicing.R --dexseq ${dexseq_sig} --rmats_dir . --fdr ${params.fdr_threshold} --output splicing_summary.tsv
    """
}

process COVERAGE_PLOTS {
    publishDir "${params.outdir}/figures", mode: 'symlink'

    input:
    path bams
    path samples

    output:
    path "coverage_*.pdf", emit: plots

    script:
    """
    for bam in *.dedup.bam; do
        sample=\$(basename \$bam .dedup.bam)
        samtools depth -a -r 20:33890000-33936000 \$bam | awk -v s=\$sample '{print s"\\t"\$0}' >> UQCC1_coverage.tsv
        samtools depth -a -r 16:54952000-54958000 \$bam | awk -v s=\$sample '{print s"\\t"\$0}' >> CRNDE_coverage.tsv
        samtools depth -a -r 3:183698000-183708000 \$bam | awk -v s=\$sample '{print s"\\t"\$0}' >> ABCC5_coverage.tsv
    done
    conda run -n splicing-env Rscript ${projectDir}/scripts/plot_coverage.R --samples ${samples} --output_dir .
    """
}

process SPLICING_HEATMAP {
    publishDir "${params.outdir}/figures", mode: 'symlink'

    input:
    path rmats_results
    path samples

    output:
    path "splicing_heatmap.pdf", emit: heatmap

    script:
    """
    conda run -n splicing-env Rscript ${projectDir}/scripts/plot_heatmap.R --rmats_dir . --samples ${samples} --fdr ${params.fdr_threshold} --output splicing_heatmap.pdf
    """
}
