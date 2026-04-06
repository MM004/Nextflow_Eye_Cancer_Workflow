process DEXSEQ_PREPARE_ANNOTATION {
    publishDir "${params.outdir}/dexseq", mode: 'symlink'

    input:
    path gtf

    output:
    path "dexseq.gff", emit: gff

    script:
    """
    SCRIPTS=\$(conda run -n splicing-env Rscript -e "cat(system.file('python_scripts', package='DEXSeq'))")
    conda run -n splicing-env python \$SCRIPTS/dexseq_prepare_annotation.py ${gtf} dexseq.gff
    """
}

process DEXSEQ_COUNT {
    tag "${sample_id}"
    publishDir "${params.outdir}/dexseq/counts", mode: 'symlink'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path gff

    output:
    path "${sample_id}_counts.txt", emit: counts

    script:
    """
    SCRIPTS=\$(conda run -n splicing-env Rscript -e "cat(system.file('python_scripts', package='DEXSeq'))")
    conda run -n splicing-env python \$SCRIPTS/dexseq_count.py \
            -p yes \
            -r pos \
            -s no \
            -f bam \
            ${gff} ${bam} ${sample_id}_counts.txt
    """
}

process DEXSEQ_ANALYSIS {
    publishDir "${params.outdir}/dexseq", mode: 'symlink'

    input:
    path counts
    path samples
    path gff

    output:
    path "dexseq_results.tsv", emit: results
    path "dexseq_results_significant.tsv", emit: significant

    script:
    """
    for f in *_counts.txt; do
        grep -v '^_' "\$f" | sed 's/"//g' > "\${f}.tmp" && mv "\${f}.tmp" "\$f"
    done
    conda run -n splicing-env Rscript ${projectDir}/scripts/run_dexseq.R \
            --counts_dir . \
            --samples ${samples} \
            --gff ${gff} \
            --fdr ${params.fdr_threshold} \
            --output dexseq_results.tsv \
            --output_sig dexseq_results_significant.tsv
    """
}
