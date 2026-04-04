process RMATS {
    publishDir "${params.outdir}/rmats", mode: 'symlink'

    input:
    path bams
    path samples
    path gtf

    output:
    path "output/*.MATS.JC.txt", emit: jc_results
    path "output/*.MATS.JCEC.txt", emit: jcec_results

    script:
    """
    awk -F'\\t' 'NR>1 && \$2=="mutant"  {printf "%s.dedup.bam,", \$1}' ${samples} | sed 's/,\$//' > b1.txt
    awk -F'\\t' 'NR>1 && \$2=="wildtype" {printf "%s.dedup.bam,", \$1}' ${samples} | sed 's/,\$//' > b2.txt

    conda run -n splicing-env rmats.py \
        --b1 b1.txt \
        --b2 b2.txt \
        --gtf ${gtf} \
        -t paired \
        --readLength 99 \
        --variable-read-length \
        --nthread ${task.cpus} \
        --od output \
        --tmp tmp
    """
}
