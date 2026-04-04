process MULTIQC {
    publishDir "${params.outdir}/multiqc_${stage}", mode: 'symlink'

    input:
    path '*'
    val stage

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data

    script:
    """
    multiqc .
    """
}
