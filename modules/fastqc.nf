process FASTQC {
      tag "${sample_id}"
      publishDir "${params.outdir}/fastqc_${stage}", mode: 'copy'

      input:
      tuple val(sample_id), path(reads)
      val stage

      output:
      path "*.html", emit: html
      path "*.zip" , emit: zip

      script:
      """
      fastqc --threads ${task.cpus} ${reads}
      """
}
