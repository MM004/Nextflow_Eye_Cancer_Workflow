process MARK_DUPLICATES {
      tag "${sample_id}"
      publishDir "${params.outdir}/dedup", mode: 'copy'

      input:
      tuple val(sample_id), path(bam), path(bai)

      output:
      tuple val(sample_id), path("*.dedup.bam"), path("*.dedup.bam.bai"), emit: bam
      path "*.metrics.txt", emit: metrics

      script:
      """
      picard MarkDuplicates \\
          I=${bam} \\
          O=${sample_id}.dedup.bam \\
          M=${sample_id}.metrics.txt \\
          REMOVE_DUPLICATES=false \\
          CREATE_INDEX=true \\
          VALIDATION_STRINGENCY=LENIENT

      mv ${sample_id}.dedup.bai ${sample_id}.dedup.bam.bai
      """
  }
