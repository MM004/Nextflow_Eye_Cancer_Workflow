process SAM_TO_SORTED_BAM {
      tag "${sample_id}"
      publishDir "${params.outdir}/aligned", mode: 'copy'

      input:
      tuple val(sample_id), path(sam)

      output:
      tuple val(sample_id), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam

      script:
      """
      samtools view -@ ${task.cpus} -bS ${sam} | \\
          samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
      samtools index ${sample_id}.sorted.bam
      """
  }
