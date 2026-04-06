process ALIGNMENT_QC {
      tag "${sample_id}"
      publishDir "${params.outdir}/alignment_qc", mode: 'symlink'

      input:
      tuple val(sample_id), path(bam), path(bai)

      output:
      path "*.flagstat", emit: flagstat
      path "*.idxstats", emit: idxstats

      script:
      """
      samtools flagstat ${bam} > ${sample_id}.flagstat
      samtools idxstats ${bam} > ${sample_id}.idxstats
      """
  }
