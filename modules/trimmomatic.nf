process TRIMMOMATIC {
      tag "${sample_id}"
      publishDir "${params.outdir}/trimmed", mode: 'symlink', pattern: '*.trimmed.fastq.gz'
      publishDir "${params.outdir}/trimmed/logs", mode: 'symlink', pattern: '*.log'

      input:
      tuple val(sample_id), path(reads)

      output:
      tuple val(sample_id), path("*R{1,2}.trimmed.fastq.gz"), emit: trimmed_reads
      path "*.unpaired.fastq.gz", emit: unpaired
      path "*.log"              , emit: log

      script:
      """
      trimmomatic PE \\
          -threads ${task.cpus} \\
          ${reads[0]} ${reads[1]} \\
          ${sample_id}_R1.trimmed.fastq.gz ${sample_id}_R1.unpaired.fastq.gz \\
          ${sample_id}_R2.trimmed.fastq.gz ${sample_id}_R2.unpaired.fastq.gz \\
          ILLUMINACLIP:${params.adapters}:2:30:10 \\
          CROP:${params.trimmed_read_length} \\
          LEADING:3 \\
          TRAILING:3 \\
          SLIDINGWINDOW:4:15 \\
          MINLEN:36 \\
          2> ${sample_id}_trimmomatic.log
      """
}
