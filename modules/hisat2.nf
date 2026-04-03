process HISAT2_EXTRACT_SPLICESITES {                                                                                                                                                                        
      publishDir "${params.outdir}/references", mode: 'symlink'                        

      input:                                                                                                                                                                                                
      path gtf

      output:
      path "splicesites.tsv", emit: splicesites

      script:
      """
      hisat2_extract_splice_sites.py ${gtf} > splicesites.tsv
      """
  }

process HISAT2_ALIGN {
      tag "${sample_id}"
      publishDir "${params.outdir}/aligned", mode: 'symlink', pattern: '*.log'
      publishDir "${params.outdir}/dedup", mode: 'symlink', pattern: '*.{bam,bai,txt}'

      input:
      tuple val(sample_id), path(reads) // ("sample_1", [sample_1_R1.trimmed.fastq.gz, sample_1_R2.trimmed.fastq.gz])     
      path index			// all the HISAT2 index files (genome.1.ht2, genome.2.ht2, ... genome.8.ht2) 
      path splicesites			// splicesites.tsv from the extract step

      output:
      tuple val(sample_id), path("*.dedup.bam"), path("*.dedup.bam.bai"), emit: bam
      path "*.metrics.txt", emit: metrics
      path "*.log"        , emit: log

      script:
      def index_base = index[0].toString().replaceAll(/\.\d\.ht2$/, '')
      """
      # Align and sort
      hisat2 \\
          -x ${index_base} \\
          -1 ${reads[0]} \\
          -2 ${reads[1]} \\
          --known-splicesite-infile ${splicesites} \\
          --dta \\
          --rg-id ${sample_id} \\
          --rg SM:${sample_id} \\
          --rg PL:ILLUMINA \\
          --rg LB:${sample_id} \\
          -p ${task.cpus} \\
          --summary-file ${sample_id}_hisat2.log \\
      | samtools sort -@ 4 -o ${sample_id}.sorted.bam

      samtools index ${sample_id}.sorted.bam

      # Mark duplicates
      picard MarkDuplicates \\
          I=${sample_id}.sorted.bam \\
          O=${sample_id}.dedup.bam \\
          M=${sample_id}.metrics.txt \\
          REMOVE_DUPLICATES=false \\
          CREATE_INDEX=true \\
          VALIDATION_STRINGENCY=LENIENT

      mv ${sample_id}.dedup.bai ${sample_id}.dedup.bam.bai

      # Delete intermediate to free space
      rm ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai
      """
  }
