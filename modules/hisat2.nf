process HISAT2_EXTRACT_SPLICESITES {                                                                                                                                                                        
      publishDir "${params.outdir}/references", mode: 'copy'                        
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
      publishDir "${params.outdir}/aligned", mode: 'copy', pattern: '*.log'

      input:
      tuple val(sample_id), path(reads) // ("sample_1", [sample_1_R1.trimmed.fastq.gz, sample_1_R2.trimmed.fastq.gz])     
      path index			// all the HISAT2 index files (genome.1.ht2, genome.2.ht2, ... genome.8.ht2) 
      path splicesites			// splicesites.tsv from the extract step

      output:
      tuple val(sample_id), path("*.sam"), emit: sam
      path "*.log"                       , emit: log

      script:
      def index_base = index[0].toString().replaceAll(/\.\d\.ht2$/, '')
      """
      hisat2 \\
          -x ${index_base} \\
          -1 ${reads[0]} \\
          -2 ${reads[1]} \\
          --known-splicesite-infile ${splicesites} \\
          --dta \\
          -p ${task.cpus} \\
          --summary-file ${sample_id}_hisat2.log \\
          -S ${sample_id}.sam
      """
  }
