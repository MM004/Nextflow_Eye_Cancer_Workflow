nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW      } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED  } from './modules/fastqc'
include { TRIMMOMATIC                } from './modules/trimmomatic'
include { MULTIQC as MULTIQC_RAW    } from './modules/multiqc'
include { MULTIQC as MULTIQC_TRIMMED } from './modules/multiqc'

workflow {
    // Parse sample sheet -> channel of [sample_id, [r1, r2]]
    Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def r1 = file("${params.data_dir}/${row.sra_accession}_1.fastq.gz")
            def r2 = file("${params.data_dir}/${row.sra_accession}_2.fastq.gz")
            tuple(row.sample_id, [r1, r2])
        }
        .set { reads_ch }

    // Phase 2: QC & Trimming
    FASTQC_RAW(reads_ch, "raw")

    TRIMMOMATIC(reads_ch)

    FASTQC_TRIMMED(TRIMMOMATIC.out.trimmed_reads, "trimmed")

    MULTIQC_RAW(FASTQC_RAW.out.zip.collect(), "raw")
    MULTIQC_TRIMMED(FASTQC_TRIMMED.out.zip.collect(), "trimmed")
}
