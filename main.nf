nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW      } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED  } from './modules/fastqc'
include { TRIMMOMATIC                } from './modules/trimmomatic'
include { MULTIQC as MULTIQC_RAW    } from './modules/multiqc'
include { MULTIQC as MULTIQC_TRIMMED } from './modules/multiqc'
include { HISAT2_EXTRACT_SPLICESITES } from './modules/hisat2'
include { HISAT2_ALIGN               } from './modules/hisat2'
include { ALIGNMENT_QC               } from './modules/alignment_qc'
include { DEXSEQ_PREPARE_ANNOTATION  } from './modules/dexseq'
include { DEXSEQ_COUNT               } from './modules/dexseq'
include { DEXSEQ_ANALYSIS            } from './modules/dexseq'
include { RMATS                      } from './modules/rmats'
include { SPLICING_SUMMARY  } from './modules/visualization'
include { COVERAGE_PLOTS    } from './modules/visualization'
include { SPLICING_HEATMAP  } from './modules/visualization'
include { MULTIQC as MULTIQC_FINAL } from './modules/multiqc'

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
    
    // Reference files
    hisat2_index = Channel.fromPath("${params.hisat2_index}*.ht2").collect()
    gtf          = file(params.gtf)
    samples_file = file(params.samples)

    // Phase 2: QC & Trimming
    FASTQC_RAW(reads_ch, "raw")

    TRIMMOMATIC(reads_ch)

    FASTQC_TRIMMED(TRIMMOMATIC.out.trimmed_reads, "trimmed")

    MULTIQC_RAW(FASTQC_RAW.out.zip.collect(), "raw")
   
    MULTIQC_TRIMMED(FASTQC_TRIMMED.out.zip.collect(), "trimmed")

    // Phase 3: Alignment
    HISAT2_EXTRACT_SPLICESITES(gtf)
    HISAT2_ALIGN(
        TRIMMOMATIC.out.trimmed_reads,
        hisat2_index,
        HISAT2_EXTRACT_SPLICESITES.out.splicesites.first()
    )
    ALIGNMENT_QC(HISAT2_ALIGN.out.bam)
    
    // Phase 4A: Seq (differential exon usage)
    DEXSEQ_PREPARE_ANNOTATION(gtf)
    DEXSEQ_COUNT(
        HISAT2_ALIGN.out.bam,
        DEXSEQ_PREPARE_ANNOTATION.out.gff.first()
    )
    DEXSEQ_ANALYSIS(
        DEXSEQ_COUNT.out.counts.collect(),
        samples_file,
        DEXSEQ_PREPARE_ANNOTATION.out.gff.first()
    )

    // Phase 4B: rMATS (alternative splicing events)
    RMATS(
        HISAT2_ALIGN.out.bam
            .flatMap { sample_id, bam, bai -> [bam, bai] }
            .collect(),
        samples_file,
        gtf
    )

    // Phase 5: Visualization & Summary
    SPLICING_SUMMARY(
        DEXSEQ_ANALYSIS.out.significant,
        RMATS.out.jc_results.collect(),
        samples_file
    )

    COVERAGE_PLOTS(
        HISAT2_ALIGN.out.bam
            .flatMap { sample_id, bam, bai -> [bam, bai] }
            .collect(),
        samples_file
    )

    SPLICING_HEATMAP(
        RMATS.out.jc_results.collect(),
        samples_file
    )

    MULTIQC_FINAL(
        FASTQC_RAW.out.zip
            .mix(FASTQC_TRIMMED.out.zip)
            .mix(TRIMMOMATIC.out.log)
            .mix(HISAT2_ALIGN.out.log)
            .mix(HISAT2_ALIGN.out.metrics)
            .mix(ALIGNMENT_QC.out.flagstat)
            .collect(),
        "final"
    )
}
