# Nextflow Eye-Cancer Workflow

Reproducing the RNA-seq analysis from Furney et al. (2013), which showed that SF3B1 mutations are
associated with differential alternative splicing in uveal melanoma.

## References

1. **Furney et al.** "SF3B1 mutations are associated with alternative splicing in uveal melanoma" —
*Cancer Discov.* 2013
   https://pubmed.ncbi.nlm.nih.gov/23861464/
2. **Harbour et al.** "Recurrent mutations at codon 625 of the splicing factor SF3B1 in uveal
melanoma" — *Nat Genet.* 2013
   https://pubmed.ncbi.nlm.nih.gov/23313955/

## Samples

8 uveal melanoma class 1 tumors from Harbour et al. (SRA062359), re-analyzed by Furney et al.:

| Sample ID | Tumor   | SRR Accession | SF3B1 Status | Mutation |
|-----------|---------|---------------|--------------|----------|
| sample_1  | MM010   | SRR628582     | mutant       | R625C    |
| sample_2  | MM065T  | SRR628583     | mutant       | R625C    |
| sample_3  | MM064T  | SRR628584     | mutant       | R625H    |
| sample_4  | MM176T  | SRR628589     | mutant       | R625C    |
| sample_5  | MM016T  | SRR628585     | wildtype     | —        |
| sample_6  | MM082T  | SRR628586     | wildtype     | —        |
| sample_7  | MM089T  | SRR628587     | wildtype     | —        |
| sample_8  | MM132T  | SRR628588     | wildtype     | —        |

## Pipeline Overview

raw FASTQ ── FastQC ──────────────────────────── MultiQC (raw)
   │
   └── Trimmomatic (99bp crop) ── FastQC ─────── MultiQC (trimmed)
         │
         └── HISAT2 (GRCh37) ── samtools sort ── Picard MarkDuplicates
                                                      │
                                                      ├── samtools flagstat/idxstats
                                                      ├── DEXSeq (differential exon usage)
                                                      ├── rMATS (alternative splicing events)
                                                      └── Visualization (coverage plots, heatmap)

## Setup

### 1. Clone the repository
```bash
git clone <repo-url>
cd Nextflow_Eye_Cancer_Workflow

### 2. Create conda environments

# Main pipeline environment (QC, trimming, alignment)
conda env create -f environment.yml
conda activate uveal-pipeline

# Splicing analysis environment (DEXSeq, rMATS — requires Python 3.12)
conda env create -f splicing-env.yml

### 3. Download references

mkdir -p references && cd references

# HISAT2 pre-built index for GRCh37
wget https://genome-idx.s3.amazonaws.com/hisat/grch37_genome.tar.gz
tar -xzf grch37_genome.tar.gz && rm grch37_genome.tar.gz

# Genome FASTA and GTF (Ensembl release 75 — last GRCh37 build)
wget https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.prima
ry_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip *.gz

cd ..

### 4. Symlink or download FASTQ data

ln -s /path/to/fastq/files data
Expected files: data/SRR6285{82-89}_{1,2}.fastq.gz

### 5. Run the pipeline

nextflow run main.nf -resume
Phase 4 processes (DEXSeq, rMATS) automatically call the splicing-env environment via conda run.

**Tools and Justifications**
┌────────────────────┬──────────────────┬────────────┬────────────────────────────────────────────┐
│        Step        │ Original Paper   │ Our Choice │               Justification                │
│                    │     (Furney)     │            │                                            │
├────────────────────┼──────────────────┼────────────┼────────────────────────────────────────────┤
│ Alignment          │ TopHat           │ HISAT2     │ TopHat deprecated; HISAT2 is its successor │
│                    │                  │            │  by the same authors                       │
├────────────────────┼──────────────────┼────────────┼────────────────────────────────────────────┤
│ Duplicate marking  │ Picard           │ Picard     │ Same tool                                  │
├────────────────────┼──────────────────┼────────────┼────────────────────────────────────────────┤
│ Differential exon  │ DEXSeq           │ DEXSeq     │ Same tool, still maintained                │
│ usage              │                  │            │                                            │
├────────────────────┼──────────────────┼────────────┼────────────────────────────────────────────┤
│ Alternative        │ MATS             │ rMATS      │ Updated version of same algorithm          │
│ splicing           │                  │ (turbo)    │                                            │
├────────────────────┼──────────────────┼────────────┼────────────────────────────────────────────┤
│ Read trimming      │ 99bp crop        │ 99bp crop  │ Matches Furney methodology                 │
├────────────────────┼──────────────────┼────────────┼────────────────────────────────────────────┤
│ Reference genome   │ GRCh37           │ GRCh37     │ Same genome for reproducibility            │
├────────────────────┼──────────────────┼────────────┼────────────────────────────────────────────┤
│ Pipeline manager   │ —                │ Nextflow   │ Course requirement; ensures                │
│                    │                  │ DSL2       │ reproducibility                            │
└────────────────────┴──────────────────┴────────────┴────────────────────────────────────────────┘
**Phase 2 Results — QC Summary**

- Read survival after trimming: ~96.7% across all samples
- Reads cropped to 99bp as per Furney methodology
- Outlier: sample_4 (MM176T) — 60% duplication rate, 52% GC content, fails FastQC GC content check.
Kept in analysis (used in original paper).

**Phase 3 Results — Alignment Summary**
┌────────────────────────┬───────┬───────────────┬───────────┬──────────────┐
│         Sample         │ Reads │ Concordant 1x │ Multi-map │ Overall Rate │
├────────────────────────┼───────┼───────────────┼───────────┼──────────────┤
│ sample_1 (MM010, mut)  │ 13.8M │     92.4%     │   5.0%    │    98.70%    │
├────────────────────────┼───────┼───────────────┼───────────┼──────────────┤
│ sample_2 (MM065T, mut) │ 16.8M │     82.6%     │   14.5%   │    98.57%    │
├────────────────────────┼───────┼───────────────┼───────────┼──────────────┤
│ sample_3 (MM064T, mut) │ 10.3M │     92.4%     │   4.3%    │    98.39%    │
├────────────────────────┼───────┼───────────────┼───────────┼──────────────┤
│ sample_4 (MM176T, mut) │ 17.5M │     85.6%     │   8.1%    │    96.53%    │
├────────────────────────┼───────┼───────────────┼───────────┼──────────────┤
│ sample_5 (MM016T, wt)  │ 14.0M │     82.9%     │   14.3%   │    98.60%    │
├────────────────────────┼───────┼───────────────┼───────────┼──────────────┤
│ sample_6 (MM082T, wt)  │ 10.4M │     91.7%     │   5.2%    │    98.47%    │
├────────────────────────┼───────┼───────────────┼───────────┼──────────────┤
│ sample_7 (MM089T, wt)  │ 16.6M │     84.1%     │   13.5%   │    98.79%    │
├────────────────────────┼───────┼───────────────┼───────────┼──────────────┤
│ sample_8 (MM132T, wt)  │ 12.7M │     89.6%     │   7.7%    │    98.67%    │
└────────────────────────┴───────┴───────────────┴───────────┴──────────────┘
All samples >96.5% alignment rate. sample_4 (MM176T) is the weakest, consistent with QC outlier
status.

**Phase 4 — Differential Splicing Results**

DEXSeq: 644,354 exonic bins tested, 222 significant at FDR < 0.1

rMATS: 247 A3SS, 76 A5SS, 345 RI, 654 SE significant events at FDR < 0.1

Reproduction of Furney et al. Table 2 target genes:
┌─────────┬───────────────────────┬────────────────────┬─────────────┐
│  Gene   │        DEXSeq         │       rMATS        │ Reproduced? │
├─────────┼───────────────────────┼────────────────────┼─────────────┤
│ UQCC    │ 12 bins (FDR 2.9e-08) │         —          │     YES     │
├─────────┼───────────────────────┼────────────────────┼─────────────┤
│ CRNDE   │ 3 bins (FDR 2.3e-04)  │ A3SS (FDR 5.2e-11) │     YES     │
├─────────┼───────────────────────┼────────────────────┼─────────────┤
│ ABCC5   │  2 bins (FDR 0.042)   │ A3SS (FDR 2.0e-05) │     YES     │
├─────────┼───────────────────────┼────────────────────┼─────────────┤
│ GUSBP11 │ 2 bins (FDR 4.0e-11)  │   SE (FDR 0.049)   │     YES     │
├─────────┼───────────────────────┼────────────────────┼─────────────┤
│ ANKHD1  │   1 bin (FDR 0.006)   │   SE (FDR 0.006)   │     YES     │
├─────────┼───────────────────────┼────────────────────┼─────────────┤
│ ADAM12  │           —           │         —          │     no      │
└─────────┴───────────────────────┴────────────────────┴─────────────┘
5 of 6 target genes confirmed. The three strongest candidates (ABCC5, CRNDE, UQCC) are all
reproduced. ADAM12 was the weakest finding in the original paper.

## Directory Structure

Nextflow_Eye_Cancer_Workflow/
├── main.nf                 # Main workflow (Nextflow DSL2)
├── nextflow.config         # Parameters and process resources
├── samples.tsv             # Sample metadata
├── environment.yml         # Conda env: QC, trimming, alignment
├── splicing-env.yml        # Conda env: DEXSeq, rMATS (Python 3.12)
├── modules/
│   ├── fastqc.nf
│   ├── trimmomatic.nf
│   ├── multiqc.nf
│   ├── hisat2.nf           # HISAT2 align + samtools sort + Picard MarkDuplicates
│   ├── alignment_qc.nf
│   ├── dexseq.nf           # DEXSeq prepare annotation, count, analysis
│   └── rmats.nf            # rMATS alternative splicing detection
├── scripts/
│   └── run_dexseq.R        # DEXSeq differential exon usage analysis
├── data/                   # Symlink to FASTQ files
├── references/             # GRCh37 genome, GTF, HISAT2 index
├── results/                # All pipeline outputs
└── work/                   # Nextflow work directory

## Target Genes

Differential splicing expected in 6 genes (Furney et al. Table 2):
┌─────────┬─────────────────────────────────────────┬──────────────────┐
│  Gene   │             Splicing Event              │ Detection Method │
├─────────┼─────────────────────────────────────────┼──────────────────┤
│ ABCC5   │ Intron 5 retention                      │ DEXSeq + rMATS   │
├─────────┼─────────────────────────────────────────┼──────────────────┤
│ CRNDE   │ Alternative 3' splice site (exon 4)     │ DEXSeq + rMATS   │
├─────────┼─────────────────────────────────────────┼──────────────────┤
│ UQCC    │ Alternative terminal exons              │ DEXSeq + rMATS   │
├─────────┼─────────────────────────────────────────┼──────────────────┤
│ GUSBP11 │ Cassette exon 7                         │ rMATS            │
├─────────┼─────────────────────────────────────────┼──────────────────┤
│ ANKHD1  │ Alternative 3' splice site (exon 3)     │ rMATS            │
├─────────┼─────────────────────────────────────────┼──────────────────┤
│ ADAM12  │ Alternative terminal exons (exon 18/19) │ rMATS            │
└─────────┴─────────────────────────────────────────┴──────────────────┘
