# Uveal Melanoma SF3B1 Differential Splicing Pipeline
Nextflow DSL2 pipeline reproducing the RNA-seq analysis from
**Furney et al. (2013)** — *SF3B1 mutations are associated with alternative
splicing in uveal melanoma* (Cancer Discovery 3(10):1122-9).

## Quick Start

```bash
# 1. Activate the base conda environment
conda activate uveal-pipeline

# 2. Run the pipeline (use screen for long runs)
screen -S nxf
nextflow run main.nf -resume

# 3. Results appear in results/
```

> **Important:** The `uveal-pipeline` conda environment must be active before
> running the pipeline. The Trimmomatic adapter path is resolved automatically
> from `$CONDA_PREFIX`.

## Requirements

### Conda Environments

Two conda environments are required:

| Environment | Python | Purpose |
|---|---|---|
| `uveal-pipeline` | 3.14 | FastQC, Trimmomatic, HISAT2, samtools, Picard, MultiQC |
| `splicing-env` | 3.12 | DEXSeq, HTSeq, rMATS, R packages (optparse, pheatmap, ggplot2) |

Create from exported YAML files:

```bash
conda env create -f environment.yml        # uveal-pipeline
conda env create -f splicing-env.yml       # splicing-env
```

> **Note:** rMATS and HTSeq require Python ≤3.12, which is why a separate
> environment is needed. The pipeline calls it via `conda run -n splicing-env`.

### Reference Data

Download GRCh37 (Ensembl release 75) into `references/`:

```bash
# Genome FASTA
wget https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

# GTF annotation
wget https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz

# Build HISAT2 index
mkdir -p references/grch37
hisat2-build Homo_sapiens.GRCh37.75.dna.primary_assembly.fa references/grch37/genome
```

### Raw Data

Download paired-end FASTQ files from SRA062359 into `data/`:

```bash
mkdir -p data
for SRR in SRR628582 SRR628583 SRR628584 SRR628585 SRR628586 SRR628587 SRR628588 SRR628589; do
    fasterq-dump --split-files -O data/ ${SRR}
    gzip data/${SRR}_1.fastq data/${SRR}_2.fastq
done
```

## Configuration

### Configurable Parameters

Edit `nextflow.config` to match your setup:

| Parameter | Description | Default |
|---|---|---|
| `params.work_dir` | Nextflow work directory (use `/tmp` if disk-limited) | `/tmp/kzilberburg_nxf` |
| `params.splicing_env` | Name of Python 3.12 conda env for DEXSeq/rMATS | `splicing-env` |
| `params.genome` | Path to GRCh37 reference FASTA | `references/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa` |
| `params.gtf` | Path to Ensembl 75 GTF | `references/Homo_sapiens.GRCh37.75.gtf` |
| `params.hisat2_index` | HISAT2 index base path | `references/grch37/genome` |
| `params.adapters` | Trimmomatic adapter FASTA | Auto-resolved from `$CONDA_PREFIX` |
| `params.fdr_threshold` | FDR cutoff for significance | `0.1` |
| `params.trimmed_read_length` | Read length after cropping | `99` |

### Key Configuration Notes

- **`work_dir`**: The pipeline generates large intermediate files (~60GB during alignment).
  If your home directory has limited quota, point this to a scratch partition like `/tmp`.
- **`adapters`**: Resolved via `${System.getenv('CONDA_PREFIX')}` — requires `uveal-pipeline`
  to be the active conda environment when launching the pipeline.
- **`splicing_env`**: Referenced by `conda run -n` inside processes. Change this if your
  Python 3.12 environment has a different name.

## Samples

8 uveal melanoma RNA-seq samples from Harbour et al. (SRA062359):

| Sample | Condition | SRA Accession | Alignment Rate |
|---|---|---|---|
| sample_1 | SF3B1 mutant | SRR628582 | 97.10% |
| sample_2 | SF3B1 mutant | SRR628583 | 97.14% |
| sample_3 | SF3B1 mutant | SRR628584 | 97.25% |
| sample_4 | SF3B1 mutant | SRR628589 | 96.53% |
| sample_5 | wildtype | SRR628585 | 97.10% |
| sample_6 | wildtype | SRR628586 | 97.05% |
| sample_7 | wildtype | SRR628587 | 96.80% |
| sample_8 | wildtype | SRR628588 | 97.23% |

## Pipeline Overview

```
raw FASTQ ──► FastQC (raw) ──────────────────────────────────┐
    │                                                        │
    ▼                                                        │
Trimmomatic (crop 99bp, adapter removal)                     │
    │                                                        │
    ├──► FastQC (trimmed) ───────────────────────────────────┤
    │                                                        │
    ▼                                                        │
HISAT2 align + samtools sort + Picard MarkDuplicates         │
    │                                                        │
    ├──► samtools flagstat + idxstats ───────────────────────┤
    │                                                        │
    ├────────────────┬──────────────┬────────────────┐       │
    ▼                ▼              ▼                ▼       │
DEXSeq count      rMATS      Coverage plots   PSI heatmap    │
    │                │              │                │       │
    ▼                │              │                │       │
DEXSeq analysis      │              │                │       │
    │                │              │                │       │
    └────────┬───────┘              │                │       │
             ▼                      │                │       │
    Splicing summary                │                │       │
             │                      │                │       │
             └──────────────────────┴────────────────┘       │
                                                             │
                      MultiQC (final) ◄──────────────────────┘
```

### Processes

| Process | Tool | Module |
|---|---|---|
| FASTQC_RAW / FASTQC_TRIMMED | FastQC | `modules/fastqc.nf` |
| TRIMMOMATIC | Trimmomatic 0.40 | `modules/trimmomatic.nf` |
| HISAT2_EXTRACT_SPLICESITES | HISAT2 | `modules/hisat2.nf` |
| HISAT2_ALIGN | HISAT2 + samtools + Picard | `modules/hisat2.nf` |
| ALIGNMENT_QC | samtools flagstat/idxstats | `modules/alignment_qc.nf` |
| DEXSEQ_PREPARE_ANNOTATION | DEXSeq (R/Bioconductor) | `modules/dexseq.nf` |
| DEXSEQ_COUNT | HTSeq/DEXSeq | `modules/dexseq.nf` |
| DEXSEQ_ANALYSIS | DEXSeq | `modules/dexseq.nf` |
| RMATS | rMATS turbo | `modules/rmats.nf` |
| SPLICING_SUMMARY | R (custom script) | `modules/visualization.nf` |
| COVERAGE_PLOTS | R/ggplot2 | `modules/visualization.nf` |
| SPLICING_HEATMAP | R/pheatmap | `modules/visualization.nf` |
| MULTIQC_RAW / MULTIQC_TRIMMED / MULTIQC_FINAL | MultiQC | `modules/multiqc.nf` |

## Results

### Quality Control (Phase 2)

- **Read survival:** 96.7% average after Trimmomatic
- **Read length:** Cropped from 101bp to 99bp (rMATS compatibility)
- Adapter content removed (TruSeq3-PE-2)

### Alignment (Phase 3)

- **Aligner:** HISAT2 with splice-site awareness (GRCh37/Ensembl 75)
- **Alignment rate:** All 8 samples >96.5%
- Duplicate marking with Picard MarkDuplicates (flagged, not removed)

### Differential Splicing (Phase 4)

**DEXSeq:** 222 significant exonic bins (FDR < 0.1)
**rMATS:** 1,322 significant alternative splicing events (FDR < 0.1)

#### Target Gene Reproduction (Furney et al. Table 1)

| Gene | DEXSeq | rMATS | Reproduced? |
|---|---|---|---|
| ABCC5 | 20 bins (FDR 1.12e-14) | SE (FDR 0.0) | Yes |
| CRNDE | 3 bins (FDR 7.88e-06) | SE (FDR 0.0) | Yes |
| UQCC1 | 31 bins (FDR 5.66e-27) | SE, A3SS (FDR 0.0) | Yes |
| GUSBP11 | 2 bins (FDR 4.05e-11) | SE (FDR 7.5e-05) | Yes |
| ANKHD1 | 0 bins | A5SS (FDR 0.015) | Partial (rMATS only) |
| ADAM12 | Not testable | Not expressed | No |

**5 of 6 target genes reproduced.** 4 confirmed by both methods.

- **ADAM12** was not expressed in these samples (absent from rMATS output entirely),
  consistent with the original study noting variable expression across cohorts.
- **GUSBP11** has two Ensembl IDs in GRCh37 (ENSG00000228315 and ENSG00000214265);
  the pipeline searches for both.

### Visualization (Phase 5)

Generated figures in `results/figures/`:

| File | Description |
|---|---|
| `coverage_ABCC5.pdf` | Read coverage across ABCC5 exons (mutant vs wildtype) |
| `coverage_CRNDE.pdf` | Read coverage across CRNDE exons |
| `coverage_UQCC1.pdf` | Read coverage across UQCC1 exons |
| `splicing_heatmap.pdf` | PSI heatmap of top differentially spliced events |
| `splicing_summary.tsv` | Combined DEXSeq + rMATS target gene summary table |

## Directory Structure

```
Nextflow_Eye_Cancer_Workflow/
├── main.nf                      # Main workflow
├── nextflow.config              # Pipeline configuration
├── samples.tsv                  # Sample metadata (sample_id, condition, sra_accession)
├── environment.yml              # uveal-pipeline conda env export
├── splicing-env.yml             # splicing-env conda env export
├── modules/
│   ├── fastqc.nf                # FastQC (reused for raw + trimmed)
│   ├── trimmomatic.nf           # Read trimming and adapter removal
│   ├── multiqc.nf               # MultiQC report aggregation
│   ├── hisat2.nf                # Splice-aware alignment + deduplication
│   ├── alignment_qc.nf          # samtools flagstat + idxstats
│   ├── dexseq.nf                # Differential exon usage analysis
│   ├── rmats.nf                 # Alternative splicing event detection
│   └── visualization.nf         # Summary tables, coverage plots, heatmaps
├── scripts/
│   ├── run_dexseq.R             # DEXSeq analysis script
│   ├── summarize_splicing.R     # Combine DEXSeq + rMATS results
│   ├── plot_coverage.R          # Per-gene coverage plots
│   └── plot_heatmap.R           # PSI heatmap across samples
├── data/                        # Raw FASTQ files (not tracked in git)
├── references/                  # Genome, GTF, HISAT2 index (not tracked)
└── results/                     # All pipeline outputs
    ├── fastqc_raw/
    ├── fastqc_trimmed/
    ├── multiqc_raw/
    ├── multiqc_trimmed/
    ├── multiqc_final/
    ├── trimmed/
    ├── aligned/
    ├── dedup/
    ├── alignment_qc/
    ├── dexseq/
    ├── rmats/
    └── figures/
```

## Troubleshooting

| Problem | Solution |
|---|---|
| Disk quota exceeded during alignment | Set `params.work_dir` to `/tmp/` or another scratch partition |
| SSH disconnect kills pipeline | Run inside `screen -S nxf` |
| R script changes not picked up by `-resume` | Add a comment change (e.g., `# v2`) to the process script block |
| rMATS fails with Python error | Ensure `splicing-env` uses Python ≤3.12 |
| DEXSeq count files have wrong format | Preprocessing step removes quotes and summary lines automatically |
| `$CONDA_PREFIX` not set | Activate `uveal-pipeline` before running: `conda activate uveal-pipeline` |

## Reference

Furney SJ, Pedersen M, Gentien D, et al. SF3B1 mutations are associated
with alternative splicing in uveal melanoma. *Cancer Discovery*. 2013;
3(10):1122-1129. doi:10.1158/2159-8290.CD-13-0330
