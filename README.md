# Uveal Melanoma SF3B1 Differential Splicing Pipeline
Nextflow DSL2 pipeline reproducing the RNA-seq analysis from
**Furney et al. (2013)** — *SF3B1 mutations are associated with alternative
splicing in uveal melanoma* (Cancer Discovery 3(10):1122-9).

Authors:

- Tanja Gesslbauer
- Jessica Kirchner
- Susanna Kummer
- Mathias Mayrgündter
- Konstantin Zilberburg

## Quick Start

```bash
# 1. Activate the base conda environment
conda activate uveal-pipeline

# 2. Run the pipeline (use screen for long runs)
screen -S nxf
nextflow run main.nf -resume

# when using screen
## Detach (leave it running) — press:
Ctrl + A, then D
## Close your terminal / disconnect SSH — job keeps running
## Later, reconnect and reattach
screen -r nxf
## List sessions if you forget the name
screen -ls

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
mv Homo_sapiens.GRCh37.75.gtf references/   # <-- add this

# Build HISAT2 index
mkdir -p references/grch37
hisat2-build Homo_sapiens.GRCh37.75.dna.primary_assembly.fa references/grch37/genome
```

### Raw Data

Download paired-end FASTQ files from SRA062359 into `data/`:

```bash
mkdir -p data

for SRR in SRR628582 SRR628583 SRR628584 SRR628585 SRR628586 SRR628587 SRR628588 SRR628589; do

    if [[ -f data/${SRR}_1.fastq.gz && -f data/${SRR}_2.fastq.gz ]]; then
        echo "[$(date +%T)] Skipping ${SRR} — already exists"
        continue
    fi

    echo "[$(date +%T)] Downloading ${SRR}..."
    fasterq-dump --split-files -O data/ ${SRR}

    echo "[$(date +%T)] Compressing ${SRR}..."
    gzip data/${SRR}_1.fastq data/${SRR}_2.fastq

    echo "[$(date +%T)] Done: ${SRR}"
    echo "---"
done

echo "[$(date +%T)] All samples complete."
```

## Configuration

### Configurable Parameters

Edit `nextflow.config` to match your setup:

| Parameter | Description | Default |
|---|---|---|
| `params.work_dir` | Nextflow work directory (use `/tmp` if disk-limited). Defaults to a user-specific path so multiple users can run the pipeline on the same server without permission conflicts. | `/tmp/${System.getProperty('user.name')}_nxf` |
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

### Process Resources

Edit the `process` block in `nextflow.config` to match your hardware:

| Process | CPUs | Memory | Notes |
|---|---|---|---|
| Default | 2 | 8 GB | Applied to all unlisted processes |
| `TRIMMOMATIC` | 4 | 16 GB | Threading ceiling is low; no benefit beyond 4 |
| `HISAT2_ALIGN` | 16 | 32 GB | Scales well with cores; `maxForks = 2` recommended on ≥32-core machines |
| `DEXSEQ_COUNT` | 1 | 32 GB | Single-threaded; memory may spike on large BAMs |
| `DEXSEQ_ANALYSIS` | 1 | 64 GB | R loads full count matrix; increase if OOM errors occur |
| `RMATS` | 8 | 32 GB | `maxForks = 2` safe on machines with ≥64 GB free RAM |

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


## 📜 License
 
This pipeline is released under the **MIT License**. See [`LICENSE`](LICENSE) for full details.
 
The bioinformatics tools invoked by this pipeline are subject to their own respective licenses.
 
---
 
## 📚 References
 
If you use this pipeline in your research, please cite the following foundational studies:
 
> **Furney SJ, et al.** (2013). SF3B1 mutations are associated with alternative splicing in uveal melanoma. *Cancer Discovery*, 3(10):1122–1129.  
> DOI: [10.1158/2159-8290.CD-13-0330](https://doi.org/10.1158/2159-8290.CD-13-0330) · PMID: [23861464](https://pubmed.ncbi.nlm.nih.gov/23861464/)
 
> **Harbour JW, et al.** (2013). Recurrent mutations at codon 625 of the splicing factor SF3B1 in uveal melanoma. *Nature Genetics*, 45(2):133–135.  
> DOI: [10.1038/ng.2523](https://doi.org/10.1038/ng.2523) · PMID: [23313955](https://pubmed.ncbi.nlm.nih.gov/23313955/)

---

Tools:

> **Anaconda Inc.** (2020). Anaconda Software Distribution. Vers. 2-2.4.0.
> URL: [https://docs.anaconda.com](https://docs.anaconda.com)

> **Anders S, Reyes A, Huber W.** (2012). Detecting differential usage of exons from RNA-seq data. *Genome Research*, 22(10):2008–2017.
> DOI: [10.1101/gr.133744.111](https://doi.org/10.1101/gr.133744.111) · PMID: [22722343](https://pubmed.ncbi.nlm.nih.gov/22722343/)

> **Bolger AM, Lohse M, Usadel B.** (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15):2114–2120.
> DOI: [10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170) · PMID: [24695404](https://pubmed.ncbi.nlm.nih.gov/24695404/)

> **Di Tommaso P, et al.** (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4):316–319.
> DOI: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820) · PMID: [28398311](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> **Ewels P, Magnusson M, Lundin S, Käller M.** (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19):3047–3048.
> DOI: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) · PMID: [27312411](https://pubmed.ncbi.nlm.nih.gov/27312411/)

> **Kim D, Paggi JM, Park C, Bennett C, Salzberg SL.** (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37(8):907–915.
> DOI: [10.1038/s41587-019-0201-4](https://doi.org/10.1038/s41587-019-0201-4) · PMID: [31375807](https://pubmed.ncbi.nlm.nih.gov/31375807/)

> **Li H, Handsaker B, Wysoker A, et al.** (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16):2078–2079.
> DOI: [10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352) · PMID: [19505943](https://pubmed.ncbi.nlm.nih.gov/19505943/)

> **Picard Toolkit.** (2019). Broad Institute, GitHub Repository.
> URL: [https://broadinstitute.github.io/picard](https://broadinstitute.github.io/picard)

> **Shen S, Park JW, Lu ZX, Lin L, Henry MD, Wu YN, Zhou Q, Xing Y.** (2014). rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *Proceedings of the National Academy of Sciences*, 111(51):E5593–E5601.
> DOI: [10.1073/pnas.1419161111](https://doi.org/10.1073/pnas.1419161111) · PMID: [25480548](https://pubmed.ncbi.nlm.nih.gov/25480548/)

---
 
## 🙏 Acknowledgements
 
This pipeline was developed with support from the [nf-core](https://nf-co.re/) community framework. We thank the authors of all tools used within the workflow.
