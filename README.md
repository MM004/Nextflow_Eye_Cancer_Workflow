# Nextflow Eye-Cancer Workflow

## Papers
### 1. SF3B1 mutations are associated with alternative splicing in uveal melanoma
- https://pubmed.ncbi.nlm.nih.gov/23861464/
### 2. Recurrent mutations at codon 625 of the splicing factor SF3B1 in uveal melanoma 
- https://pubmed.ncbi.nlm.nih.gov/23313955/

## Create Download Script for Data in 2.Paper

## Recrate Workflow of 1.Paper
 # 1. Clone the repository
  git clone <repo-url>
  cd Nextflow_Eye_Cancer_Workflow

  # 2. Create conda environment
  conda env create -f environment.yml
  conda activate uveal-pipeline

  # 3. Download references
  mkdir -p references && cd references
  wget https://genome-idx.s3.amazonaws.com/hisat/grch37_genome.tar.gz
  wget https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
  wget https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
  tar -xzf grch37_genome.tar.gz
  gunzip *.gz
  cd ..

  # 4. Symlink or download FASTQ data into data/

  # 5. Run the pipeline
  nextflow run main.nf -resume

  Pipeline Overview

  1. QC: FastQC on raw reads
  2. Trimming: Trimmomatic (crop to 99bp, adapter removal, quality trim)
  3. QC: FastQC on trimmed reads + MultiQC reports
  4. Alignment: HISAT2 to GRCh37
  5. Deduplication: Picard MarkDuplicates
  6. Splicing analysis: DEXSeq + rMATS
  7. Visualization: Coverage plots, splicing heatmap
