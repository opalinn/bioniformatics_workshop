# Project 6. Differential RNA expression analysis

Authors: Basova Victoria, Malysheva Polina 

## Table of content

[Part 1. Obtaining data](#part-1-obtaining-data)

[Part 2. Data preparation for analysis](#part-2-data-preparation-for-analysis)

[Part 3. Data analysis using DESeq2](#part-3-data-analysis-using-deseq2)

## Part 1. Obtaining data

The data was downloaded from the SRA archive via ftp server. SRA IDs: SRR941816, SRR941817, SRR941818, SRR941819.

Reference genome and annotation file [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/)

## Part 2. Data preparation for analysis

### Data filtering

FastQC was used to test the quality of the reads:

```bash
fastqc -o /path/to/output /practice_6/filtered_raw_data/SRR94181*
```

Fastp was used to filter low quality reads and trim adapter sequences:

```bash
fastp -i SRR941819.fastq -o /home/polinam/BI_practice/practice_6/filtered_raw_data/SRR941819_f.fastq
```

### Aligning with HISAT2

Build index:

```bash
hisat2_extract_splice_sites.py GCF_000146045.2_R64_genomic.gtf > genome.ss
hisat2_extract_exons.py GCF_000146045.2_R64_genomic.gtf > genome.exon
hisat2-build -p 16 GCF_000146045.2_R64_genomic.fna genome
```

Align reads and get .bam file for each data(in single-end mode):

```bash
hisat2 -p 16 -x genome -U filtered_raw_data/SRR941816_f.fastq | samtools sort > SRR941816_out.bam 
```

### Quantifying with featureCounts

Use featureCounts to obtain counts for each gene:

```bash
featureCounts -g gene_id -a GCF_000146045.2_R64_genomic.gtf -o srr18_feature_counts.txt SRR941818_out.bam
```

## Part 3. Data analysis using DESeq2 

R script is [here](https://figshare.com/articles/software/Scripts_for_RNA-seq_project/14239304)
