# Project 7. Dead Man’s Teeth. Introduction to metagenomics analysis

Authors: Basova Victoria, Malysheva Polina 

## Table of content

[Part 1. Obtaining data](#part-1-obtaining-data)

[Part 2. 16S rRNA analysis using dada2 script](#part-2-16s-rrna-analysis-using-dada2-script)

## Part 1. Obtaining data

The data was downloaded from the SRA archive via ftp server. SRA IDs: SRR957750, SRR957753, SRR957756, SRR957760, SRR986773, SRR986774, SRR986778, SRR986779, SRR986782

## Part 2. 16S rRNA analysis using dada2 script

Link to script [here](https://benjjneb.github.io/dada2/tutorial.html)

1. QC of raw data: `plotQualityProfile`

2. Trim low-quality reads: `filterAndTrim`, options trimLeft = 32, truncLen=140,  truncQ = 5. Some reads are low-quality, amplicon sequencing data often keeps the artificial sequences (barcodes). We trimmed it using options (primer + adapter = 35bp, amplicon = 145 bp, so m = 35 and n = 140)

3. Learn the Error Rates: `learnErrors`. Train the error model on our data. The DADA2 algorithm uses a parametric error model, and each amplicon data set has its own set of error rates. 

4. Dereplication: `derepFastq` and `makeSequenceTable`. Combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. Output: abundance table

5. Bimera removal: `removeBimeraDenovo`. Removal of chimeric sequences, a sequence fragment that is erroneously created during the PCR amplification process by the artificial joining of two or more unrelated DNA fragments.

6. Taxonomy assignment: `assignTaxonomy` and `addSpecies`. Taxonomic classification of our sequences. First, we do this on the higher level through naive Bayesian classifier in dada2, then we go for assignment of species-level annotations to reads in pairs. SILVA Database downloaded from [here](https://zenodo.org/records/1172783https://zenodo.org/records/1172783)

7. Extract data to MicrobiomeAnalyst: made ASV table and metadata for online server

## Part 3. Whole-metagenomic analysis

For this part we took data from G12 individual. 

1. Shotgun sequence data profiling using Kraken2

```bash
conda install kraken2
```
```bash
kraken2-build --standard --db $g12_ind

kraken2 --db /g12_ind --output out_file --report rep_file /data/SRR957742.fastq
```

2. Visualization of the Kraken results as a Sankey diagram

Use online-server [Pavian metagenomics data explorer](https://fbreitwieser.shinyapps.io/pavian/#)
