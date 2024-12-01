# Project 3. E.coli outbreak investigation 

Authors: Basova Victoria, Malysheva Polina 

## Table of content

[Part 1. Exploring the dataset](#part-1-exploring-the-dataset)

[Part 2. K-mer profile and genome size estimation (optional)](#part-2-k-mer-profile-and-genome-size-estimation-optional)

[Part 3. Assembling E. coli X genome from paired reads](#part-3-assembling-e-coli-x-genome-from-paired-reads)

## Part 1. Exploring the dataset

We have three libraries from the TY2482 sample 

- SRR292678 - paired end, insert size 470 bp (forward reads, reverse reads, 400 Mb each)

- SRR292862 – mate pair, insert size 2 kb, (forward reads, reverse reads, 200 Mb each)

- SRR292770 – mate pair, insert size 6 kb, (forward reads, reverse reads, 200 Mb each)

**FastQC**

For all datasets we observe excellent quality of reads. Sequencing was successful. 

_SRR292678_

- Total Sequences	5499346

- Sequences flagged as poor quality	0

- Sequence length	90

- %GC	49

_SRR292770_

- Total Sequences	5102041

- Sequences flagged as poor quality	0

- Sequence length	49

- %GC	50

_SRR292862_

- Total Sequences	5102041

- Sequences flagged as poor quality	0

- Sequence length	49

- %GC	50

## Part 2. K-mer profile and genome size estimation (optional)

## Part 3. Assembling E. coli X genome from paired reads

Use SPAdes for this task:

```
/bin/spades.py	-1	/data/SRR292862/SRR292862_S2_L001_R1_001.fastq.gz	-2	/data/SRR292862/SRR292862_S2_L001_R2_001.fastq.gz	-o	/spades_output_srr292862	
```
After running on single pair-end _SRR292678_, the output directory contained **contigs.fasta**, **scaffolds.fasta**, and other files. 

We need to check quality of assembled genomes via QUAST for _SRR292678_ (first run, single pair-end)

```
python quast.py -o /data/quast_out_srr292678 contigs.fasta
```
It produced some files. HTML report is in this folder (quast_srr292678.html)

## Part 5. 

## Part 5. Genome Annotation

We used prokka to create a genomic annotation

Install

```
conda install -c conda-forge -c bioconda -c defaults prokka
```

Run

```
 prokka --outdir prokka_SRR292678/ --compliant spades_output_SRR292678/scaffolds.fasta
```

`--compliant` for force Genbank/ENA/DDJB compliance (default OFF)

## Part 6. Finding the closest relative of E. coli X

We need to locate 16S rRNA in the assembled E. coli X genome to find relatives of E. coli X. So we used barnap to predict 16S rRNA genes. 

