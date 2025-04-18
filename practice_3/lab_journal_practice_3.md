# Project 3. E.coli outbreak investigation 

Authors: Basova Victoria, Malysheva Polina 

## Table of content

[Part 1. Exploring the dataset](#part-1-exploring-the-dataset)

[Part 2. K-mer profile and genome size estimation (optional)](#part-2-k-mer-profile-and-genome-size-estimation-optional)

[Part 3. Assembling E. coli X genome from paired reads](#part-3-assembling-e-coli-x-genome-from-paired-reads)

[Part 4. Genome Annotation](#part-4-genome-annotation)

[Part 5. Finding the closest relative of E. coli X](#part-5-finding-the-closest-relative-of-e-coli-x)

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

To calculate the possible number of k-mers and predict the genome size after assembling we used jellyfish:

```
jellyfish count -m 31 -C -s 1G -o kmers_31.jf SRR292678_F.fastq SRR292678_R.fastq

jellyfish histo -o kmers_31.histo kmers_31.jf
```

N: Depth of coverage, M: Kmer peak, K: Kmer-size, L: avg read length T: Total bases

N = (M*L)/(L-K+1)

Genome_size = T/N

N = (66.9 * 90)/(90 – 31 + 1) = 6021 / 60 = 100,35 

Genome size is 4.93 Mbp (494.9 Mbp / 100)


## Part 3. Assembling E. coli X genome 

Use SPAdes for this task:

```
/bin/spades.py	-1	/data/SRR292862/SRR292862_S2_L001_R1_001.fastq.gz	-2	/data/SRR292862/SRR292862_S2_L001_R2_001.fastq.gz	-o	/spades_output_srr292862	
```
After running on single pair-end _SRR292678_, the output directory contained **contigs.fasta**, **scaffolds.fasta**, and other files. 

We need to check quality of assembled genomes via QUAST for _SRR292678_ (first run, single pair-end)

```
python quast.py -o /data/quast_out_srr292678 contigs.fasta
```
It produced some files. HTML reports are in this folder. 

After we ran SPAdes again by consolidating three libraries: SRR292678 as a paired ends,  SRR292862 and SRR292770 as a mate pairs.

QUAST 	report for single SRR292678:

- N50	105346

- N90	21026

- contigs	206

QUAST	report for all 3 libs:

- N50	335515

- N90	79998

- contigs	105

## Part 4. Genome Annotation

We used prokka to create a genomic annotation.

```
 prokka --outdir prokka_SRR292678/ --compliant spades_output_SRR292678/scaffolds.fasta
```

`--compliant` for force Genbank/ENA/DDJB compliance (default OFF)

## Part 5. Finding the closest relative of E. coli X

We need to locate 16S rRNA in the assembled E. coli X genome to find relatives of E. coli X. So we used barnap to predict 16S rRNA genes and save sequences of genes to .txt:  

```
barrnap prokka_SRR292678/PROKKA_12012024.fna --outseq filename.fasta
```

