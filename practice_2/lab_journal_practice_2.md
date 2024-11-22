# Project 2. Why did I get the flu? Deep sequencing, error control, p-value, viral evolution.

Authors: Basova Victoria, Malysheva Polina 

## Table of content

[Part 0. Snakemake](#part-0-snakemake)

[Part 1. Get and inspect data](#part-1-get-and-inspect-data)

[Part 2. Alignment](#part-2-alignment)

[Part 3. Look for common and rare variants with VarScan](#part-3-look-for-common-and-rare-variants-with-varscan)

[Part 4. Work with the control sample sequencing data](#part-4-work-with-the-control-sample-sequencing-data)

[Part 5. Compare the control results to patient's results](#part-5-compare-the-control-results-to-patients-results)

## Part 0. Snakemake
All steps executed during this study are available in `Snakefile`.

### How to use

- create a directory containing reference `.fasta` file, samples `.fastq.gz` files and `Snakefile`
- create a virtual environment and install all requirements
- create a `config.yaml` file with actual names of reference and samples
Example:
```yaml
reference: sequence.fasta
samples:
  - sample_0
  - sample_1
  - sample_2
  - sample_3
```
- run snakemake 
```bash
snakemake -p --cores <preferable num of cores>
```

### Results
- The program creates a `.vcf` file for each sample and a `.txt` file with a table, which consists position, nucleotide substitution and frequence from each `.vcf` file
- The program also creates all intermediate files in a working directory

### Requirements
- Snakemake
- bwa
- samtools
- seqkit
- bcftools
- varscan

## Part 1. Get and inspect data

Download files:

-  SRR1705851 archive 

Inspecting data:

For SRR:

- 358265 reads

Seqkit Stats:

- type: DNA

- num_seqs: 358265

- sum_len: 52717864

- min_len: 35

- avg_len: 147.1

- max_len: 151

FastQC Report

- Per base sequence quality. We can see that the quality of the reads is high, but

![image](https://github.com/user-attachments/assets/0a686443-04bc-4436-bbfc-3ed0bdb0f92f)

- Per base sequence content. There are some problems here because of a large number of overrepresented sequences

![image](https://github.com/user-attachments/assets/8e5f52cf-ed66-4d75-ac97-8e1ed753dea5)

- Overrepresented sequences. Here we see a long list of such sequences. It is variants of hemagglutinin (HA) gene

![image](https://github.com/user-attachments/assets/bab31fcc-1ae4-4d0c-8e9c-ce62925e749d)

- Adapter Content. No adapters found

## Part 2. Alignment

All steps to create a .mpileup file are described in details in [lab journal #1](https://github.com/opalinn/bioniformatics_workshop/blob/main/practice_1/lab_journal_practice_1.md)

## Part 3. Look for common and rare variants with VarScan

Since our variants can be quite rare, we set the depth limit with the -d flag equal to 32000. That's more than the average coverage (30206.4)

```
samtools mpileup -d 30200 -f reference_ha_gene.fasta alignment_sorted.bam > p2.mpileup
```

After looking for common variants with --min-var-freq 0.95 we found 5 mutations. With a parameter --min-var-freq 0.001 we found 18 rare mutations. 

To get information from the .vcf file about the place where the mutation occurred, we use one of the following two methods: 

```
bcftools query -f '%POS\t%REF\t%ALT\t[%FREQ]\n' ali_1_001.vcf > ali_1_freq.txt
```
```
cat p2.vcf | awk 'NR>24 {print $1, $2, $4, $5}'
```

## Part 4. Work with the control sample sequencing data

Next, we needed to analyze three control samples. The work with them is the same as with the previous files: from quality check to aligning the reads to a reference and extract some data from .vcf files.

## Part 5. Compare the control results to patient's results

This step was done using excel spreadsheets. The file with the table presents mutation frequencies for 4 samples (patient and 3 controls), the mean of them, the standard deviation and the value equal to the sum of the mean and three standard deviations were calculated. 







