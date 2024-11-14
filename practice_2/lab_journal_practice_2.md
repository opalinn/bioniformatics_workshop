# Project 2. Why did I get the flu? Deep sequencing, error control, p-value, viral evolution.

Authors: Basova Victoria, Malysheva Polina 

## Table of content

[Part 1. Get and check data](#part-1-get-and-check-data)

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

Get reference sequence

Index the reference file, align your roommate’s viral data to the reference sequence and make an mpileup. 

 samtools flagstat alignment.bam
 
359374 + 0 in total (QC-passed reads + QC-failed reads)

356345 + 0 primary

0 + 0 secondary

3029 + 0 supplementary

0 + 0 duplicates

0 + 0 primary duplicates

358118 + 0 mapped (99.65% : N/A)

355089 + 0 primary mapped (99.65% : N/A)

0 + 0 paired in sequencing

0 + 0 read1

0 + 0 read2

0 + 0 properly paired (N/A : N/A)

0 + 0 with itself and mate mapped

0 + 0 singletons (N/A : N/A)

0 + 0 with mate mapped to a different chr

0 + 0 with mate mapped to a different chr (mapQ>=5)

