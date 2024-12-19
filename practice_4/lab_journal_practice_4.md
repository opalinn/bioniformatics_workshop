# Project 4. Tardigrades: from genestealers to space marines 

Authors: Basova Victoria, Malysheva Polina 

## Table of content

[Part 1. Obtaining data](#part-1-obtaining-data)

[Part 2. Structural annotation](#part-2-structural-annotation)

[Part 3. Physical localization](#part-3-physical-localization)

[Part 4. Localization prediction](#part-4-localization-prediction)

[Part 5. BLAST search and Pfam prediction](#part-5-blast-search-and-pfam-prediction)

## Part 1. Obtaining data

We will be using an assembled genome of the Ramazzottius varieornatus (YOKOZUNA-1 strain). Accession number GCA_001949185.1

## Part 2. Structural annotation

Files with al proteins and .gff file were made with AUGUSTUS. 

## Part 3. Physical localization 

Our goal is to find out which proteins from the mass spectrometry results correspond to the proteins of the R. varieornatus. 
To do this, we created a local database using protein data from the genome using blast:

```
makeblastdb -in data/augustus_aa.fasta -dbtype prot -out tardigrad_prot_db
```

After we perfomed a search: 

```
blastp -db tardigrad_prot_db -query mass_spec_proteins.fasta -outfmt 6  -out blast_res_on_massspec.txt
```

Next, we need to extract the complete amino acid sequences from the original protein file. To do this, we use seqtk

```
seqtk subseq augustus_proteins.fasta blast_res_on_massspec.txt > output_file.fasta
```

## Part 4. Localization prediction

We predicted where these proteins are located in the cell based on their sequences using two web servers: WoLF PSORT and TargetP. 

WoLF PSORT predicts cellular localization of proteins based on the presence of a signal peptide on their N-terminus. 

TargetP also predicts the subcellular localization of eukaryotic proteins. The location assignment is based on the predicted presence of any of the N-terminal presequences: chloroplast transit peptide (cTP), mitochondrial targeting peptide (mTP) or secretory pathway signal peptide (SP). 

## Part 5. BLAST search and Pfam prediction



















