# Project 5. H+, or how to build a perfect human 

Authors: Basova Victoria, Malysheva Polina 

## Table of content

[Part 1. Obtaining data](#part-1-obtaining-data)

[Part 2. File conversion](#part-2-file-conversion)

[Part 3. Origins, haplogroups](#part-3-origins-haplogroups)

[Part 4. Annotation - sex and eye colour](#part-4-annotation---sex-and-eye-colour)

[Part 5. Annotation of all SNPs, selection of clinically relevant ones using VEP (Variant Effect Predictor)](#part-5-annotation-of-all-snps-selection-of-clinically-relevant-ones-using-vep-variant-effect-predictor)

## Part 1. Obtaining data

We will work with raw 23andMe data from Mikhail Raiko (the proper human genome reference version is GRCh37) and [Manu Sporny](https://github.com/msporny/dna) (the proper human genome reference version is GRCh36)

## Part 2. File conversion 

For further analysis it is necessary to convert the raw 23andMe file to .vcf file using `PLINK v1.90b6.21 64-bit (19 Oct 2020)`

Installation

```bash
conda install bioconda::plink
```

Using 

```bash
plink --23file snp_teacher.txt --recode vcf --out snp_teacher_clean --output-chr MT --snps-only just-acgt
```

## Part 3. Origins, haplogroups

[mthap version 0.19c (2023-08-15)](https://dna.jameslick.com/mthap/)  online tool is used for haplogroup determination 

## Part 4. Annotation - sex and eye colour

A set of markers for eye and hair color are taken from the [article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3694299/): rs12896399, rs12913832, rs1545397, rs16891982, rs1426654, rs885479, rs6119471, rs12203592 

Command line to search for SNPs of interest:

```bash
grep -w -f snp_eye_color.txt snp_teacher.txt > eye_teacher_snp.txt
```

## Part 5. Annotation of all SNPs, selection of clinically relevant ones using VEP (Variant Effect Predictor)
