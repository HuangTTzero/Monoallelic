# Random Monoallelic Expression Identification

This repository contains scripts and tools for identifying random monoallelic expression (RME) in sequencing data, focusing on allele-specific analysis using SNP-masked references and sequence assignment to two genomes. The project involves steps from SNP filtering to sequence processing, culminating in the generation of FASTQ files for downstream analysis with CellRanger.

## 1. SNP-Masked Reference Creation

### 1.1 Filter SNPs for High Confidence Variants

This script filters the latest VCF file to retain high-confidence SNPs from various strains (e.g., 129P2 and PWK).

```bash
perl filter_SNPs_129P2.pl
perl filter_SNPs_PWK.pl
```
### 1.2 Obtain PWK-Specific SNPs for Allele-Specific Read Assignment

This step extracts PWK-specific SNPs to be used for assigning reads to specific alleles.

```bash
Rscript get_pwk_specific.R
python snp2align.py
```
### 1.3 Merge SNPs from 129P2 and PWK to Create a Combined VCF File

This script merges the SNPs from the 129P2 and PWK strains into a single VCF file, which will be used for masking the reference genome.

```bash
bash get_merge_snp.sh
```
### 1.4 Create the SNP-Masked Reference Genome Using SNPs from Both Strains

This step generates a SNP-masked reference genome by applying the combined SNP data from both strains.

```bash
bash get_masked_ref.sh
```
## 2. Sequence Assignment to Two Genomes

Once the SNP-masked reference is generated, you will assign sequencing reads to two different genomes (e.g., PWK and 129P2) based on the SNP data.
```bash
bash split_ATAC.sh
bash split_RNA.sh
```
## 3. Generate FASTQ Files from BAM Files for Allelic Specific CellRanger Output

This step converts BAM files into FASTQ files, with sequencing reads assigned to specific alleles. These FASTQ files are required for running CellRanger for downstream analysis of allelic expression.
```bash
bash bam2fastq_ATAC.sh
bash bam2fastq_RNA.sh
```
## 4. Data Analysis Code
The data analysis scripts for generating figures and processing data can be found in the Figure folder of this repository. These scripts include visualizations and analyses that support the results of this project.

## 5. Software for BAM File Assignment
The scripts used for assigning reads in BAM files to specific alleles come from the following repositories:

https://github.com/sandberg-lab/Smart-seq3/tree/master/allele_level_expression
https://github.com/jinxu9/AlleleSpecificATACseq

We have made simple modifications to these original scripts to adapt them for use with our hybrid progenitor genomes (129P2 and PWK). The modified versions are included in the software folder.
