# olea-europaea-shade-uvb-rnaseq

Repository with the scripts necessary to analyze olive tree RNA-seq data from a Rubio de Casas Lab experiment.

Parting from the Illumina PE sequencing raw reads (not publicly available yet), the scripts can be executed line by line in this order:

## 1. RNAseq_Olivo_salmon.sh

This bash script sets a conda environment to perform QC, trimming, and pseudoalignment to genome of the RNA-seq libraries.

## 2. RNAseq_olivo_DGEA.R

This R script performs differential gene expression analysis, and generates plots and DEG lists for downstream analyses.

## 3. RNAseq_olivo_WGCNA.R

This R script performs co-expression network analysis, extracting co-expressed gene modules.

## 4. RNAseq_olivo_GOTEA.R

This R script performs GO term enrichment analysis from the DEG and gene module lists, producing plots and tables to summarize the enriched GO term information.

## RNAseq_olivo_functions.R

An auxiliary R script with functions used by other R scripts in the repository.

