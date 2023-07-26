# bti2023

## Pipeline Description and Function
This code provides a pipeline for the analysis of reads from the M82 and Pennellii tomatoes and identification of differentially expressed genes. It is intended to be run from a Linux command line. Python programs included for the analysis of specific data are also meant to be run from the command line, and are followed by code to prepare data for the programs and then run the programs from the command line. All python programs are included with multi-processing capability and as a single-thread version. 

This pipeline first takes in the raw reads from the M82 and Pennellii tomatoes, cleans them by removing adaptor, low quality, polyA, and rRNA sequences, then maps the clean reads to the reference genomes for both tomatoes. It then takes the BAM files produced by the mapping to the reference genome and, based on the average mismatch count for each read id for each parent, assigns the read ids to a corresponding parent. 

The pipeline provides two ways to extract the reads for each parent based on the read ids: using a provided python script or the gatk package. It then assigns the reads mapped to each parent to genes of that parent using the package FeatureCounts, and provides code to extract the mean fragment length for later normalization using VDM or FPKM.

## Note 
The final python script is an addendum used for the liftoff package to remove genes that have multiple pairings, and was not used in the read analysis pipeline. 
