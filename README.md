# Variant-Calling-Amplicon
Variant Calling Pipeline (Amplicon/Targeted Sequencing)

This repository contains a Bash pipeline that automates variant calling from paired-end amplicon or targeted sequencing data, using industry-standard tools like Fastp, BWA, GATK, and SAMtools. It processes both control and case samples, performs joint genotyping, and outputs high-quality VCF files for downstream analysis.

**What the Pipeline Does**

Trims reads using fastp
Aligns reads to the reference genome using BWA-MEM
Sorts and indexes BAM files with samtools
Adds read group info using GATK
Calculates coverage using bedtools
Runs GATK HaplotypeCaller per sample
Performs joint genotyping using GATK GenotypeGVCFs
Outputs final VCF file with variants across all samples

