#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=PCA


# This script uses a previously-generated LD-filtered .bed file to conduct a genomic PCA in plink2

cd $SLURM_SUBMIT_DIR


#########
# SETUP #
#########


# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Load the admixture conda environment -- has plink2 installed
conda activate admixture

# Path to the LD-pruned file prefix
# This file removed SNPs with R2 > 0.1 in 50 SNP windows, sliding every 10 SNPs
infile="/work/bs66/sierra_hybrids/admixture/extracted.admixture_prep_LD-50-10-0.1"


################
# RUN ANALYSIS #
################

# PCA, with allele loadings written
plink2 --bfile $infile --pca biallelic-var-wts --out PCA --allow-extra-chr
