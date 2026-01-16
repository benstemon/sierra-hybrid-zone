#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=excess_het_search
#SBATCH --array=0


# This is an array batch script
# Rather than the typical array submission, this has the array hard-coded
# Simply submit with sbatch scriptname.sh
# This script calculates heterozygosity with vcftools


##########################################################################################
# SETUP #
##########################################################################################

cd $SLURM_SUBMIT_DIR

# Source environments
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Load vcftools
module load vcftools

# Paths to input vcfs
# Imputed vcf has imputed genotypes and was built on SNPs with <5% missing data
# Standard vcf has no imputation and has max missingness of 50%
# get these from the online supplemental data repo
imputedvcf="gemmasnps-.05missing.imputed.vcf.gz"
standardvcf="finalfilter.vcf.gz"

# Paths to hybrid ID files (one line per sample ID)
hybrids="data/hybridsonly-25-75.txt"

# Path to SNPs to analyze
insnps="data/insnps_AFD1.txt"

##########################################################################################


# Use --hardy to calculate p-value for each site with HWE test
# All vcfs have MAF 0.05 already, so that is baseline
# Could make MAF more stringent if we think it is an issue
# Set up a little case in block to generate each of these simultaneously
case $SLURM_ARRAY_TASK_ID in
    0) vcftools --gzvcf $imputedvcf --hardy --keep $hybrids --positions $insnps --out hwe-imputed-hybrids-25-75 ;;
esac




