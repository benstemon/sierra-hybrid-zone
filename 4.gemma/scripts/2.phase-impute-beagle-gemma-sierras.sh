#!/bin/sh

#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p wessinger-48core
#SBATCH --job-name=impute


cd $SLURM_SUBMIT_DIR

# This script runs beagle on the VCF file to impute missing genotypes for GEMMA association mapping


##################
##### SETUP #####
##################


# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Activate conda environment with bcftools installed and load java
conda activate mapping_etc
module load java 

# Path to .vcf file to prune
invcf="/work/bs66/sierra_hybrids/VCFs-sierras/finalfilter.vcf.gz"

# Prefix for output (filters and imputation)
myprefix="gemmasnps-.05missing"

# Path to the compiled local beagle software
beagle="/work/bs66/software/beagle/beagle.22Jul22.46e.jar"


##################
##### FILTER #####
##################


#apply missingness filter (no more than 5% missing data)
bcftools filter $invcf --threads 12 -e 'F_MISSING>0.05' -Oz -o $myprefix.vcf.gz


##################
##### IMPUTE #####
##################

#run beagle to phase and impute genotypes with default settings
java -jar $beagle nthreads=12 gt=$myprefix.vcf.gz out=$myprefix.imputed


