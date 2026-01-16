#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=prep_bgchm

# This script prepares the input files for bgchm

cd $SLURM_SUBMIT_DIR

##################
##### SETUP #####
##################


# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Activate conda environment with bcftools and pandas installed
conda activate mapping_etc

# Path to the input .vcf and the reference genome index file
invcf="gemmasnps-.05missing.imputed.vcf.gz"

# Path to the file containing chr and ps information for SNPs to be included
insnps="data/subsampled-maf0.1-bgcinput.txt"

# Path to the parentfile and samplelist for the vcf-to-bgchm python script
# One line per sample name
# Sample list is ALL samples, including parents
# parent list is just the parents (both p0 and p1)
samplelist="data/samplelist-bgchm-sierras.txt"
parentlist="data/parentlist-bgchm-sierras.txt"

# OPTIONAL: if using --p0file and/or p1file option in vcf2bgchm script, paths to
# P0 and P1 files. One line per sample name
p0file="data/p0file_newberryi-sierras.txt"
p1file="data/p1file_davidsonii-sierras.txt"


###########################
##### FILE GENERATION #####
###########################

# Use SNP information to filter the imputed vcf
bcftools view -T $insnps $invcf -Oz -o bgchmsnps-subsampled-maf0.1.vcf.gz

# Use custom python script to convert the .vcf to bgchm format (known genotypes)
# Columns are SNPs, rows are samples, and genotypes are 0, 1, and 2
# Also need separate files for hybrids, parent 1, and parent 2
python vcf2bgchm-gzipupdate_v2.py --vcf bgchmsnps-subsampled-maf0.1.vcf.gz --samples $samplelist --pure_parents $parentlist --p0file $p0file --p1file $p1file


