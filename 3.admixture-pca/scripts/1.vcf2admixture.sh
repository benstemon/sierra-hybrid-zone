#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=vcf2admixture


### This script converts a VCF input file to a plink bed file and filters SNPs based on LD thresholds. Then it runs ADMIXTURE.

cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Load the admixture conda environment
conda activate admixture


# Path to filtered vcf and other parameters for filtering LD, and threads
# Parameters for LD pruning pulled from Walsh et al. 2023 (https://github.com/joanam/scripts/blob/master/ldPruning.sh)
invcf="finalfilter.vcf.gz"
windowsize=50
slidesize=10
r2thresh=0.1
numthreads=8

processname="admixture_prep_LD-50-10-0.1"


##################################
# GENERATE ADMIXTURE INPUT FILES #
##################################

# Generate initial .bed file and produce list of SNPS to be filtered/retained based on LD thresholds
plink2 --vcf $invcf --make-bed --indep-pairwise $windowsize $slidesize $r2thresh --max-alleles 2 --out $processname --snps-only --double-id --allow-extra-chr --set-missing-var-ids @:#\$r:\$a --new-id-max-allele-len 63


# Generate the .bed file with proper SNPs extracted
plink2 --bed $processname.bed --bim $processname.bim --fam $processname.fam --extract $processname.prune.in --out extracted.$processname --make-bed --allow-extra-chr


# Reformat chromosome names to integers in bim file to make compatible with Admixture
mv extracted.$processname.bim extracted.$processname.bim.original
cat extracted.$processname.bim.original | sed 's/scaffold_//g' > extracted.$processname.bim



################################
######## RUN ADMIXTURE #########
################################

for K in 1 2 3 4 5;
do
    admixture --cv extracted.$processname.bed $K -j$numthreads | tee log${K}.out
done




