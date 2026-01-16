#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=vcf2bimbam


cd $SLURM_SUBMIT_DIR


##################
##### SETUP #####
##################

# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Path to qctool executable
qctool="qctool"

# Path to GEMMA executable
gemma="gemma-0.98.5"

# Path to a phenotype file -- can be anything
pheno1="GEMMA_INPUTFILE.modified_hue.txt"


###################
##### FILTER #####
###################

# Remove parent taxa and one sample with missing phenotypic data
conda activate mapping_etc

bcftools view gemmasnps-.05missing.imputed.vcf.gz -S ^excluded_samples_GWAS.txt -Oz -o gemmasnps-.05missing.imputed-NOPARENTS.vcf.gz

tabix gemmasnps-.05missing.imputed-NOPARENTS.vcf.gz


#####################
#### SNPs in LD ####
#####################

# Switch to conda env with plink2 installed
conda deactivate
conda activate admixture

# Generate initial list of SNPs to be filtered
# Windowsize 50, slidesize 10, r2 threshold 0.1
plink2 --vcf gemmasnps-.05missing.imputed-NOPARENTS.vcf.gz --make-bed --indep-pairwise 50 10 0.1 --max-alleles 2 --out LDpruned-NOPARENTS --snps-only --double-id --allow-extra-chr --set-missing-var-ids @:# --new-id-max-allele-len 63


###################
# CONVERT, FILTER #
###################

# Switch to conda env with proper libraries installed
conda deactivate
conda activate gwas

# Convert to bimbam format
$qctool -g gemmasnps-.05missing.imputed-NOPARENTS.vcf.gz \
 -excl-positions LDpruned-NOPARENTS.prune.out \
 -ofiletype bimbam_dosage \
 -og gemmasnps-.05missing-NOPARENTS-LDpruned.bimbam

# Generate some quick statistics for the input .vcf
wc -l "LDpruned-NOPARENTS.prune.in" > variant-counts-NOPARENT-LDSNPs-included.txt


# To use dosage from GP field
#$qctool -g gemmasnps-.05missing.vcf.gz -ofiletype bimbam_dosage \
# -vcf-genotype-field GP -permissive \
# -og gemmasnps-.05missing-dosage.bimbam


#####################
##### GT-Matrix #####
#####################

# Specify bimbam input file
bimbam="gemmasnps-.05missing-NOPARENTS-LDpruned.bimbam"

# Generate genomic relationship matrix
# This can be created with any phenotype matrix -- it does not change
$gemma -g $bimbam -p $pheno1 -gk 1 -outdir ./ -o GTMATRIX-NOPARENT-LDpruned


