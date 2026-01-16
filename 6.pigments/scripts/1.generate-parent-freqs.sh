#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=AFD_parents_sierras

cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with pixy installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


module load vcftools


invcf="gemmasnps-.05missing.imputed.vcf.gz"

# davidsonii
vcftools --gzvcf $invcf \
  --keep davidsonii-fst-pops.txt \
  --freq \
  --out dav_allele_freq

# newberryi
vcftools --gzvcf $invcf \
  --keep newberryi-fst-pops.txt \
  --freq \
  --out new_allele_freq


