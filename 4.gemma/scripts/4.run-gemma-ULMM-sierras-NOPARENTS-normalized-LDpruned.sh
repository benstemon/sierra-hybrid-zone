#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=ULMM
#SBATCH --array=0-25

# NOTE: This is an array script. The array is the number of files in the input_phenotype_matrices file -1. (25)


cd $SLURM_SUBMIT_DIR


##################
##### SETUP #####
##################

# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Load gwas conda environment
conda activate gwas

# Path to GEMMA executable
gemma="/work/bs66/software/GEMMA/gemma-0.98.5"

# Specify bimbam input file
# Remember -- the LD pruning file is ONLY for the relatedness matrix
# NOT for the ULMM
bimbam="gemmasnps-.05missing-NOPARENTS.bimbam"

# Specify the relatedness matrix
GTmatrix="GTMATRIX-NOPARENT-LDpruned.cXX.txt"

# Specify the covariates file
# 1st columns all 1s, then the second PC axis
#covariates="GEMMA-sierras-PC1only-covariate.txt"

#specify the phenotype file -- in same order as .vcf headers
phenofiles=(phenotypes_noparents_normalized/GEMMA_INPUTFILE*)
phenotype="${phenofiles[$SLURM_ARRAY_TASK_ID]}"


################
##### ULMM #####
################

# Get the name of the phenotype from the file
phenoname="${phenotype/phenotypes_noparents_normalized\/GEMMA_INPUTFILE./}"
phenoname="${phenoname/.txt/}"

# Run Gemma ULMM model -- Wald, LRT, and Score test, with covariates
$gemma -g $bimbam \
 -p $phenotype -n 2 \
 -k $GTmatrix \
 -lmm 4 \
 -outdir ./output_NOPARENTS_normalized_LDpruned \
 -o ulmm.NOPARENTS_normalized_LDpruned.$phenoname


