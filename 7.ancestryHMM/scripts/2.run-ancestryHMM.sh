#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=run_aHMM

# This script runs ancestryHMM

cd $SLURM_SUBMIT_DIR

##################
##### SETUP ######
##################


# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


# Activate conda environment with ancestryHMM installed
conda activate ancHMM


# Paths to the .ahmm file, the .ploidy file, and output directory
# Need to get these from the supplemental data repo
infile="allsamples-0.75-sierras-pruned.ahmm"
ploidyfile="ahmm-0.75-sierras.ploidy"
outdir="outfiles_0.75-sierras"


####################
##### ANALYSIS #####
####################

# Make the output directory if it doesn't already exist
mkdir -p $outdir
cd $outdir

# Run ancestryHMM
# a = # of ancestral pops and overall ancestry proportion of each in the sample
# I have ~61% newberryi, ~39% davidsonii. new = panel 0 and dav = panel 1
# p parameters are pulse for ancestry 0 and 1, respectively. Then admixture time and prop
# negative numbers allow for an estimate of admixture time and current proportion
ancestry_hmm -i $infile -s $ploidyfile -a 2 0.61 0.39 -p 0 -100 0.61 -p 1 -100 0.39 -e 1e-3

# Also run the Viterbi model to obtain tract lengths
ancestry_hmm -i $infile -s $ploidyfile -a 2 0.61 0.39 -p 0 -100 0.61 -p 1 -100 0.39 -e 1e-3 -v

