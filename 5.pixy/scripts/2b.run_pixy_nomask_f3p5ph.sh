#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=pixy-parents-nomask

cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with pixy installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc


#paths to (1) the input vcf, (2) the desired out-directory, and (3) the populations file
invcf="filtered_allsites_PARENTS.vcf.gz"
popfile="pixy_parentpops_sierras.txt"


#Run pixy f3'5'h zoom 1
pixy --stats pi fst dxy \
 --vcf $invcf \
 --populations $popfile \
 --window_size 100 \
 --n_cores 8 \
 --output_folder "pixyout-nomask_f3p5ph_zoom1" \
 --output_prefix f3p5ph_zoom1 \
 --chromosomes "scaffold_1086" \
 --interval_start 36500000 \
 --interval_end 44000000
