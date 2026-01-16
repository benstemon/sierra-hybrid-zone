#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=bgchm-parallel-sierras
#SBATCH --array=0-9%10
#SBATCH --mail-type=END
#SBATCH --output=/dev/null

# NOTE: This is an array script. The array is the number of jobs you need to submit, contingent on the batch size in SNPs. As an example, if you have 25,000 SNPs and have set batch_size to 10,000, you would need to submit three jobs. I have 49019 SNPs, and am doing batches of 5000. THE ARRAY MUST START WITH 0.


cd $SLURM_SUBMIT_DIR

# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Activate conda environment with R and bgchm installed
conda activate r-env


# Specify a file that has all of the SNPs for the analysis
snpsfile="genotypes-uncoded-allsamples.tsv"


# Specify the batch size
batch_size=5000


# Find total SNPs
totalsnps=$(( $(awk -F'\t' '{print NF; exit}' "$snpsfile") - 1 ))


# Find the starting and ending SNP, based on the array submission
snpstart=$(( SLURM_ARRAY_TASK_ID * batch_size + 1 ))
snpend=$(( (SLURM_ARRAY_TASK_ID + 1) * batch_size ))
batchnumber=$((SLURM_ARRAY_TASK_ID + 1))

# Adjust snpend if it exceeds totalsnps
if [ $snpend -gt $totalsnps ]; then
  snpend=$totalsnps
fi


# Print these to the slurmfile for record
echo $snpstart
echo $snpend
echo $batchnumber


# Submit the Rscript with parameters for snpstart, snpend, and batch number
Rscript --vanilla 4b.run-bgchm-parallel-sierras.R "$snpstart" "$snpend" "$batchnumber"

