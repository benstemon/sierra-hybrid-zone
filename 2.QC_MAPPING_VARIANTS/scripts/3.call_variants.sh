#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=variants_bcftools


cd $SLURM_SUBMIT_DIR


##################
##### SETUP #####
##################

########################
# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Activate conda environment with packages installed
# Needs bcftools (v1.15.1 works)
conda activate mapping_etc

# Load other modules
module load vcftools



# Hard path to the base directory, reference genome, and mapped, filtered reads
# Ensure basedir has "/" at end
basedir="sierra_hybrids/"
refgenome="annot_Pdavidsonii_genome.1mb.fasta"
mapped_filtered_reads="mapped_filtered_bams_sierras"
newsamplenames="samplename_changer.txt"


# Change this to match number of cores allocated
numthreads=16

########################


##################
# CALL GENOTYPES #
##################

# Make list of mapped, filtered bams and out directory
bamlist=($mapped_filtered_reads/*.bam)
mkdir -p ${basedir}VCFs-sierras
cd ${basedir}VCFs-sierras


# Call variants (GT likelihoods), including indels, and normalize
bcftools mpileup --threads $numthreads -Ou \
 -a FORMAT/AD,FORMAT/DP \
 --min-BQ 20 --min-MQ 30 \
 -f $refgenome ${bamlist[@]} | 
 bcftools call --threads $numthreads -mv -f GQ,GP -Ou |
 bcftools norm --threads $numthreads -f $refgenome -Oz -o unfiltered_vcf.gz


# Change sample names in the header, to truncated version
bcftools reheader -s $newsamplenames --threads $numthreads unfiltered_vcf.gz -o tmp.vcf.gz
mv unfiltered_vcf.gz rmthis.vcf.gz
mv tmp.vcf.gz unfiltered_vcf.gz
# Double check that this file was generated correctly before removing the previous .vcf


# Index the vcf
tabix unfiltered_vcf.gz


