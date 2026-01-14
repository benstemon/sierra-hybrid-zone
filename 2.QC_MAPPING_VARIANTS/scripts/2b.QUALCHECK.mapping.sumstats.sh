#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=bam_sumstats


cd $SLURM_SUBMIT_DIR


# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


# Activate conda environment with samtools installed
conda activate mapping_etc


# Hard paths to reference genome, input, and output directories
refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
bamdir="/work/bs66/sierra_hybrids/mapped_filtered_bams_sierras"
outdir="/work/bs66/sierra_hybrids/sumstats_mapped_filtered_bams_sierras"


# .bam file summary stats
cd $bamdir
for i in *.bam;
do
    # Generate coverage file
    samtools coverage -m -A -w 40 $i > "$outdir/${i/_mapped_filtered.bam/.coverage.txt}"
    
    # Calculate mean depth of coverage
    samtools depth $i | \
     awk '{sum+=$3} END {print sum/NR}' > "$outdir/${i/_mapped_filtered.bam/.meandepth.txt}"
done


# Make genome-wide summaries of coverage plots
cd $outdir
echo -e "File\tCoverage" > PERCENT_COVERAGE_SUMMARY.txt
echo -e "File\tCoverage" > MEAN_DEPTH_SUMMARY.txt
for i in *coverage.txt;
do
    # Percent coverage per sample, per chromosome 
    printf "$i\t" >> PERCENT_COVERAGE_SUMMARY.txt
    awk '/Percent covered:/ {printf "%s\t", $NF} END {printf "\n"}' $i >> PERCENT_COVERAGE_SUMMARY.txt
    
    # Mean depth of coverage per sample, per chromosome
    printf "$i\t" >> MEAN_DEPTH_SUMMARY.txt
    awk '/Mean coverage:/ {printf "%s\t", $NF} END {printf "\n"}' $i >> MEAN_DEPTH_SUMMARY.txt
done

