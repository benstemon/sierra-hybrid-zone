#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=mapping_pipeline

#NOTE: this is an arrayed batch script.
#With 358 total samples...
#sbatch --array=0-357%7 2.mapping_pipeline_array.sh


cd $SLURM_SUBMIT_DIR


#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with packages installed
#needs samtools v1.15.1 and bamutil v1.0.15 (I also installed bwa but that should be OK from module)
#version of samtools installed on cluster is old (no markdup and old fixmate options)

conda activate mapping_etc


# Hard path to the filtered reads and reference genome
# Ensure final "/" at the end of basedir
basedir="sierra_hybrids/"
filtered_reads="filtered_reads_sierras"
refgenome="annot_Pdavidsonii_genome.1mb.fasta"
numthreads=4 # Change this to match number of cores allocated


# Path to the out directories for filtered bams and summary stats.
outdir=${basedir}mapped_filtered_bams_sierras
mkdir -p $outdir

statsdir=${basedir}sumstats_mapped_filtered_bams_sierras
mkdir -p $statsdir




# Identify target files
cd $filtered_reads

r1s=(*R1_001.fastq.gz)
read1="${r1s[$SLURM_ARRAY_TASK_ID]}"
read2="${read1/L001_R1/L001_R2}"


# Mapping and cleaning pipeline:
#a. map with bwa mem
#b. fixmate -m fills in mate coords and add score tags to tell markdup which reads to keep
#c. sort alignment (put mapped reads in physical order)
#d. markdup marks and removes (-r) duplicate reads
#e. view filters to retain only reads with mapping quality >20 (99% mapping confidence)
#f. bam clipOverlap softclips read overlaps

bwa mem -t $numthreads -M $refgenome $read1 $read2 | \
 samtools fixmate -@ $numthreads -m -u -O bam - - | \
 samtools sort -@ $numthreads -u | \
 samtools markdup -r -@ $numthreads -u -s - - | \
 samtools view -@ $numthreads -h -q 20 -u | \
 bam clipOverlap --in -.ubam \
 --out $outdir/"${read1/trimmed_L001_R1_001.fastq.gz/mapped_filtered.bam}" --stats

# Index reads
samtools index -b $outdir/"${read1/trimmed_L001_R1_001.fastq.gz/mapped_filtered.bam}"


