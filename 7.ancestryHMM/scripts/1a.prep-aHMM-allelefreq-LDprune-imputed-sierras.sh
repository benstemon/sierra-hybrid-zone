#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=prep_aHMM

# This script prepares the input file for ancestry HMM

cd $SLURM_SUBMIT_DIR

##################
##### SETUP #####
##################


# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Activate conda environment with bcftools installed
conda activate mapping_etc

# Path to the input .vcf and the reference genome index file
invcf="gemmasnps-.05missing.imputed.vcf.gz"
reference_fai="annot_Pdavidsonii_genome.1mb.fasta.fai"

# Set the SNP window size, window slide, and r^2 filtering threshold
# In admixture I used 50 SNP windows, sliding 10 SNPs, and r^2 threshold of 0.1
# Here I will consider 50kb windows, and r2 threshold of 0.2
windowsize=50
r2thresh=0.2

# Path to the parent pop index for the python script to find SNPs with allele freq diffs
# REMEMBER THE INDEX IS BASE-0!
parentpops="data/parentpop-index-sierras.txt"

# Set the shared prefix for this run
pref="0.75-sierras"

# Path to the ancestry HMM samples file
ahmm_samples="data/aHMM_samples-sierras.txt"


###########################
##### FILE GENERATION #####
###########################

# Generate header file
bcftools view -h $invcf > tmp.txt
awk 'FNR==NR {contig_info = contig_info "##contig=<ID=" $1 ",length=" $2 ">\n"; next} /^##INFO/ && !printed {gsub(/\n$/, "", contig_info); print contig_info; printed=1} {print}' $reference_fai tmp.txt > imputed_vcf_header.txt
rm tmp.txt

# Make a new output directory
mkdir -p infiles_$pref

# Filter down to ancestry-informative sites
python filter_allele_diff_snps.py -f $parentpops -i $invcf -s infiles_$pref/outsnps-$pref.txt -o infiles_$pref/$pref-snps.vcf.gz -m imputed_vcf_header.txt -a 0.75 --p1count 10 --p2count 10

# Check how many variants this .vcf has
bcftools +counts infiles_$pref/$pref-snps.vcf.gz > infiles_$pref/variant-counts-$pref-snps.txt

# Filter out to two additional .vcfs, one for each of the parent populations
#NEWBERRYI -- coded as 1 in the file
awk '$2 == 1 {print $1}' $parentpops > infiles_$pref/newberryi_list.txt
bcftools view infiles_$pref/$pref-snps.vcf.gz -S infiles_$pref/newberryi_list.txt -Oz -o infiles_$pref/$pref-snps-newberryi.vcf.gz

#DAVIDSONII -- coded as 0 in the file
awk '$2 == 0 {print $1}' $parentpops > infiles_$pref/davidsonii_list.txt
bcftools view infiles_$pref/$pref-snps.vcf.gz -S infiles_$pref/davidsonii_list.txt -Oz -o infiles_$pref/$pref-snps-davidsonii.vcf.gz


###########################
##### LD PRUNING (R2) #####
###########################


# Switch to python directory with plink2 installed
conda activate admixture

# Convert each of the .vcf files to plink .bed format
#ALLSAMPLES
plink2 --vcf infiles_$pref/$pref-snps.vcf.gz --make-bed --out infiles_$pref/allsamples-$pref --double-id --allow-extra-chr --set-missing-var-ids @:#[RDH]\$1,\$2
#NEWBERRYI
plink2 --vcf infiles_$pref/$pref-snps-newberryi.vcf.gz --make-bed --out infiles_$pref/newberryi-$pref --double-id --allow-extra-chr --set-missing-var-ids @:#[RDH]\$1,\$2
#DAVIDSONII
plink2 --vcf infiles_$pref/$pref-snps-davidsonii.vcf.gz --make-bed --out infiles_$pref/davidsonii-$pref --double-id --allow-extra-chr --set-missing-var-ids @:#[RDH]\$1,\$2


# Prune data in LD to get SNP lists
#NEWBERRYI
plink2 --bfile infiles_$pref/newberryi-$pref --indep-pairwise $windowsize kb $r2thresh --out infiles_$pref/newberryi-$pref --allow-extra-chr --bad-ld
#DAVIDSONII
plink2 --bfile infiles_$pref/davidsonii-$pref --indep-pairwise $windowsize kb $r2thresh --out infiles_$pref/davidsonii-$pref --allow-extra-chr --bad-ld


# Combine unique retained sites into a single file
cat infiles_$pref/newberryi-$pref.prune.in infiles_$pref/davidsonii-$pref.prune.in | sort | uniq > infiles_$pref/newberryi.davidsonii-$pref.uniq.prune.in


# Extract SNPs from the full file
plink2 --bfile infiles_$pref/allsamples-$pref --extract infiles_$pref/newberryi.davidsonii-$pref.uniq.prune.in --make-bed --out infiles_$pref/allsamples-$pref-pruned --allow-extra-chr


# Convert back to a .vcf
plink2 --bfile infiles_$pref/allsamples-$pref-pruned --recode vcf --out infiles_$pref/allsamples-$pref-pruned --allow-extra-chr

# Check how many variants this final .vcf has
conda deactivate
bcftools +counts infiles_$pref/allsamples-$pref-pruned.vcf > infiles_$pref/variant-counts-$pref-pruned-snps.txt

# reheader the vcf to fix plink weirdness in output
bcftools reheader -h imputed_vcf_header.txt -o TMP.vcf infiles_$pref/allsamples-$pref-pruned.vcf

# replace the original pruned vcf with the reheadered version
mv TMP.vcf infiles_$pref/allsamples-$pref-pruned.vcf


###########################
##### aHMM CONVERSION #####
###########################


# Convert the filtered .vcf file to ancestry HMM format
# -g 1 means to use the genotype field rather than allele depths
# -m 0 means no filter on distance between SNPs
# --min_diff 0 means no threshold for allele freq diffs between parents (we already filtered on that)
python vcf2ahmm-BWSedits.py -g 1 -m 0 --min_diff 0 -v infiles_$pref/allsamples-$pref-pruned.vcf -s $ahmm_samples > infiles_$pref/allsamples-$pref-pruned.ahmm

# rename the outputted ploidy file
mv ahmm.ploidy infiles_$pref/ahmm-$pref.ploidy

