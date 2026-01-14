#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=vcf_sumstats

echo "This script uses CLI with three parameters: input vcf, output directory, and number of snps to subsample a .vcf and extract a series of summary statistics."

cd $SLURM_SUBMIT_DIR

##################
##### SETUP #####
##################

# Check if all required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <invcf> <outdir> <nsnp>"
    exit 1
fi

invcf="$1"
outdir="$2"
nsnp="$3"

########################
# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Activate conda environment with bcftools installed
conda activate mapping_etc

# Load other modules
module load vcftools

# Create output directory if it doesn't exist
mkdir -p $outdir

# Create the .vcf index if it doesn't exist
if [ ! -f "$invcf.tbi" ]; then
    tabix "$invcf"
fi

#####################
### Subsample VCF ###
#####################

# Check how many unfiltered variants
bcftools +counts $invcf > variant-counts-"${invcf%%.vcf.gz}".txt

# Subsample VCF of interest to specified number of random variants
# This samples with replacement but is fast
subsetname="subset-${invcf%%.vcf.gz}-${nsnp}snps.vcf"
bcftools view $invcf --no-header --exclude-types indels | shuf -n $nsnp > $subsetname

# Append header back to the subsampled VCF
tabix -H $invcf > tmpfile
cat tmpfile $subsetname > tmpfile2
mv tmpfile2 $subsetname
rm tmpfile

#####################
### Summary stats ###
#####################

outprefix="$outdir/stats.${invcf%%.vcf.gz}"

# Allele frequency
vcftools --vcf $subsetname --freq2 --out $outprefix --max-alleles 2

# Mean depth/individual
vcftools --vcf $subsetname --depth --out $outprefix

# Mean depth/site
vcftools --vcf $subsetname --site-mean-depth --out $outprefix

# Site quality
vcftools --vcf $subsetname --site-quality --out $outprefix

# Missing data/individual
vcftools --vcf $subsetname --missing-indv --out $outprefix

# Missing data/site
vcftools --vcf $subsetname --missing-site --out $outprefix

# Heterozygosity and inbreeding coefficient/individual
vcftools --vcf $subsetname --het --out $outprefix

# Per site heterozygosity and homozygosity
vcftools --vcf $subsetname --hardy --out $outprefix

echo "Script execution completed."
