#!/bin/sh

#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p wessinger-48core
#SBATCH --job-name=filter-variants


cd $SLURM_SUBMIT_DIR


##################
##### SETUP #####
##################

########################
# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

# Activate conda environment with bcftools installed
conda activate mapping_etc

# Load other modules
module load vcftools

# Hard path to the unfiltered .vcf file
unfiltered_vcf="/work/bs66/sierra_hybrids/VCFs-sierras/unfiltered_vcf.gz"

# Change this to match number of cores allocated
numthreads=12

########################


##################
### filter VCF ###
##################


# Filters on QUALITY, DEPTH (~mean/2, mean*2), MAF, and MISSINGNESS of reads.
# In addition, filter clusters of indels and SNPs within 10 bp of indels
# and exclude all monomorphic SNPs (only variant alleles called)
vcftools --gzvcf $unfiltered_vcf \
 --min-alleles 2 \
 --min-meanDP 2.5 \
 --max-meanDP 10 \
 --minDP 2 \
 --max-missing 0.50 \
 --maf 0.05 \
 --recode --recode-INFO-all --stdout | 
 bcftools view - -m2 -M2 -v snps -Ou |
 bcftools filter - \
 --IndelGap 10 \
 --SnpGap 10 \
 --soft-filter LowQual \
 --exclude 'AC==AN || %QUAL<20' \
 --set-GTs . \
 --threads $numthreads \
 -Oz -o basefilter.vcf.gz

# Additional filter on overly heterozygous sites.
bcftools query -f '%CHROM\t%POS\t%INFO/DP[\t%GT]\n' basefilter.vcf.gz |
 awk -F '\t' '{
   hom_ref=0;
   het=0;
   hom_alt=0;
   for (i=4; i<=NF; i++) {
     if ($i == "0/0" || $i == "0|0") hom_ref++;
     else if ($i == "0/1" || $i == "1/0" || $i == "0|1" || $i == "1|0") het++;
     else if ($i == "1/1" || $i == "1|1") hom_alt++;
   }
   total=hom_ref + het + hom_alt;
   if (total > 0) {
     heterozygosity = (het / total);
     if (heterozygosity >= 0.75) print $1 "\t" $2;
   }
 }' > het.75-sites-removed.txt


bcftools view basefilter.vcf.gz -T ^het.75-sites-removed.txt \
 -Oz -o finalfilter.vcf.gz


# Index the final filtered .vcf. The basefilter .vcf can be removed if desired.
#rm basefilter.vcf.gz
tabix finalfilter.vcf.gz




