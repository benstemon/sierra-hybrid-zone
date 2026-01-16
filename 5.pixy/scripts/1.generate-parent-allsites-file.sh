#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=variants_allsites


cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with bcftools installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

module load vcftools

# Information for important files and parameters
refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
bamdir="/work/bs66/sierra_hybrids/mapped_filtered_bams_sierras"
newparents="newberryi-fst-pops.txt"
davparents="davidsonii-fst-pops.txt"
numthreads=8



# A. Generate list of bamfiles to include for variant calling
# Clear the output file
> parent_bamlist.txt

# Iterate through both input files
for input_file in $newparents $davparents; do
    while read -r sample; do
        # Find the matching file in $bamdir
        match=$(ls "$bamdir" | grep "^${sample}.*\.bam$")
        if [[ -n "$match" ]]; then
            for file in $match; do
                # Write the full path of the matched file to the output file
                echo "$bamdir/$file" >> parent_bamlist.txt
            done
        else
            echo "No match found for $sample" >> parent_bamlist.txt
        fi
    done < "$input_file"
done

# Make new file for renaming .vcf sample names
cat $newparents $davparents > newnames.txt


# B. Make all-sites .vcf file with the created parent_bamlist.txt, otherwise the same
bcftools mpileup --threads $numthreads -Ou \
 -a FORMAT/AD,FORMAT/DP \
 --min-BQ 20 --min-MQ 30 \
 -f $refgenome -b parent_bamlist.txt | 
 bcftools call --threads $numthreads -m -f GQ,GP -Ou |
 bcftools norm --threads $numthreads -f $refgenome -Oz -o unfiltered_vcf.gz

# Change sample names in the header, to truncated version
bcftools reheader -s newnames.txt --threads $numthreads unfiltered_vcf.gz -o tmp.vcf.gz
mv unfiltered_vcf.gz rmthis.vcf.gz
mv tmp.vcf.gz unfiltered_vcf.gz
# Double check that this file was generated correctly before removing the previous .vcf


# Index the vcf
tabix unfiltered_vcf.gz



# C. filter Variant and Invariant sites in vcf
# Invariant sites -- this is easy
 vcftools --gzvcf unfiltered_vcf.gz \
 --max-non-ref-ac 0 \
 --recode --recode-INFO-all --stdout | 
 bcftools view - \
 --threads $numthreads \
 -Oz -o tmp-invariants.filtered.vcf.gz

# Variant sites
# Filters on QUALITY, DEPTH (~mean/2, mean*2), MAF, and MISSINGNESS of reads.
# In addition, filter clusters of indels and SNPs within 10 bp of indels
# and exclude all monomorphic SNPs (only variant alleles called)
vcftools --gzvcf unfiltered_vcf.gz \
 --non-ref-ac 1 \
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
 -Oz -o tmp-variants.filtered.vcf.gz


#index both vcfs
tabix tmp-invariants.filtered.vcf.gz
tabix tmp-variants.filtered.vcf.gz

#concatenate the two vcfs
bcftools concat tmp-invariants.filtered.vcf.gz \
 tmp-variants.filtered.vcf.gz \
 --threads $numthreads --allow-overlaps -Oz \
 -o filtered_allsites_PARENTS.vcf.gz

#remove temporary files
rm tmp-invariants.filtered.vcf.gz*
rm basefilter.vcf.gz*
rm tmp-variants.filtered.vcf.gz*

#index final merged file
tabix filtered_allsites_PARENTS.vcf.gz


