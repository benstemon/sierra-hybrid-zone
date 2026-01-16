#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=perm_array
#SBATCH --array=0-99
#SBATCH --export=ALL,TRAIT_IDX=0
#SBATCH --output=logs/perm_array.%A_%a_trait${TRAIT_IDX}.out

cd $SLURM_SUBMIT_DIR

# Setup environment
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate gwas

# Paths
gemma="gemma-0.98.5"
bimbam="../gemmasnps-.05missing-NOPARENTS.bimbam"
GTmatrix="../GTMATRIX-NOPARENT-LDpruned.cXX.txt"

# Get phenotype file from list
mapfile -t PHENOFILES < pheno_file_list.txt
phenotype_file="${PHENOFILES[$TRAIT_IDX]}"

# Extract phenotype name
filename=$(basename "$phenotype_file")
phenoname="${filename#GEMMA_INPUTFILE.}"
phenoname="${phenoname%.txt}"

# Make scratch directory
scratch_dir="scratch/temp_${phenoname}_$SLURM_ARRAY_TASK_ID"
mkdir -p $scratch_dir

# Extract and shuffle phenotype column
cut -f2 "$phenotype_file" > "$scratch_dir/pheno_col.txt"
cut -f1 "$phenotype_file" > "$scratch_dir/ids.txt"
shuf "$scratch_dir/pheno_col.txt" > "$scratch_dir/shuf_values.txt"
paste "$scratch_dir/ids.txt" "$scratch_dir/shuf_values.txt" > "$scratch_dir/shuffled_pheno.txt"
shuffled="$scratch_dir/shuffled_pheno.txt"

# Run GEMMA
$gemma -g $bimbam \
    -p $shuffled -n 2 \
    -k $GTmatrix \
    -lmm 1 \
    -outdir $scratch_dir \
    -o temp_$SLURM_ARRAY_TASK_ID > /dev/null 2>&1

# Extract top 1% p-values
assoc_file="$scratch_dir/temp_${SLURM_ARRAY_TASK_ID}.assoc.txt"
outfile="perm_output/pval_top1pct_${phenoname}_${SLURM_ARRAY_TASK_ID}.txt"



# Filter out header and nan values, then count usable SNPs
# Need to filter nan because there are some edge cases causing (I think) numerical instability or model fitting issues -- differs across permutations so unsure of the reson why. When ordering by size they get pushed to the top.
if [ -f "$assoc_file" ]; then
    num_top=$(awk 'NR > 1 && $12 != "nan" {count++} END {printf "%.0f", count * 0.01}' "$assoc_file")
    
    # Extract valid p-values, sort, and write top 1%
    awk 'NR > 1 && $12 != "nan" {print $12}' "$assoc_file" | sort -g | head -n "$num_top" > "$outfile"
else
    echo "NA" > "$outfile"
fi


# Cleanup scratch
rm -r $scratch_dir
