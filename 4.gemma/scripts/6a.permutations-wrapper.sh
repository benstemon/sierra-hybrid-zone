#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=wrapper


# Wrapper to submit one SLURM array job per trait

num_traits=$(wc -l < pheno_file_list.txt)

for i in $(seq 0 $((num_traits - 1))); do
    sbatch --export=ALL,TRAIT_IDX=$i \
        --output=logs/perm_array.%A_%a_trait${i}.out \
        1b.permutations-alltraits.sh
done
