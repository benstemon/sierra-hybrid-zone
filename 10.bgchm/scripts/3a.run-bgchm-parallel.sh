#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=bgchm-HI-subset
#SBATCH --mail-type=END


cd $SLURM_SUBMIT_DIR

# Source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


##################
##### SETUP ######
##################
##########################################################################################

# Activate conda environment with R and bgchm installed
conda activate r-env

# Hard-coded input and output file names
# These files are generated with previous scripts
input1="genotypes-recoded-hybrids.tsv"
input2="genotypes-recoded-p0.tsv"
input3="genotypes-recoded-p1.tsv"
output1="1000snp-recoded-hybrids.tsv"
output2="1000snp-recoded-p0.tsv"
output3="1000snp-recoded-p1.tsv"

##########################################################################################

########################
##### SUBSET SNPS ######
########################
##########################################################################################


# Temporary file to store the header and random column indices
tempfile=$(mktemp)

# Get the number of columns in the file
num_cols=$(awk -F'\t' 'NR==1{print NF}' "$input1")

# Generate a list of column indices (1-based) and shuffle them, then pick 1000 columns
shuf -i 2-$num_cols -n 1000 | sort -n > "$tempfile"

# Prepend 1 (to always include the first column)
sed -i '1i1' "$tempfile"

# Create a comma-separated list of columns to select (1-based index)
columns=$(paste -sd, "$tempfile")

# Function to select columns and output to a new file
select_columns() {
    input=$1
    output=$2
    columns=$3
    awk -v cols="$columns" 'BEGIN {
        split(cols, a, ",");
        for (i in a) col_idx[a[i]];
    }
    {
        out = "";
        for (i=1; i<=NF; i++) {
            if (i in col_idx) {
                out = out (out == "" ? "" : "\t") $i;
            }
        }
        print out;
    }' FS='\t' OFS='\t' "$input" > "$output"
}

# Select the columns from each input file and write to the corresponding output file
select_columns "$input1" "$output1" "$columns"
select_columns "$input2" "$output2" "$columns"
select_columns "$input3" "$output3" "$columns"

# Clean up
rm "$tempfile"
##########################################################################################


# Submit the R script to estimate cline SDs and HI
Rscript --vanilla 3b.run-bgchm-parallel.R


