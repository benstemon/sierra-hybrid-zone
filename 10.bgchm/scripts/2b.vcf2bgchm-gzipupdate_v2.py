import pandas as pd
import argparse
import gzip
from io import StringIO

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process input files and output files for bgchm.")
parser.add_argument("--vcf", type=str, help="Path to the input .vcf file")
parser.add_argument("--samples", type=str, help="Path to the sample list file")
parser.add_argument("--pure_parents", type=str, help="Path to the pure parent list file")
parser.add_argument("--p0list", type=str, help="Comma-separated list of unique identifiers for P0 samples. Choose either this option or p0file.")
parser.add_argument("--p1list", type=str, help="Comma-separated list of unique identifiers for P1 samples. Choose either this option or p1file.")
parser.add_argument("--p0file", type=str, help="Path to a file containing P0 sample names. Choose either this option or p0list.")
parser.add_argument("--p1file", type=str, help="Path to a file containing P1 sample names. Choose either this option or p1list.")
args = parser.parse_args()

# Function to read sample names from a file
def read_sample_file(file_path):
    with open(file_path, 'r') as f:
        return [line.strip() for line in f.readlines() if line.strip()]

# Function to read a gzipped VCF file
def read_vcf_gz(file_path):
    with gzip.open(file_path, 'rt') as f:
        lines = [line for line in f if not line.startswith('##')]
    vcf = pd.read_csv(StringIO(''.join(lines)), sep='\t', comment='#', header=None)
    return vcf

# Read input files
vcf = read_vcf_gz(args.vcf)
samples = pd.read_csv(args.samples, header=None)
pure_parents = pd.read_csv(args.pure_parents, header=None)

# Convert p0 and p1 lists/files to lists of sample names
if args.p0file:
    p0list = read_sample_file(args.p0file)
elif args.p0list:
    p0list = args.p0list.split(",")
else:
    p0list = []

if args.p1file:
    p1list = read_sample_file(args.p1file)
elif args.p1list:
    p1list = args.p1list.split(",")
else:
    p1list = []




# Extract genotype data and transpose
genos = vcf[vcf.columns[9:]]
genos.columns = samples[0]
for c in genos.columns:
    for i in genos.index:
        genos.at[i, c] = genos.at[i, c].split(":")[0]
genos_T = genos.T

# Create names for columns
nc = []
for i in vcf.index:
    chr = vcf.loc[i][0]
    ps = str(vcf.loc[i][1])
    variant = chr + ":" + ps
    nc.append(variant)
genos_T.columns = nc

# Write all samples to file
genos_T.to_csv("genotypes-uncoded-allsamples.tsv", sep="\t")

# Recode genotypes based on bgchm specifications
genos_coded_T = genos_T.replace("0|0", 0)
genos_coded_T = genos_coded_T.replace("0|1", 1)
genos_coded_T = genos_coded_T.replace("1|0", 1)
genos_coded_T = genos_coded_T.replace("1|1", 2)
genos_coded_T = genos_coded_T.replace(".|.", None)

# Write recoded genotypes (all samples) to file
genos_coded_T.to_csv("genotypes-recoded-allsamples.tsv", sep="\t")

# Generate the sub-matrices with parental reference sets
p0 = []
p1 = []

for s in pure_parents[0]:
    if any(pop in s for pop in p0list):
        p0.append(s)
    elif any(pop in s for pop in p1list):
        p1.append(s)


# Generate and write P1 recoded genotypes
genos_coded_T.loc[p0].to_csv("genotypes-recoded-p0.tsv", sep="\t")

# Generate and write P2 recoded genotypes
genos_coded_T.loc[p1].to_csv("genotypes-recoded-p1.tsv", sep="\t")

# Generate and write hybrid recoded genotypes
genos_coded_T.drop(pure_parents[0]).to_csv("genotypes-recoded-hybrids.tsv", sep="\t")
