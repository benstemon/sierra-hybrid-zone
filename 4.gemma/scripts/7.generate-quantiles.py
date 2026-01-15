import glob
import numpy as np
import os
from collections import defaultdict

# Create output directory
os.makedirs("quantiles_pvals", exist_ok=True)

# Gather all permutation result files
files = glob.glob("perm_output/pval_top1pct_*.txt")

# Collect p-values by phenotype
phenotype_pvals = defaultdict(list)

# Extract phenotype name and for each phenotype assemble p values
for fname in files:
    base = os.path.basename(fname)
    parts = base.replace("pval_top1pct_", "").replace(".txt", "").rsplit("_", 1)
    if len(parts) != 2:
        continue
    phenotype, _ = parts
    with open(fname) as f:
        for line in f:
            try:
                p = float(line.strip())
                phenotype_pvals[phenotype].append(p)
            except ValueError:
                continue

# Process each phenotype
for phenotype, pvals in phenotype_pvals.items():
    if not pvals:
        continue
    quantiles = np.percentile(pvals, [0.0001, 0.001, 0.01, 0.1, 1, 2.5, 50, 97.5])
    top_snp_pval = min(pvals)

    # Write output
    outfile = f"quantiles_pvals/{phenotype}_significance_quantiles.txt"
    with open(outfile, "w") as f:
        f.write("quantile\tvalue\n")
        f.write(f"0.0001%\t{quantiles[0]}\n")
        f.write(f"0.001%\t{quantiles[1]}\n")
        f.write(f"0.01%\t{quantiles[2]}\n")
        f.write(f"0.1%\t{quantiles[3]}\n")
        f.write(f"1%\t{quantiles[4]}\n")
        f.write(f"2.5%\t{quantiles[5]}\n")
        f.write(f"50%\t{quantiles[6]}\n")
        f.write(f"97.5%\t{quantiles[7]}\n")
        f.write(f"Top SNP\t{top_snp_pval}\n")
