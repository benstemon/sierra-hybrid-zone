bcftools query \
  --regions-file DFRlike_GEMMAhits.txt \
  -f '%CHROM\t%POS[\t%GT]\n' \
  gemmasnps-.05missing.imputed-NOPARENTS.vcf.gz \
  > ANR_genotypes.txt


bcftools query -l gemmasnps-.05missing.imputed-NOPARENTS.vcf.gz > sample_names.txt


# Step 1: Get sample names into one line with tabs
header="site\t$(paste -sd '\t' sample_names.txt)"
body=$(awk '{print $1"_"$2"\t"substr($0, index($0, $3))}' DFRlike_genotypes.txt)
(echo -e "$header"; echo -e "$body") > final_DFRlike_genotype_matrix.tsv