# Local Ancestry Inference

Additionally uses data from:
imputed .vcf
reference genome.fai file
[1.phenotype_curation](~/1.phenotype_curation)


#### 1. Filter to ancestry-informative sites, LD-prune, and convert to aHMM format
* see [`1a.prep-aHMM-allelefreq-LDprune-imputed-sierras.sh`](scripts/1a.prep-aHMM-allelefreq-LDprune-imputed-sierras.sh), [`1b.filter_allele_diff_snps.py`](scripts/1b.filter_allele_diff_snps.py), and [`1c.vcf2ahmm-BWSedits.py`](scripts/1c.vcf2ahmm-BWSedits.py)
* The first script [`1b.filter_allele_diff_snps.py`](scripts/1b.filter_allele_diff_snps.py) calls the other two (python scripts to convert genotypes from .vcf format into AncestryHMM format.
* The first python script filters to SNPs with a given AFD between parent pops. In this case, it is 75%
* The second python script [`1c.vcf2ahmm-BWSedits.py`](scripts/1c.vcf2ahmm-BWSedits.py) is from ancHMM github (found [here](https://github.com/russcd/Ancestry_HMM/blob/master/scripts/vcf2ahmm.py)), with some edits to work with this specific phased data format ("|")
* also uses plink2, bcftools, and ahmm

#### 2. Run AncestryHMM
* see [`2.run-ancestryHMM.sh`](scripts/2.run-ancestryHMM.sh)
* uses files generated from previous step

#### 3. Summarize posteriors and generate plots
* see [`3.plot_ancHMM.R`][scripts/3.plot_ancHMM.R]
* uses posterior results from step 2: can download these from the online supplemental data repo
* also uses [`combined phenotype dataset`](~/1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv)
