# Genomic signatures of reproductive isolation are decoupled from floral divergence in a long-standing hybrid zone

This code corresponds to Stone et al. 2026: "Genomic signatures of reproductive isolation are decoupled from floral divergence in a long-standing hybrid zone". [`DOI: arxiv doi`](link to paper)

The link to supplemental data on FigShare is here: [`name`](link to data repo)

The raw reads associated with this study are stored in the SRA repository under project number PRJNAxxxxx: [`url to SRA data`](link to SRA data)

Note that this study uses the *Penstemon davidsonii* genome as a reference for mapping. Information about that genome can be found in [`Ostevik et al. 2024`](https://academic.oup.com/g3journal/article/14/3/jkad296/7503391), and the genome and associated materials are stored in the SRA here: [`https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1010203`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1010203)

## Directory setup
Directories are generally in order of analysis pipeline, though some cross-contamination cannot be helped. Most directories also include scripts for data visualization, mostly because in many cases visualization is integrated with analysis.

### [Phenotype curation](1.phenotype_curation/)
Scripts for organizing and analyzing phenotypic variation

### [Quality control](2.QC_MAPPING_VARIANTS/)
Scripts for assessing read quality and filtering raw reads, then mapping reads, filtering, and assessing mapping quality

### [Admixture and Genomic PCA](3.admixture-pca/)
Scripts for running admixture and a genomic PCA

### [GWAS](4.gemma/)
Scripts for imputing genotypes and performing GWAS (among other related things)

### [Genomic divergence and differentiation](5.pixy)
Scripts for calculating dxy, Fst, pi and others

### [Pigment analysis](6.pigments)
Scripts for assessing relationships between genotype and pigment production

### [Local ancestry inference](7.ancestrHMM)
Scripts for inferring local ancestry

### [Heterozygosity](8.heterozygosity)
Scripts for estimating heterozygosity and heterozygosity excess

### [Pollinator Visual Modeling](9.visualmodels)
Scripts for conducting pollinator visual modeling for honeybees and humnmingbirds

### [Bayesian Genomic Clines](10.bgchm)
Scripts for estimating bayesian genomic clines with bgchm

