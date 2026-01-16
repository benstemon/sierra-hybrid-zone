# Genomic signatures of reproductive isolation are decoupled from floral divergence in a long-standing hybrid zone

This code corresponds to Stone et al. 2026: "Genomic signatures of reproductive isolation are decoupled from floral divergence in a long-standing hybrid zone". [`DOI: arxiv doi`](link to paper)

The link to supplemental data on FigShare is here: [`name`](link to data repo)

The raw reads associated with this study are stored in the SRA repository under project number PRJNAxxxxx: [`url to SRA data`](link to SRA data)

Note that this study uses the *Penstemon davidsonii* genome as a reference for mapping. Information about that genome can be found in [`Ostevik et al. 2024`](https://academic.oup.com/g3journal/article/14/3/jkad296/7503391), and the genome and associated materials are stored in the SRA here: [`https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1010203`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1010203)

## Directory setup
Directories are generally in order of analysis pipeline, though some cross-contamination cannot be helped

### [Phenotype curation](1.phenotype_curation/)
Scripts for organizing and analyzing phenotypic variation

### [Quality control](2.QC_MAPPING_VARIANTS/)
Scripts for assessing read quality and filtering raw reads, then mapping reads, filtering, and assessing mapping quality

### [Admixture and Genomic PCA](3.admixture-pca/)
Scripts for running admixture and a genomic PCA and some basic visualization

### [GWAS](4.gemma/)
Scripts for imputing genotypes, performing GWAS, and plotting results (among other related things)
