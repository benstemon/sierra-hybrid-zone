# Pollinator Visual Modeling

Additionally uses data from:
[1.phenotype_curation](~/1.phenotype_curation/)
[6.pigments](~/6.pigments/)
[10.bgchm](~/10.bgchm/)
raw floral and leaf reflectance spec data (obtain from supplemental data repo)

#### 1. Run visual models and plot colorspace model
* see [`1.run-visual-models.r`](scripts/1.run-visual-models.r)
* Runs the model on one set of parameters (The set used in final publication), but is set up such that a list of parameter values to try can be submitted and all combinations of those parameters will be run.
* Uses raw spec data for flowers and leaf tissue: these will need to be downloaded from the supplemental data repo

#### 2. Plotting spectra and visual model results
* see [`2.plotting-spec-and-visual-models.r`](scripts/2.plotting-spec-and-visual-models.r)