# Estimate DFR-like, DFR, and ANR gene tree

#### 1. Code for estimating the gene tree from unaligned sequences
All unaligned protein sequences putatively orthologous to genes of interest included in [data/unaligned-subsets](data/unaligned-subsets)
<br>

Then, generate combined unaligned file.
```
# Generate combined unaligned .fa with all protein sequences
cd data/unaligned-subsets
cat *.fa > ../unaligned-allsamples.fa
```
<br>

Next comes alignment with muscle
```
cd ..
mkdir genetree
muscle -in unaligned-allsamples.fa -out genetree/aligned-allsamples.fa
```
<br>

Manually change gene names in aligned-allsamples.fa to reduce clutter.
<br>

Next, use trimal to trim alignments
```
trimal -in genetree/aligned-allsamples.fa -out genetree/aligned-trimmed.fa -automated1
```
<br>

Finally, estimate gene tree, with Arabidopsis as outgroup, 1000 ultrafast bootstraps, and 1000 SH-aLRT test reps
```
iqtree -s genetree/aligned-trimmed.fa -m MFP -B 1000 -alrt 1000 --seqtype AA -o cinnamoyl_coa_reductase_Arabidopsis_thaliana
```