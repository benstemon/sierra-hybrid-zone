## QC and mapping, deduplication, overlap clipping, and summary stats

### Quality control pipeline
1. Run QC of raw sequencing reads (fastqc) and summarize (multiqc)
2. Merge Illumina lanes by forward and reverse reads
3. 
    a. Trim adapters, quality filter, and enable base correction in overlapped regions (fastp)
    b. Run QC on trimmed, filtered data (fastqc)
    c. Summarize results (multiqc)


#### 1. QC on raw reads with fastqc, summarized with multiqc. Then trim and quality filter with fastp, and re-check QC.
* see [`1.run_QC.sh`](scripts/1.run_QC.sh)
* Default options for quality filtering
* Base correction for overlapping reads enabled
* poly-x trimming on 3' ends enabled
* Limit read length to 30 bp
* Enable auto-detection of adapters


#### 2a. Map and filter reads
* see [`2a.mapping_pipeline_array.sh`](scripts/2a.mapping_pipeline_array.sh)
* Duplicates marked and removed with samtools markdup
* Reads with low mapping quality (q<20)
* Overlapping paired end reads clipped with bamutil clipOverlap
* Reads are mapped to the *P. davidsonii* reference genome with bwa-mem


#### 2b. Additional QC to check summary stats of mapped reads
* see [`2b.QUALCHECK.mapping.sumstats`](scripts/2b.QUALCHECK.mapping.sumstats)
* Some basic summary statistics generated with samtools coverage and samtools stats

#### 3. Call variants
* see [`3.call_variants.sh`](3.call_variants.sh)
* Calls variants with bcftools
* Also see [`samplename_changer.txt`](data/samplename_changer.txt), which is a basic text file with info for name substitution.

#### 4. Filter variants
* see [`4.filter_variants.sh`](scripts/4.filter_variants.sh)
* Calls variants with bcftools

#### Bonus script
* see [`QUALCHECK.vcf.sumstats.sh`](scripts/QUALCHECK.vcf.sumstats.sh) for a simple way to get some basic summary stats from vcfs, including number of SNPs, depth/site, etc. while optionally subsetting the data.



