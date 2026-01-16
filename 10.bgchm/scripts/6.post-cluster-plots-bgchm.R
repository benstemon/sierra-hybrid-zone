# Plotting genomic clines stuff
# There's probably a lot more in here than needed but I tried so many things
# And it's all integrated, so, sorry about that

library(tidyverse)
library(ggh4x)
library(ggpmisc)
library(bgchm)
library(ggsignif)
library(scales)
library(zoo)

setwd("~/10.bgchm/")

# Data I/O -- connecting across analyses
################################################################################
# Read in parent names
parentnames <- read_delim("data/parentlist-bgchm-sierras.txt",
                          delim = "\t", col_names = F) %>%
  rename(newID = X1) %>%
  mutate(newID = str_remove(newID, "BWS_"))

# Read in sample names, removing parent samples
samplenames <- read_delim("data/samplelist-bgchm-sierras.txt",
                          delim = "\t", col_names = F) %>%
  rename(newID = X1) %>%
  mutate(newID = str_remove(newID, "BWS_")) %>%
  rowwise() %>%
  mutate(pop = as.numeric(str_split(newID, "-")[[1]][1])) %>%
  ungroup() %>%
  filter(!newID %in% parentnames$newID)

# Read in the HI data and link to sample names
HI_df <- read_delim("data/hybrid_index-sierras.txt", delim = "\t") %>%
  cbind(., samplenames) %>%
  arrange(.[, 1]) %>%
  mutate(index = row_number()) %>%
  rename(HI_50 = X50.,
         HI_2.5 = X2.5.,
         HI_5 = X5.,
         HI_95 = X95.,
         HI_97.5 = X97.5.)

# Next, load in admixture data and check relationship between admixture and HI
admixture_data <- read.csv("~/1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv")

combined_data <- left_join(HI_df, admixture_data, by = c("newID", "pop"))




# Read in SNP ID information -- will come in handy soon
snpids <- read_delim("data/subsampled-maf0.1-bgcinput.txt",
                     delim = "\t", col_names = F)

# Read in SNP AFD and MAF information
# get this from online repo
AFDs <- read_delim("AFD-table-parents.txt",
                   delim = " ") %>%
  select(chr, ps, AFD) %>%
  # get this from online repo supplemental data
  left_join(read_csv("data/MAF-snps-NOPARENTS.csv"),
            by = c("chr", "ps"))

# Read in gradients and centers
# Each row is a SNP
# all "result" files must be obtained from online data repo supplemental
gradients <- read_delim("result-cline-gradients.csv", delim = ",") %>%
  select(-1) %>%
  rename("gradient_50" = V1, "gradient_5" = V2, "gradient_95" = V3)

centers <- read_delim("result-cline-centers.csv", delim = ",") %>%
  select(-1) %>%
  rename("center_50" = V1, "center_5" = V2, "center_95" = V3)

# Combine to make all clines object, and merge with AFD + MAF
allclines <- cbind(gradients, centers, snpids) %>%
  mutate(snpnumber = row_number()) %>%
  rename("chr" = X1, "ps" = X2) %>%
  mutate(chr = as.numeric(str_remove(chr, "scaffold_"))) %>%
  left_join(., AFDs, by = c("chr", "ps"))

# Add trait information to clines
traitclines <- read_csv("data/data-fullinfo-traitsnps-maf0.1-2scaff0.1-bgchm.csv") %>%
  select(-bpstart, -bpend) %>%
  left_join(allclines, by = c("chr", "ps")) %>%
  mutate(snpclass = "trait") %>%
  select(-AFD.y, -maf.y) %>%
  rename(AFD = AFD.x, maf = maf.x)

backgroundclines <- anti_join(allclines, traitclines, by = c("chr", "ps")) %>%
  mutate(phenotype = NA,
         modelID = NA,
         snpclass = "background")

# Merge these to make new allclines object
# also -- add AFDbins
allclines <- bind_rows(traitclines, backgroundclines) %>%
  relocate(snpnumber, chr, ps) %>%
  group_by(snpnumber) %>%
  arrange(snpnumber) %>%
  ungroup() %>%
  mutate(AFDbin = case_when(
    AFD >= 0.00 & AFD < 0.05  ~ "0-05",
    AFD >= 0.05 & AFD < 0.10  ~ "05-10",
    AFD >= 0.10 & AFD < 0.15  ~ "10-15",
    AFD >= 0.15 & AFD < 0.20  ~ "15-20",
    AFD >= 0.20 & AFD < 0.25  ~ "20-25",
    AFD >= 0.25 & AFD < 0.30  ~ "25-30",
    AFD >= 0.30 & AFD < 0.35  ~ "30-35",
    AFD >= 0.35 & AFD < 0.40  ~ "35-40",
    AFD >= 0.40 & AFD < 0.45  ~ "40-45",
    AFD >= 0.45 & AFD < 0.50  ~ "45-50",
    AFD >= 0.50 & AFD < 0.55  ~ "50-55",
    AFD >= 0.55 & AFD < 0.60  ~ "55-60",
    AFD >= 0.60 & AFD < 0.65  ~ "60-65",
    AFD >= 0.65 & AFD < 0.70  ~ "65-70",
    AFD >= 0.70 & AFD < 0.75  ~ "70-75",
    AFD >= 0.75 & AFD < 0.80  ~ "75-80",
    AFD >= 0.80 & AFD < 0.85  ~ "80-85",
    AFD >= 0.85 & AFD < 0.90  ~ "85-90",
    AFD >= 0.90 & AFD < 0.95  ~ "90-95",
    AFD >= 0.95 & AFD <= 1.00 ~ "95-100"
  ))


# WRITE THE ALLCLINES DATA.FRAME
# This is also available from online data supplemental repo
write_delim(allclines, file = "result-allclines-sierras.txt", delim = "\t")


# Filter to steep clines (those with CIs fully > 1)
steepclines <- allclines %>%
  filter(gradient_5 > 1)
write_delim(steepclines, file = "result-steepclines-sierras.txt", delim = "\t")


# Filter to shallow clines
shallowclines <- allclines %>%
  filter(gradient_95 < 1)
write_delim(shallowclines, file = "result-shallowclines-sierras.txt", delim = "\t")


# Also write object with neutral clines
neutralclines <- allclines %>%
  anti_join(., shallowclines, by = c('chr', 'ps')) %>%
  anti_join(., steepclines, by = c('chr', 'ps'))
write_delim(neutralclines, file = "result-neutralclines-sierras.txt", delim = "\t")

################################################################################


# (1) Plot Hybrid index, and correlation between HI and admixture proportion
################################################################################
# Plot hybrid index estimates with 90% equal-tail probability intervals
ggplot(data = combined_data, aes(x = index, y = HI_50)) +
  geom_point(size = 0.4) +
  geom_segment(data = combined_data,
               aes(x = index, xend = index, y = HI_5, yend = HI_95),
               alpha = 0.5) +
  labs(x = "Individual (sorted by HI)", y = "Hybrid index (HI)") +
  theme_bw() +
  ggtitle("Sierras: Hybrid Index + 90% equal-tail probability intervals")
ggsave("plot-hybrid-index-sierras.png",
       height = 4, width = 4, dpi = 400, units = "in")


# Make regression of admixture proportion and hybrid index
ggplot(data = combined_data, aes(x = K1, y = HI_50)) +
  geom_point() +
  stat_poly_line() +
  stat_poly_eq(label.x = "right") +
  labs(x = bquote("Admixture Proportion " * italic("(P. newberryi)")),
       y = "Hybrid Index",
       title = "Sierras: Admixture proportion vs. Hybrid Index") +
  theme_bw() +
  theme(title = element_text(size = 9))
ggsave("plot-regression-hybridindex-admixtureproportion-sierras.png",
       height = 4, width = 4, dpi = 400, units = "in")

################################################################################


# (2) Identify steep and shallow clines
################################################################################
which.max(steepclines$gradient_50)
which.min(shallowclines$gradient_50)

# Some basic stats
mean(allclines$gradient_50)
max(steepclines$gradient_50)
min(shallowclines$gradient_50)
nrow(steepclines)
nrow(steepclines)/nrow(allclines)
nrow(shallowclines)
nrow(shallowclines)/nrow(allclines)

# Test plots...
# gencline_plot(center = c(0.5, 0.5), v = c(0.2,5), pdf = F, cvec = c("blue", "red"))

#gencline_plot(center = steepclines$center_50[5475], v = steepclines$gradient_50[5475],
#              pdf = F)
#gencline_plot(center = shallowclines$center_50[6345], v = shallowclines$gradient_50[6345],
#              pdf = F)


# Plot only the top n steepest
# Uncomment for specific clines, also need to re-comment some other parts
nsteep = 20
# New df, where there are 101 points (each percentage 0-100 of HI)
# and for each SNP, you have (a) HI, (b) constant u, and (c) phi
top_steepest_clines <- steepclines %>%
  arrange(desc(gradient_50)) %>%
  slice_head(n = nsteep) %>%
  slice(rep(row_number(), each = 101)) %>%
  mutate(h = rep(seq(0, 1, by = 0.01), times = nsteep)) %>%
  mutate(u = log(center_50/(1-center_50)) * gradient_50) %>%
  mutate(phi = (h^gradient_50)/(h^gradient_50 + (1-h)^gradient_50 * exp(u)))

ggplot(data = top_steepest_clines, aes(x = h, y = phi, group = snpnumber)) +
  geom_line() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = "dashed", col = "red") +
  scale_x_continuous(limits = c(0.0, 1.0), expand = c(0.005, 0.005)) +
  scale_y_continuous(limits = c(0.0, 1.0), expand = c(0.005, 0.005)) +
  labs(x = "Hybrid Index", y = "Ancestry Probability") +
  theme_bw() +
  ggtitle("Sierras: Top 20 steepest clines")
ggsave("plot-top20-steepest-clines-sierras.png",
       height = 4, width = 4, dpi = 400, units = "in")

# These are the locations of the top 20 snps
allclines[unique(top_steepest_clines$snpnumber),]

################################################################################


# (3) Rolling window steep-to-shallow ratio (8 main chromosomes)
# Right now includes AFD 0.3 filter, beware
################################################################################
windowsize = 500000
ratiomat <- plot_allclines %>%
  filter(chr %in% c("1086", "1087", "2532", "2533", "2684", "2685", "2686", "2687")) %>%
  mutate(centertype = case_when(center_95 < 0.45 ~ "dav",
                                center_5 > 0.55 ~ "new",
                                TRUE ~ "neutral")) %>%
  mutate(slopetype = case_when(gradient_95 < 1 ~ "shallow",
                               gradient_5 > 1 ~ "steep",
                               TRUE ~ "neutral")) %>%
  # filter(AFD > 0.3) %>% # See what happens if I filter out uninformative SNPs
  mutate(window = floor(ps / windowsize) * windowsize) %>%
  group_by(chr, window) %>%
  summarise(n_total = n(),
            dav = sum(centertype == "dav") / n_total,
            new = sum(centertype == "new") / n_total,
            neutral_center = sum(centertype == "neutral") / n_total,
            shallow = sum(slopetype == "shallow") / n_total,
            steep = sum(slopetype == "steep") / n_total,
            neutral_slope = sum(slopetype == "neutral") / n_total,
            .groups = "drop")

# Make new lenghts for chromosomes
chr_lengths <- allclines %>%
  filter(chr %in% ratiomat$chr) %>%
  group_by(chr) %>%
  summarise(chr_length = max(ps)) %>%
  arrange(chr) %>%
  mutate(offset = cumsum(lag(chr_length, default = 0)),
         chr_color = factor(row_number() %% 2))

# Make new labels for chromosomes
chr_labels <- chr_lengths %>%
  filter(chr %in% ratiomat$chr) %>%
  mutate(mid = offset + chr_length / 2)

chr_shades <- data.frame(chr = unique(ratiomat$chr)) %>%
  arrange(chr) %>%
  mutate(shade = rep(c("light", "dark"), length.out = n()))

# SLOPE MAT
# Refresh with cumulative chr lengths
# Then pivot longer and make a smoothed version
ratiomat_slope <- ratiomat %>%
  left_join(chr_lengths %>% select(chr, offset, chr_color), by = "chr") %>%
  mutate(cum_pos = window + offset) %>%
  pivot_longer(cols = c("neutral_slope", "shallow", "steep"),
               names_to = "slopetype", values_to = "value") %>%
  left_join(chr_shades, by = "chr") %>%
  mutate(slope_color_group = paste0(slopetype, "_", shade)) %>%
  arrange(chr, cum_pos) %>%
  group_by(chr, slopetype) %>%
  mutate(value_smooth = rollmean(value, k = 5, fill = NA, align = "center")) %>%
  ungroup()
slope_colors <- c(shallow_light = "#1f78b4",
                  shallow_dark  = "#1f78b4",
                  steep_light   = "#ff7f00",
                  steep_dark    = "#ff7f00",
                  neutral_slope_light = "#33a02c",
                  neutral_slope_dark  = "#33a02c" )

# CENTER MAT
# Refresh with cumulative chr lengths
# Then pivot longer and make a smoothed version
ratiomat_center <- ratiomat %>%
  left_join(chr_lengths %>% select(chr, offset, chr_color), by = "chr") %>%
  mutate(cum_pos = window + offset) %>%
  pivot_longer(cols = c("neutral_center", "dav", "new"),
               names_to = "centertype", values_to = "value") %>%
  left_join(chr_shades, by = "chr") %>%
  mutate(center_color_group = paste0(centertype, "_", shade)) %>%
  arrange(chr, cum_pos) %>%
  group_by(chr, centertype) %>%
  mutate(value_smooth = rollmean(value, k = 5, fill = NA, align = "center")) %>%
  ungroup()
center_colors <- c(dav_light = "#1f78b4",
                   dav_dark  = "#1f78b4",
                   new_light   = "#ff7f00",
                   new_dark    = "#ff7f00",
                   neutral_center_light = "#33a02c",
                   neutral_center_dark  = "#33a02c" )


# # Read in test locations for FT genes
# ftgenes <- read_csv("~/Desktop/special_ft_genes.csv") %>%
#   left_join(chr_lengths %>% select(chr, offset), by = "chr") %>%
#   mutate(ps = (bpstart+bpend)/2,
#          cum_pos = ps + offset) %>%
#   filter(genemodel %in% c("mRNA7437", "mRNA10281", "mRNA10313", "mRNA14864","mRNA14865"))

# Read in test locations for FT genes
huegenes <- read_csv("../huegenes-locations.csv") %>%
  left_join(chr_lengths %>% select(chr, offset), by = "chr") %>%
  mutate(ps = (bpstart+bpend)/2,
         cum_pos = ps + offset) %>%
  filter(genemodel == 'mRNA2598')


# Plot slope, highlighting hue region
ggplot(ratiomat_slope, aes(x = cum_pos/1e6, y = value_smooth, color = slope_color_group, group = interaction(chr, slopetype))) +
  annotate("rect", xmin = 34150604/1e6, xmax = 46470553/1e6,
           ymin = -Inf, ymax = Inf,
           fill = "grey80",alpha = 0.3) +
  geom_line() +
  scale_x_continuous(breaks = chr_labels$mid / 1e6,
                     labels = chr_labels$chr,
                     expand = c(0.01, 0, 0.01, 0)) +
  scale_y_continuous(limits = c(0,1), labels = label_number(accuracy = 0.1)) +
  geom_vline(data = huegenes, aes(xintercept = cum_pos / 1e6), color = "red", linetype = "dashed") +
  geom_vline(xintercept = chr_lengths$offset / 1e6, color = "grey50", linetype = "dashed", linewidth = 0.3) +
  labs(x = NULL, y = "Cline slope proportion") +
  scale_color_manual(values = slope_colors, guide = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10.5, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10.5),
        axis.title.x = element_text(size = 11), 
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = grid::unit(c(0.1, 1, 0.1, 1), "mm"))
ggsave(filename = "plot-genomic-clines-steep-vs-shallow-RATIOS-+huegenes-MAIN8.pdf",
       height = 2, width = 8, dpi = 900, units = "in")


# Find top steep windows
top_steep_windows <- ratiomat_slope %>%
  filter(slopetype == "steep",
         n_total > 10) %>%
  arrange(desc(value)) %>%
  slice_head(n = 10)
write.csv(top_steep_windows, file = "result-top-steep-windows-RATIOS.csv")


# Plot center
ggplot(ratiomat_center, aes(x = cum_pos/1e6, y = value_smooth, color = center_color_group, group = interaction(chr, centertype))) +
  geom_line() +
  scale_x_continuous(breaks = chr_labels$mid / 1e6,
                     labels = chr_labels$chr,
                     expand = c(0.01, 0, 0.01, 0)) +
  scale_y_continuous(limits = c(0,1), labels = label_number(accuracy = 0.1)) +
  labs(x = NULL, y = "Cline center proportion") +
  scale_color_manual(values = center_colors, guide = "none") +
  theme_bw() +
  geom_vline(data = huegenes, aes(xintercept = cum_pos / 1e6), color = "red", linetype = "dashed") +
  theme(axis.text.x = element_text(size = 8.0, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8.0),
        panel.grid.minor.x = element_blank(),
        plot.margin = grid::unit(c(0.1, 1, 0.1, 1), "mm"))
ggsave(filename = "plot-genomic-clines-center-skew-RATIOS-AFDfilter+huegenes.pdf",
       height = 2, width = 8, dpi = 900, units = "in")

################################################################################

# (4) Test whether scaffolds are enriched for steep/shallow clines
# With permutation test
################################################################################
# Function to perform randomization test -- is a specific chr enriched for steep clines?
scaffold_enrichment_test <- function(allclines, steepclines, n_perm = 1000) {
  # Step 1: Observed steep cline count per chromosome
  observed_counts <- steepclines %>%
    mutate(chr = as.character(chr)) %>%
    count(chr, name = "observed")
  
  # Step 2: Get chromosome pool for sampling
  chr_pool <- allclines$chr
  
  # Step 3: For each permutation, shuffle steep clines and count per chromosome
  null_counts_all <- replicate(n_perm, {
    shuffled_chr <- sample(chr_pool, size = nrow(steepclines), replace = TRUE)
    tibble(chr = shuffled_chr) %>%
      count(chr, name = "null") %>%
      complete(chr = unique(allclines$chr), fill = list(null = 0)) %>%
      arrange(chr) %>%
      pull(null)
  }, simplify = "matrix")
  
  rownames(null_counts_all) <- sort(as.character(unique(allclines$chr)))
  
  # Step 4: Create summary data.frame
  null_df <- as_tibble(null_counts_all, rownames = "chr") %>%
    pivot_longer(-chr, names_to = "perm", values_to = "null_count") %>%
    group_by(chr) %>%
    summarise(null_mean = mean(null_count),
              p_value = mean(null_count >= observed_counts$observed[match(chr, observed_counts$chr)]),
              .groups = "drop")
  
  # Combine with observed
  result <- observed_counts %>%
    full_join(null_df, by = "chr") %>%
    replace_na(list(observed = 0, null_mean = 0, p_value = 1))
  
  return(result)
}


# Perform 1000 permutations
steep_enrichment <- scaffold_enrichment_test(allclines, steepclines, n_perm = 1000)
write_csv(steep_enrichment, file = "result-STEEP-enrichment-permutation.csv")
shallow_enrichment <- scaffold_enrichment_test(allclines, shallowclines, n_perm = 1000)
write_csv(shallow_enrichment, file = "result-SHALLOW-enrichment-permutation.csv")


################################################################################

# (5) Test whether scaffolds are enriched for offset cline centers
# With permutation test
################################################################################
# Function to perform randomization test -- is a specific chr enriched for offset centers?
scaffold_enrichment_test <- function(allclines, steepclines, n_perm = 1000) {
  # Step 1: Observed steep cline count per chromosome
  observed_counts <- steepclines %>%
    mutate(chr = as.character(chr)) %>%
    count(chr, name = "observed")
  
  # Step 2: Get chromosome pool for sampling
  chr_pool <- allclines$chr
  
  # Step 3: For each permutation, shuffle steep clines and count per chromosome
  null_counts_all <- replicate(n_perm, {
    shuffled_chr <- sample(chr_pool, size = nrow(steepclines), replace = TRUE)
    tibble(chr = shuffled_chr) %>%
      count(chr, name = "null") %>%
      complete(chr = unique(allclines$chr), fill = list(null = 0)) %>%
      arrange(chr) %>%
      pull(null)
  }, simplify = "matrix")
  
  rownames(null_counts_all) <- sort(as.character(unique(allclines$chr)))
  
  # Step 4: Create summary data.frame
  null_df <- as_tibble(null_counts_all, rownames = "chr") %>%
    pivot_longer(-chr, names_to = "perm", values_to = "null_count") %>%
    group_by(chr) %>%
    summarise(null_mean = mean(null_count),
              p_value = mean(null_count >= observed_counts$observed[match(chr, observed_counts$chr)]),
              .groups = "drop")
  
  # Combine with observed
  result <- observed_counts %>%
    full_join(null_df, by = "chr") %>%
    replace_na(list(observed = 0, null_mean = 0, p_value = 1))
  
  return(result)
}


# Perform 1000 permutations
newberryi_enrichment <- scaffold_enrichment_test(plot_allclines, plot_rightclines, n_perm = 1000)
write_csv(newberryi_enrichment, file = "result-NEWBERRYI-introgression-enrichment-permutation.csv")
davidsonii_enrichment <- scaffold_enrichment_test(plot_allclines, plot_leftclines, n_perm = 1000)
write_csv(davidsonii_enrichment, file = "result-DAVIDSONII-introgression-enrichment-permutation.csv")


################################################################################

# 6a. Test whether trait loci are enriched for steep/shallow clines compared to background
# Also look at cline centers
# Special emphasis on comparing clines in the Hue Region vs. the Inflection region vs. background
################################################################################

# Filter to sites, and make two levels of comparison:
# traitclass is specifically the gene models picked out as "trait loci"
# region is the genomic region encompassed by these genomic regions
# And may include loci tagged as "background" in traitclass
plotmat <- allclines %>%
  mutate(traitclass = case_when(phenotype == "modified_hue" ~ "Hue",
                                is.na(phenotype) ~ "Background",
                                TRUE ~ "Other Floral")) %>%
  mutate(traitclass = factor(traitclass, levels = c("Background", "Other Floral", "Hue"))) %>%
  mutate(region = case_when(chr == "2532" & ps >= 37854627 & ps <= 42260154 ~ "Inflection Region",
                            chr == "1086" & ps >= 34150604 & ps <= 46470553 ~ "Hue Region",
                            TRUE ~ "Rest of Genome")) %>%
  arrange(traitclass)

# AFD VS SLOPE, BACKGROUND VS. FLORAL LOCI
ggplot(plotmat, aes(x = AFDbin, y = gradient_50, colour = snpclass)) +
  geom_boxplot(outliers = F) +
  labs(x = "Allele frequency difference bin",
       y = "Cline slope",
       title = "Cline slopes for background vs. trait SNPs\nacross all Allele frequency bins") +
  scale_color_manual(name = "SNP type", labels = c("Background", "Floral loci"),
                     values = c("grey", "orange")) +
  scale_y_continuous(limits = c(0.25, 2.25)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = c(.14, .84),
        legend.title = element_text(face = "bold"),
        legend.background=element_blank(),
        legend.key=element_blank())
ggsave("plot-boxplot-slope-vs-SNPtype-vs-AFDbin.png",
       height = 3.5, width = 4.875, units = "in", dpi = 400)
# (Ratio is ~1.3)

# AFD VS SLOPE, BACKGROUND VS. FLORAL LOCI, SEPARATE HUE
ggplot(plotmat, aes(x = AFDbin, y = gradient_50, colour = traitclass)) +
  geom_boxplot(outliers = F) +
  labs(x = "Allele frequency difference bin",
       y = "Cline slope",
       title = "Cline slopes for background vs. trait SNPs\nacross all Allele frequency bins") +
  scale_color_manual(name = "SNP type", labels = c("Background", "Other floral loci", "Hue"),
                     values = c("grey", "skyblue", "firebrick")) +
  scale_y_continuous(limits = c(0.25, 2.25)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = c(.14, .84),
        legend.title = element_text(face = "bold"),
        legend.background=element_blank(),
        legend.key=element_blank())
ggsave("plot-boxplot-slope-vs-SNPtype-vs-AFDbin-HUESEPARATE.png",
       height = 4.5, width = 5.75, units = "in", dpi = 400)

# AFD VS CENTER, BACKGROUND VS FLORAL LOCI, SEPARATE HUE
ggplot(plotmat, aes(x = AFDbin, y = center_50, colour = traitclass)) +
  geom_boxplot(outliers = F) +
  labs(x = "Allele frequency difference bin",
       y = "Cline center",
       title = "Cline centers for background vs. trait SNPs\nacross all Allele frequency bins") +
  scale_color_manual(name = "SNP type", labels = c("Background", "Other floral loci", "Hue"),
                     values = c("grey", "skyblue", "firebrick")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0,0.2,0), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = c(.14, .84),
        legend.title = element_text(face = "bold"),
        legend.background=element_blank(),
        legend.key=element_blank())
ggsave("plot-boxplot-center-vs-SNPtype-vs-AFDbin-HUESEPARATE.png",
       height = 4.5, width = 5.75, units = "in", dpi = 400)

# AFD VS SLOPE -- REGIONS
ggplot(plotmat, aes(x = AFDbin, y = gradient_50, colour = region)) +
  geom_boxplot(outliers = F) +
  labs(x = "Allele frequency difference bin",
       y = "Cline slope",
       title = NULL) +
  scale_color_manual(name = NULL, values = c("darkred", "aquamarine3", "grey20")) +
  scale_y_continuous(limits = c(0.25, 2.25)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = c(.14, .84),
        legend.title = element_text(face = "bold"),
        legend.background=element_blank(),
        legend.key=element_blank())
ggsave("plot-boxplot-slope-vs-SNPtype-vs-AFDbin-REGIONS.png",
       height = 3, width = 5.75, units = "in", dpi = 600)

# AFD VS CENTER -- REGIONS
ggplot(plotmat, aes(x = AFDbin, y = center_50, colour = region)) +
  geom_boxplot(outliers = F) +
  labs(x = "Allele frequency difference bin",
       y = "Cline center",
       title = NULL) +
  scale_color_manual(name = NULL, values = c("darkred", "aquamarine3", "grey20")) +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0,0.25,0.5,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = c(.14, .84),
        legend.title = element_text(face = "bold"),
        legend.background=element_blank(),
        legend.key=element_blank())
ggsave("plot-boxplot-center-vs-SNPtype-vs-AFDbin-REGIONS.png",
       height = 3, width = 5.75, units = "in", dpi = 600)


# VIOLINS - CENTERS - FLORAL VS HUE
ggplot(plotmat, aes(x = traitclass, y = center_50, colour = traitclass)) +
  geom_violin(draw_quantiles=c(.05, .5, .95)) +
  labs(x = NULL,
       y = NULL,
       title = NULL) +
  scale_color_manual(name = "SNP type", labels = c("Background", "Other Floral", "Hue"),
                     values = c("darkorange", "steelblue3", "darkred")) +
  scale_y_continuous(limits = c(0, 1.25), expand = c(0,0), breaks = c(0, 0.5, 1),
                     labels = label_number(accuracy = 0.1)) +
  geom_signif(comparisons = list(c("Background", "Other Floral"),
                                 c("Hue", "Background"),
                                 c("Other Floral", "Hue")),
              y_position = c(0.97, 1.13, 1.05), textsize = 5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor = element_line(colour = "grey80"))
ggsave("plot-boxplot-center-vs-SNPtype-NOafdbin-HUESEPARATE.png",
       height = 4.5, width = 4.5, units = "in", dpi = 400)

# VIOLINS SLOPES - FLORAL VS HUE
ggplot(plotmat, aes(x = traitclass, y = gradient_50, colour = traitclass)) +
  geom_violin(draw_quantiles=c(.05, .5, .95)) +
  labs(x = NULL,
       y = NULL,
       title = NULL) +
  scale_color_manual(name = "SNP type", labels = c("Background", "Other Floral", "Hue"),
                     values = c("darkorange", "steelblue3", "darkred")) +
  scale_y_continuous(limits = c(0, 2.68), expand = c(0,0), breaks = c(0,1,2),
                     labels = label_number(accuracy = 0.1)) +
  geom_signif(comparisons = list(c("Background", "Other Floral"),
                                 c("Hue", "Background"),
                                 c("Other Floral", "Hue")),
              y_position = c(2.1, 2.40, 2.25), textsize = 5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor = element_line(colour = "grey80"))
ggsave("plot-boxplot-slope-vs-SNPtype-NOafdbin-HUESEPARATE.png",
       height = 4.5, width = 4.5, units = "in", dpi = 400)

# VIOLINS SLOPES - REGIONS -- wilcoxon test (Mann-Whitney U test)
ggplot(plotmat, aes(x = region, y = gradient_50, colour = region)) +
  geom_violin(draw_quantiles=c(.05, .5, .95)) +
  labs(x = NULL,
       y = "Cline gradient",
       title = NULL) +
  scale_color_manual(values = c("darkred", "aquamarine3", "grey20")) +
  scale_y_continuous(limits = c(0, 2.68), expand = c(0,0), breaks = c(0,1,2),
                     labels = label_number(accuracy = 0.1)) +
  geom_signif(comparisons = list(c("Hue Region", "Inflection Region"),
                                 c("Rest of Genome", "Hue Region"),
                                 c("Inflection Region", "Rest of Genome")),
              y_position = c(2.1, 2.40, 2.25), textsize = 5, color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor = element_line(colour = "grey80"))
ggsave("plot-violinplot-slope-vs-SNPtype-NOafdbin-REGIONS.pdf",
       height = 4.5, width = 6, units = "in", dpi = 900)

# VIOLINS CENTERS - REGIONS
ggplot(plotmat, aes(x = region, y = center_50, colour = region)) +
  geom_violin(draw_quantiles=c(.05, .5, .95)) +
  labs(x = NULL,
       y = "Cline center",
       title = NULL) +
  scale_color_manual(values = c("darkred", "aquamarine3", "grey20")) +
  scale_y_continuous(limits = c(0, 1.25), expand = c(0,0), breaks = c(0, 0.5, 1),
                     labels = label_number(accuracy = 0.1)) +
  geom_signif(comparisons = list(c("Hue Region", "Inflection Region"),
                                 c("Rest of Genome", "Hue Region"),
                                 c("Inflection Region", "Rest of Genome")),
              y_position = c(0.97, 1.13, 1.05), textsize = 5, color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor = element_line(colour = "grey80"))
ggsave("plot-violinplot-center-vs-SNPtype-NOafdbin-REGIONS.pdf",
       height = 4.5, width = 6, units = "in", dpi = 900)



################################################################################

# 6b. Then make explicit comparison with Mann-Whitney U test, again stratified by AFD
################################################################################
# MW test -- one per AFDbin
MW_traits <- allclines %>%
  group_by(AFDbin) %>%
  filter(n_distinct(snpclass) == 2) %>%
  summarise(n_trait = sum(snpclass == "trait"),
            n_background = sum(snpclass == "background"),
            p_value_slope = wilcox.test(gradient_50 ~ snpclass)$p.value,
            median_slope_trait = median(gradient_50[snpclass == "trait"]),
            median_slope_background = median(gradient_50[snpclass == "background"]),
            p_value_center = wilcox.test(center_50 ~ snpclass)$p.value,
            median_center_trait = median(center_50[snpclass == "trait"]),
            median_center_background = median(center_50[snpclass == "background"]),
            .groups = "drop")
write.csv(MW_traits, file = "result-mann-whitney-trait-vs-background.csv")
################################################################################

