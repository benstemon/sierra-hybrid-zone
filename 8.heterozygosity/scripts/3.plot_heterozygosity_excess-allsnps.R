library(tidyverse)
library(ggh4x)

setwd("~/8.heterozygosity-selection/data/")

# Reading in HWE files
################################################################################
# Read in hybrids file
# Can be used to amalgamate parents and hybrids
# Must download this from the supplemental data repo
finalplot <- data.frame()
for(i in list.files(pattern = "\\imputed-hybrids.*\\.hwe$")){
  
  # Read in file, make column for ID
  infile <- read_delim(i, delim = "\t") %>%
    select(CHR, POS, `OBS(HOM1/HET/HOM2)`, `E(HOM1/HET/HOM2)`) %>%
    rename(chr = "CHR", ps = "POS") %>%
    mutate(o_hom1 = as.integer(str_split_i(`OBS(HOM1/HET/HOM2)`, "/", 1)),
           o_het  = as.integer(str_split_i(`OBS(HOM1/HET/HOM2)`, "/", 2)),
           o_hom2 = as.integer(str_split_i(`OBS(HOM1/HET/HOM2)`, "/", 3)),
           e_hom1 = as.integer(str_split_i(`E(HOM1/HET/HOM2)`, "/", 1)),
           e_het  = as.integer(str_split_i(`E(HOM1/HET/HOM2)`, "/", 2)),
           e_hom2 = as.integer(str_split_i(`E(HOM1/HET/HOM2)`, "/", 3))) %>%
    mutate(species = str_extract(i, "hybrids|newberryi|davidsonii"),
           chr = str_remove(chr, "scaffold_")) %>%
    select(-c(`OBS(HOM1/HET/HOM2)`, `E(HOM1/HET/HOM2)`))
  
  # Bind to final plot
    finalplot <- rbind(finalplot, infile)
}
remove(infile)
################################################################################


# Testing whether hue locus has excess heterozygosity relative to genomic background
################################################################################

# Specify hue locus
testplot <- finalplot %>%
  mutate(locus = case_when(chr == "1086" & ps >= 34150604 & ps <= 46470553 ~ "hue",
                           TRUE ~ "Background")) %>%
  filter(e_het > 0) %>%
  # Recalculating expectations... vcftools rounds heavily
  mutate(ohetp = o_het/(o_hom1+o_het+o_hom2),
         ehetp = e_het/(o_hom1+o_het+o_hom2),
         excess_het = ohetp - ehetp)

# Issue is... Individuals are not all the same ancestry class... have different
# genome-wide ancestry proportions, which skew expectations for allele frequencies...
# I need to:
# 1. Bin by admixture proportion
# 2. Only include fixed SNPs

# Make summary matrix -- dotplot for this many points is unreasonable
plotmat <- testplot %>%
  group_by(locus) %>%
  summarise(heterozygosity_excess = median(`excess_het`, na.rm = TRUE),
            lower_quartile = quantile(`excess_het`, 0.25, na.rm = TRUE),
            upper_quartile = quantile(`excess_het`, 0.75, na.rm = TRUE),
            min_whisker = quantile(`excess_het`, 0.05, na.rm = TRUE),
            max_whisker = quantile(`excess_het`, 0.95, na.rm = TRUE))

# Boxplot
ggplot(plotmat, aes(x = locus, y = heterozygosity_excess, colour = locus)) +
  geom_boxplot(aes(ymin = min_whisker, lower = lower_quartile, middle = heterozygosity_excess,
                   upper = upper_quartile, ymax = max_whisker),
               stat = "identity", width = 0.5, position = position_dodge(width = 0.6)) +
  labs(x = NULL,
       y = expression("Heterozygosity excess")) +
  # scale_color_manual(name = "Hybrids", labels = c("Rest of 1086 vs. 2532", "Hue vs. Inflection Region"),
  #                    values = c("grey60", "darkred")) +
  # scale_y_continuous(limits = c(0, 0.25)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
ggsave("plot-boxplot-background-vs-hue-excess_het-fixedsnps-hybrids25-75.png", dpi = 600,
       height = 3, width = 3, units = "in")

# Violin
ggplot(testplot, aes(x = locus, y = excess_het, colour = locus)) +
  geom_violin() +
  labs(x = NULL,
       y = expression("Heterozygosity excess")) +
  # scale_color_manual(name = "Hybrids", labels = c("Rest of 1086 vs. 2532", "Hue vs. Inflection Region"),
  #                    values = c("grey60", "darkred")) +
  # scale_y_continuous(limits = c(0, 0.25)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
ggsave("plot-violin-background-vs-hue-excess_het-fixedsnps-hybrids25-75.png", dpi = 600,
       height = 3, width = 3, units = "in")

# Violin heterozygosity
ggplot(testplot, aes(x = locus, y = ohetp, colour = locus)) +
  geom_violin() +
  labs(x = NULL,
       y = expression("Heterozygosity (Observed)")) +
  # scale_color_manual(name = "Hybrids", labels = c("Rest of 1086 vs. 2532", "Hue vs. Inflection Region"),
  #                    values = c("grey60", "darkred")) +
  # scale_y_continuous(limits = c(0, 0.25)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
ggsave("plot-violin-background-vs-hue-observedhet-fixedsnps-hybrids25-75.png", dpi = 600,
       height = 3, width = 3, units = "in")

# Find highest value
t <- testplot %>%
  arrange(desc(excess_het))

################################################################################

# Cross-genome plot of heterozygosity excess
################################################################################
# Make cumulative count for chr lengths and add color
chr_lengths <- read_csv("chrlengths.csv") %>%
  arrange(chr) %>%
  mutate(chr = as.character(chr),
         offset = cumsum(lag(chr_length, default = 0)),
         chr_color = factor(row_number() %% 2))

# Make labels for chromosomes
chr_labels <- chr_lengths %>%
  mutate(mid = offset + chr_length / 2)

chr_colors <- c("1" = "grey", "0" = "black")


# MAke new plot
t <- testplot %>%
  left_join(chr_lengths %>% select(chr, offset, chr_color), by = "chr") %>%
  filter(chr %in% c("1086", "1087", "2532", "2533", "2684", "2685", "2686", "2687")) %>%
  mutate(cum_pos = ps + offset)

# Compute chromosome boundaries (no first, no last)
chr_boundaries <- chr_lengths %>%
  mutate(end_pos = offset + chr_length) %>%
  slice_head(n = nrow(.) - 1)

# Read in location of ANR/DFR
huegenes <- read_csv("~/10.bgchm/data/huegenes-locations.csv") %>%
  mutate(chr = as.character(chr)) %>%
  left_join(chr_lengths %>% select(chr, offset), by = "chr") %>%
  mutate(ps = (bpstart+bpend)/2,
         cum_pos = ps + offset) %>%
  filter(genemodel == 'mRNA2598')


# Cross-genome plot
ggplot(t, aes(x = cum_pos/1e6, y = excess_het, color = chr_color)) +
  annotate("rect", xmin = 34150604/1e6, xmax = 46470553/1e6,
           ymin = -Inf, ymax = Inf,
           fill = "grey80",alpha = 0.3) +
  geom_point(size = 0.3, pch = 16) +
  # geom_hline(yintercept = 0, colour = "red", lty = 2) +
  labs(x = NULL,
       y = "Excess heterozygosity") +
  scale_x_continuous(breaks = chr_labels$mid / 1e6, labels = chr_labels$chr,
                     expand = c(0.01,0,0.01,0)) +
  scale_color_manual(values = chr_colors, guide = "none") +
  # Line at beginning of chr
  geom_vline(data = chr_lengths,
             aes(xintercept = offset / 1e6),
             color = "grey50", linetype = "dashed", linewidth = 0.3) +
  # line at the very end of the chr
  geom_vline(aes(xintercept = (max(chr_lengths$offset + chr_lengths$chr_length) / 1e6)),
             color = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_vline(data = huegenes,
             aes(xintercept = cum_pos / 1e6), color = "red", linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8.0, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8.0),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin=grid::unit(c(0.1,1,0.1,1), "mm"))

ggsave("plot-cross-genome-excess_het-fixedsnps-hybrids25-75.pdf", dpi = 900,
       height = 2, width = 8, units = "in")




# Where is the highest value of excess heterozygosity?
# What is the highest value of heterozygosity, and where is it?
max_het <- t %>%
  arrange(desc(ohetp)) %>%
  slice_head(n = 100)

write.csv(max_het, "result-top-het-sites.csv")


################################################################################\
