library(tidyverse)
library(ggh4x)
library(ggtern)
library(ggridges)
library(patchwork)
library(ggpmisc)
library(ggbeeswarm)
setwd("~/6.pigments/")

# 1. Plot relationship between pigment composition and hue
############################################################
# Phenotype data
phenotypes <- read.csv("~/1.phenotype_curation/data/allcombined_phenotype_admixture_subpopulations_spreadsheet-sierras.csv") %>%
  select(-total_brightness)

# Pigment data
pigments <- read_delim("data/pigment_intensity_values.txt") %>%
  filter(str_detect(newID, "129-")) %>%
  select(-species)

# Combine phenotype and pigment data
combplot <- left_join(pigments, phenotypes, by = 'newID')

# Logit transform proportion delphinidin
n <- nrow(combplot)
prop_adj <- (combplot$prop_delph * (n - 1) + 0.5) / n
combplot$logit_prop <- log(prop_adj / (1 - prop_adj))

# Checking residuals
# Fit models
fit_raw <- lm(modified_hue ~ prop_delph, data = combplot)
fit_logit <- lm(modified_hue ~ logit_prop, data = combplot)

# Extract residuals
combplot$resid_raw <- resid(fit_raw)
combplot$resid_logit <- resid(fit_logit)

# Plot residuals side by side
a <- ggplot(combplot) +
  geom_histogram(aes(x = resid_raw), bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Residuals: Raw Proportion", x = "Residual", y = "Count") +
  theme_minimal(base_size = 14)

b <- ggplot(combplot) +
  geom_histogram(aes(x = resid_logit), bins = 30, fill = "darkgreen", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Residuals: Logit-Transformed Proportion", x = "Residual", y = "Count") +
  theme_minimal(base_size = 14)

a | b

# Residuals for raw values look fine -- don't need to transform


# Plot
ggplot(combplot, aes(x = prop_delph, y = modified_hue)) +
  geom_point(aes(color = rgb), size = 3, alpha = 0.9, pch = 16) +
  scale_color_identity() +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", linewidth = 1) +
  stat_poly_eq(use_label(c("eq.label", "adj.R2", "p", "n")),
               label.x = 0.9,
               label.y = 0.9,
               size = 4.5,
               parse = TRUE) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  labs(x = "Proportion Delphinidin",
       y = "Modified Hue",
       title = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"))

# ggsave("../8.pigments/PLOTS/PLOT-pigment-vs-hue.pdf", dpi = 900,
#        height = 8, width = 8)
############################################################

# 2. Set up combined GT and pigment data frames
# Get f3'5'h GT and DFRlike/DFR GT data
################################################################################
# f3'5'h genotype datas
pcrgts <- read_csv("data/pcr-f3p5ph-genotypes.csv") %>%
  select(sample, genotype) %>%
  setNames(c("newID", "PCRgt"))

# DFRlike/DFR data:
# (a) Find AFDs in parent species of significant GEMMA hits at DFRlike/DFR
tmpgts <- read_delim("data/DFRlike_GEMMAhits.txt", col_names = F) %>%
  mutate(X1 = str_remove(X1, "scaffold_"))

# These must be obtained from supplemental data repo
AFDparents <- read_delim("AFD-table-parents.txt") %>%
  filter(chr == "1086",
         ps %in% tmpgts$X2)

# Write parent AFDs for DFRlike-DFR GEMMA SNPs to file
# All of the SNPs are either fixed or very nearly fixed
# newberryi is second column. When not fixed in newberryi, 1/52 alleles
# when not fixed in davidsonii, 1/94 alleles
write_delim(AFDparents, "DFRlike-SNPs-AFDs-parents.tsv", delim = '\t')


# Create new DFRlike-DFR matrix with genotype calls for whole locus given SNPs
DFRlikegts <- read_delim("data/final_DFRlike_genotype_matrix.tsv") %>%
  pivot_longer(cols = -site,
               names_to = "sample",
               values_to = "gt") %>%
  mutate(gt = case_when(gt %in% c("0/0", "0|0") ~ 0,
                        gt %in% c("0/1", "1/0", "0|1", "1|0") ~ 1,
                        gt %in% c("1/1", "1|1") ~ 2,
                        TRUE ~ NA)) %>%
  group_by(sample) %>%
  summarise(total_gt = sum(!is.na(gt)),
            new_DFRlike = mean(gt == 2, na.rm = TRUE),
            het_DFRlike = mean(gt == 1, na.rm = TRUE),
            dav_DFRlike = mean(gt == 0, na.rm = TRUE)) %>%
  mutate(sample = str_remove_all(sample, "BWS_")) %>%
  mutate(gt_DFRlike = case_when(het_DFRlike >= 0.5 ~ "H",
                            new_DFRlike > 0.5 ~ "N",
                            dav_DFRlike > 0.5 ~ "D")) %>%
  rename(newID = "sample") %>%
  mutate(pop = str_extract(newID,  ".*(?=-)"))

# Write combined matrix to disk
write_csv(DFRlikegts, file = "result-DFRlikegts.csv")


# Combine the GT and phenotype/pigment data
combplot <- left_join(phenotypes, pigments, by = 'newID') %>%
  left_join(., pcrgts, by = 'newID') %>%
  left_join(., DFRlikegts, by = 'newID')
write_csv(combplot, file = "result-COMBINED-GT-PIGMENT-TABLE.csv")

################################################################################


# 3. Statistical modeling of effect of GT on pigment and hue
# HUE
################################################################################
# f3'5'h alone
m1 <- lm(data = combplot, modified_hue ~ PCRgt)

# additive model (independent [main] effects)
m2 <- lm(data = combplot, modified_hue ~ PCRgt + gt_DFRlike)

# interaction model -- main effects AND interaction effects [how effect of e.g. f3'5'h depends on genotype at DFRlike, and vice versa]
m3 <- lm(data = combplot, modified_hue ~ PCRgt * gt_DFRlike)

# ANOVA
anova(m1, m2, m3)
capture.output(anova(m1, m2, m3), file = "result-anova-DFRlike-DFR.txt")

# AIC
AIC(m1, m2, m3)
capture.output(AIC(m1, m2, m3), file = "result-AIC-DFRlike-DFR.txt")

# model summaries -- write best model
summary(m1)
summary(m2)
summary(m3)
capture.output(summary(m2), file = "result-summary-DFRlike-DFR.txt")

# Write mean values of hue and delphinidin given genotype
hue_summary <- combplot %>%
  group_by(PCRgt, gt_DFRlike) %>%
  summarize(mean_hue = mean(modified_hue, na.rm = TRUE),
            mean_delph = mean(prop_delph, na.rm = TRUE),
            sd_hue = sd(modified_hue, na.rm = TRUE),
            n = n())
write_csv(hue_summary, file = "result-mean-hue-pigment-per-GT.csv")

################################################################################


