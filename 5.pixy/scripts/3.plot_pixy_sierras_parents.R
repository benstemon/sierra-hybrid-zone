library(tidyverse)
library(zoo)
library(patchwork)

# Change to directory of interest
setwd("5.pixy/")

# GENOME-WIDE:
################################################################################
# Read in files, manipulate
# These must be obtained from supplemental data repo
pi <- read_delim("sierraparents_nomask_100kb_pi.txt", delim = "\t") %>%
  select(-c(no_sites, count_diffs, count_comparisons, count_missing)) %>%
  mutate(pop = case_when(pop == "newberryi" ~ "newberryi_pi",
                         pop == "davidsonii" ~ "davidsonii_pi",
                         TRUE ~ pop)) %>%
  pivot_wider(., names_from = pop, values_from = avg_pi)
dxy <- read_delim("sierraparents_nomask_100kb_dxy.txt", delim = "\t")
fst <- read_delim("sierraparents_nomask_100kb_fst.txt", delim = "\t")

# First, find genome-wide averages for dxy
genomewide_dxy <- sum(dxy$count_diffs, na.rm = TRUE) / sum(dxy$count_comparisons, na.rm = TRUE)
cat(genomewide_dxy, file = "genome-wide-dxy-100kb.txt")

# Join divergence metrics into single df
divergence_metrics <- full_join(dxy, fst,
                                by = c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2")) %>%
  full_join(., pi, by = c("chromosome", "window_pos_1", "window_pos_2")) %>%
  mutate(ps = (window_pos_1 + window_pos_2)/2) %>%
  rename(chr = chromosome) %>%
  arrange(chr)

# To set up plotting similar to GWAS, add cumulative bp for each chr
data_cum <- divergence_metrics %>%
  group_by(chr) %>%
  summarise(max_bp = max(ps)) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
  select(chr, bp_add)
################################################################################

# PLOT GENOME-WIDE ... NEW WAY
################################################################################
# Pick a window size for smoothing (adjust as needed)
roll_window <- 40

# Calculate rolling averages for each metric
rolling_metrics <- divergence_metrics %>%
  arrange(chr, ps) %>%
  group_by(chr) %>%
  mutate(avg_dxy_filled   = na.approx(avg_dxy, x = ps, na.rm = FALSE),
         fst_filled       = na.approx(avg_wc_fst, x = ps, na.rm = FALSE),
         newb_pi_filled   = na.approx(newberryi_pi, x = ps, na.rm = FALSE),
         david_pi_filled  = na.approx(davidsonii_pi, x = ps, na.rm = FALSE),
         avg_dxy_roll   = rollmean(avg_dxy_filled, k = roll_window, fill = NA, align = "center"),
         fst_roll       = rollmean(fst_filled, k = roll_window, fill = NA, align = "center"),
         newb_pi_roll   = rollmean(newb_pi_filled, k = roll_window, fill = NA, align = "center"),
         david_pi_roll  = rollmean(david_pi_filled, k = roll_window, fill = NA, align = "center")) %>%
  ungroup() %>%
  inner_join(data_cum, by = "chr") %>%
  mutate(bp_cum = ps + bp_add)


# Scale FST to plot on a second y-axis
# Pick a factor to bring FST into similar range as dxy/pi
fst_scale <- max(rolling_metrics$avg_dxy_roll, 
                 rolling_metrics$newb_pi_roll,
                 rolling_metrics$david_pi_roll, na.rm = TRUE) / max(rolling_metrics$fst_roll, na.rm = TRUE) +
  0.01

rolling_metrics <- rolling_metrics %>%
  mutate(fst_scaled = fst_roll * fst_scale)  %>%
  mutate(chr = str_remove(chr, "scaffold_"))

# Set parameters for plotting
axis_set <- rolling_metrics %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))


# Read in gff with gene information
# These must be obtained from supplemental data repo
gff3 <- read_delim("single_isoform_davidsonii_FUNCTIONAL-INCLUDED.gff", col_names = F) %>%
  select(X1, X3, X4, X5) %>%
  setNames(c("chr", "genemodel", "bpstart", "bpend")) %>%
  filter(genemodel == "gene") %>%
  filter(chr %in% unique(dxy$chromosome)) %>%
  mutate(midpoint = ceiling((bpstart+bpend)/2))

# Join to cumulative count
gene_windows <- gff3 %>%
  inner_join(data_cum, by = "chr") %>%
  mutate(bp_cum = midpoint + bp_add) %>%
  mutate(chr = str_remove(chr, "scaffold_"))

# PLOT 1: FST, DXY, and PI
p1 <- ggplot(rolling_metrics, aes(x = bp_cum)) +
  geom_line(aes(y = avg_dxy_roll, color = "DXY"), size = 0.6) +
  geom_line(aes(y = newb_pi_roll, color = "newberryi π"), size = 0.6) +
  geom_line(aes(y = david_pi_roll, color = "davidsonii π"), size = 0.6) +
  geom_line(aes(y = fst_scaled, color = "FST"), size = 0.6) +
  # scaffold delimiters
  geom_vline(xintercept = data_cum$bp_add, color = "grey50", linetype = "dashed", linewidth = 0.3) +
  scale_y_continuous(name = expression(D[XY] ~ "/ π"),
                     sec.axis = sec_axis(~ . / fst_scale, name = expression(F[ST]))) +
  scale_color_manual(values = c("DXY" = "darkgrey",
                                "newberryi π" = "#c9368c",
                                "davidsonii π" = "#6245ba",
                                "FST" = "black"),
                     labels = c("DXY" = expression(D[XY]),
                                "newberryi π" = "newberryi π",
                                "davidsonii π" = "davidsonii π",
                                "FST" = expression(F[ST])),
                     breaks = c("FST", "DXY", "davidsonii π", "newberryi π")) +
  scale_x_continuous(label = axis_set$chr,
                     breaks = axis_set$center,
                     expand = c(0, 0)) +
  labs(x = NULL, color = NULL) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 5, r = 10, b = 0, l = 10))

# PLOT 2: GENE DENSITY HISTOGRAMS
p2 <- ggplot(gene_windows, aes(x = bp_cum)) +
  geom_histogram(binwidth = 1000000, fill = "darkgrey", color = "black") +  # binwidth in bp, adjust as needed
  geom_vline(xintercept = data_cum$bp_add, color = "grey50", linetype = "dashed", linewidth = 0.3) +
  scale_x_continuous(label = axis_set$chr,
                     breaks = axis_set$center,
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Gene count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(t = 0, r = 10, b = 5, l = 10))

combined_plot <- p1 / p2 + plot_layout(heights = c(1.5, 1))
combined_plot 
ggsave("PLOT-100kb-combined-fst-dxy-geneplot.pdf", device = cairo_pdf, height = 4.75, width = 11)

################################################################################



# HUE -- SETUP
################################################################################
# Read in GWAS file
# These must be obtained from supplemental data repo
gwaspeaks <- read_delim("filtered_ulmm.NOPARENTS_normalized_LDpruned.modified_hue.assoc.txt") %>%
  na.exclude() %>%
  mutate(chr = str_extract_all(rs, "(?<=:)(.*?)(?=:)")) %>%
  mutate(chr = as.numeric(str_remove_all(chr, "scaffold_"))) %>%
  mutate(ps = as.numeric(str_extract_all(rs, "(?<=:)[^:]+$"))) %>%
  filter(chr %in% "1086")

# Also, a priori list gene models of interest...
# This includes flavonoid/anthocyanin genes, MYBs, bHLHs, and DIV
# MYBs and BHLHs included:

# IPN2 mRNA2487
# PHL7 mRNA2495
# DIVARICATA mRNA2706
# BHLH62 mRNA2800
# GAM1 mRNA2938
# BHLH25 mRNA2943
# MYB46 mRNA3016, mRNA3018
genes_of_interest <- c("mRNA2487", "mRNA2495", "mRNA2598", "mRNA2653", "mRNA2706", 
                       "mRNA2745", "mRNA2746", "mRNA2800", "mRNA2815", "mRNA2816", 
                       "mRNA2817", "mRNA2818", "mRNA2838", "mRNA2938", "mRNA2941", 
                       "mRNA2943", "mRNA3016", "mRNA3018")

# Flavonoid and Anthocyanin biosynthesis ONLY
# ANR mRNA2598 (annotated as DFR)
# CYP93B2 mRNA2653
# C4H2 mRNA2745, mRNA2746
# UGT87A1 / UGT87A2 mRNA2815, mRNA2816, mRNA2817, mRNA2818
# F3'5'H mRNA2838
# F6H2-1-1 mRNA2941
genes_of_interest <- c("mRNA2598", "mRNA2653", "mRNA2745", "mRNA2746", "mRNA2815",
                       "mRNA2816", "mRNA2817", "mRNA2818", "mRNA2838", "mRNA2941")

#bHLH TF of interest
genes_of_interest <- c("mRNA2593")


# READ IN DIVERGENCE STATISTICS
# For 100bp windows
# These must be obtained from supplemental data repo
setwd("pixyout-nomask_f3p5ph_zoom1/")

pifile <- list.files(pattern = "*pi.txt")
dxyfile <- list.files(pattern = "*dxy.txt")
fstfile <- list.files(pattern = "*fst.txt")

pi <- read_delim(pifile, delim = "\t") %>%
  select(-c(no_sites, count_diffs, count_comparisons, count_missing)) %>%
  mutate(pop = case_when(pop == "newberryi" ~ "newberryi_pi",
                         pop == "davidsonii" ~ "davidsonii_pi",
                         TRUE ~ pop)) %>%
  pivot_wider(., names_from = pop, values_from = avg_pi)
dxy <- read_delim(dxyfile, delim = "\t")
fst <- read_delim(fstfile, delim = "\t")


# READ IN GENOME ANNOTATIONS
annotations <- read_delim("~/project_storage/project_comparative_genome/genomes/davidsonii/single_isoform/single_isoform_davidsonii_FUNCTIONAL-INCLUDED.gff",
                          delim = "\t", col_names = F) %>%
  select(-c(2,6,8)) %>%
  rename(chr=X1, model=X3, bpstart=X4, bpend=X5, strand = X7, info=X9) %>%
  mutate(chr = str_remove(chr, "scaffold_")) %>%
  filter(chr %in% "1086") %>%
  filter(model == "mRNA")


# Set boundaries for gene of interest, that's our focal point. Also set windowsize
# Original window was 36500000 - 44000000 (for pixy) -- ~50kb window
# We can just view the whole window
# Also set significance threshold
windows=c(36500000,44000000)
sigthresh=-log10(6.240259548938e-09)

# New annotations object for easy filtering
plotting_annotations <- annotations %>%
  filter(bpstart > windows[1] & bpend < windows[2]) %>%
  mutate(bpstart = bpstart/1e6, bpend = bpend/1e6) %>%
  mutate(interest = case_when(str_detect(info, paste(genes_of_interest, collapse = "|")) ~ "yes",
                              TRUE ~ "no"))

# New GWAS object for easy filtering and reducing re-reading it in
plotting_gwas <- gwaspeaks %>%
  select(c(ps, p_wald)) %>%
  filter(ps > windows[1] & ps < windows[2]) %>%
  pivot_longer(., cols = p_wald, names_to = "metric", values_to = "p_wald") %>%
  mutate(ps = ps/1e6) 

# Join pi, dxy into one data.frame
# Also calculate a rolling mean for each value
divergence_metrics <- full_join(dxy, fst,
                                by = c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2")) %>%
  full_join(., pi, by = c("chromosome", "window_pos_1", "window_pos_2")) %>%
  mutate(ps = (window_pos_1 + window_pos_2)/2) %>%
  rename(chr = chromosome) %>%
  arrange(chr) %>%
  pivot_longer(., cols = c(avg_dxy, newberryi_pi, davidsonii_pi),
               names_to = "metric", values_to = "value") %>%
  filter(ps > windows[1] & ps < windows[2]) %>%
  group_by(metric) %>%
  mutate(rolling_mean = rollmean(value, k = 5, fill = "extend", align = "center")) %>%
  ungroup() %>%
  mutate(ps = ps/1e6)
################################################################################

# Hue region 1: ANR
#####
# Now we are focusing on specific gene clusters -- x coordinates
divergence_anr <- divergence_metrics %>%
  mutate(window_pos_1 = window_pos_1/1e6,
         window_pos_2 = window_pos_2/1e6) %>%
  filter(window_pos_1 >=37.305 & window_pos_2 <=37.325)

# Let's just plot the region and see what it looks like
# Make a scaling factor because we want to plot these on the same plot
scalefactor <- (max(na.omit(divergence_anr$rolling_mean)) / max(-log10(plotting_gwas$p_wald)))

ggplot(divergence_anr, aes(x = ps, y = rolling_mean, col = metric)) +
  geom_point(data = plotting_gwas, size = 0.75,
             aes(x = ps, y = (-log10(p_wald)*scalefactor),
                 color = (-log10(p_wald) > sigthresh))) +
  labs(title = NULL,
       x = NULL) +
  scale_x_continuous(limits = c(min(divergence_anr$window_pos_1), max(divergence_anr$window_pos_2)),
                     expand = c(0.01, 0)) +
  scale_y_continuous(name = "Rolling Mean (Divergence Metrics)",
                     sec.axis = sec_axis(~./scalefactor, name = expression(-log[10](p)), ),
                     limits = c(0,0.065)) +
  geom_line() +
  geom_segment(data = plotting_annotations %>% filter(strand == "+"),
               aes(x = bpstart, xend = bpend, y = 0.13, yend = 0.13, col = interest),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               inherit.aes = FALSE) +
  geom_segment(data = plotting_annotations %>% filter(strand == "-"),
               aes(x = bpend, xend = bpstart, y = 0.06, yend = 0.06, col = interest),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               inherit.aes = FALSE,
               linewidth = 1) +
  scale_color_manual(values = c("davidsonii_pi" = "#6245ba", "newberryi_pi" = "#c9368c",
                                "avg_dxy" = "orange", "TRUE" = "darkred",
                                "FALSE" = "lightgrey", "yes" = "red", "no" = "lightgrey"),
                     breaks = c("davidsonii_pi", "newberryi_pi", "avg_dxy"),
                     labels = c(expression(italic("davidsonii") * " " * pi),
                                expression(italic("newberryi") * " " * pi),
                                expression(italic(d)[xy])),
                     name = NULL) +
  annotate("text", x = 37.3123, y = 0.063, label = "ANR", size = 5, fontface = "bold.italic", colour = "red") +
  geom_hline(yintercept = sigthresh*scalefactor, colour = "red", linetype = 2) +
  #geom_blank(aes(y=max(na.omit(rolling_mean)) * 1.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y.right = element_text(angle = 90))
ggsave(file = "PLOT-PIXY-hue_ANR.pdf",
       width = 6, height = 4, dpi = 1000)
#####

# Hue region 2: F3'5'H
#####
# Now we are focusing on specific gene clusters -- x coordinates
divergence_f3p5ph <- divergence_metrics %>%
  mutate(window_pos_1 = window_pos_1/1e6,
         window_pos_2 = window_pos_2/1e6) %>%
  filter(window_pos_1 >=40.690 & window_pos_2 <=40.705)

# Let's just plot the region and see what it looks like
# Make a scaling factor because we want to plot these on the same plot
scalefactor <- (max(na.omit(divergence_anr$rolling_mean)) / max(-log10(plotting_gwas$p_wald)))

ggplot(divergence_f3p5ph, aes(x = ps, y = rolling_mean, col = metric)) +
  geom_point(data = plotting_gwas, size = 0.75,
             aes(x = ps, y = (-log10(p_wald)*scalefactor),
                 color = (-log10(p_wald) > sigthresh))) +
  labs(title = NULL,
       x = NULL) +
  scale_x_continuous(limits = c(min(divergence_f3p5ph$window_pos_1), max(divergence_f3p5ph$window_pos_2)),
                     expand = c(0.01, 0)) +
  scale_y_continuous(name = "Rolling Mean (Divergence Metrics)",
                     sec.axis = sec_axis(~./scalefactor, name = expression(-log[10](p)), ),
                     limits = c(0, 0.065)) +
  geom_line() +
  geom_segment(data = plotting_annotations %>% filter(strand == "+"),
               aes(x = bpstart, xend = bpend, y = 0.06, yend = 0.06, col = interest),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               inherit.aes = FALSE,
               linewidth = 1) +
  geom_segment(data = plotting_annotations %>% filter(strand == "-"),
               aes(x = bpend, xend = bpstart, y = 0.07, yend = 0.07, col = interest),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               inherit.aes = FALSE) +
  scale_color_manual(values = c("davidsonii_pi" = "#6245ba", "newberryi_pi" = "#c9368c",
                                "avg_dxy" = "orange", "TRUE" = "darkred",
                                "FALSE" = "lightgrey", "yes" = "red", "no" = "lightgrey"),
                     breaks = c("davidsonii_pi", "newberryi_pi", "avg_dxy"),
                     labels = c(expression(italic("davidsonii") * " " * pi),
                                expression(italic("newberryi") * " " * pi),
                                expression(italic(d)[xy])),
                     name = NULL) +
  annotate("text", x = 40.6985, y = 0.063, label = "F3'5'H", size = 5, fontface = "bold.italic", colour = "red") +
  geom_hline(yintercept = sigthresh*scalefactor, colour = "red", linetype = 2) +
  #geom_blank(aes(y=max(na.omit(rolling_mean)) * 1.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y.right = element_text(angle = 90))
ggsave(file = "PLOT-PIXY-hue_F3p5ph.pdf",
       width = 6, height = 4, dpi = 1000)
#####
