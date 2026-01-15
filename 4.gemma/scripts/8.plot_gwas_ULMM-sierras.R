# This script first makes general plots of the ULMM GEMMA results

library(tidyverse)
library(patchwork)
library(ggpubr)

#set directory
# These must be obtained from supplemental data repo
setwd("~/4.gemma/results")

# 1. Code to iteratively plot manhattan plots for all traits:
################################################################################
# Path to directory with results for ULMM models:
filelist = list.files(getwd(), pattern = "*.assoc.txt")

# Path to directory with permutation results:
# These must be obtained from supplemental data repo
permfiles = list.files("quantiles_pvals/", full.names = T)

# Loop to generate generalized manhattan plots for each trait
for (i in 1:length(filelist)){
  
  # Set up list for plotting, if on first iteration
  if(i == 1){
    plot_list <- list()
  }
  
  # Specify file
  infile <- filelist[i]
  
  # Read in GWAS results
  resultfile <- read_delim(infile, delim = "\t") %>%
    na.exclude() %>%
    mutate(chr = str_extract_all(rs, "(?<=:)(.*?)(?=:)")) %>%
    mutate(chr = as.numeric(str_remove_all(chr, "scaffold_"))) %>%
    mutate(ps = as.numeric(str_extract_all(rs, "(?<=:)[^:]+$")))
  
  # Specify the phenotype being examined
  #phenotype <- str_extract(infile, "(?<=normalized\\.).*?(?=\\.assoc)") %>%
  #  str_remove(., "avg_")
  phenotype <- str_extract(infile, "(?<=LDpruned\\.).*?(?=\\.assoc)") %>%
    str_remove(., "avg_")
  
  # Read in the permutation file
  perminfo <- read_delim(permfiles[which(str_detect(string = permfiles, pattern = phenotype))], delim = "\t")
  
  # Format data to be plotted on a single x axis
  #https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
  # Essentially, just adding cumulative bp to each snp position
  data_cum <- resultfile %>%
    group_by(chr) %>%
    summarise(max_bp = max(ps)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>%
    select(chr, bp_add)
  
  resultfile <- resultfile %>%
    inner_join(data_cum, by = "chr") %>%
    mutate(bp_cum = ps + bp_add)
  
  # Set parameters for plotting
  axis_set <- resultfile %>%
    group_by(chr) %>%
    summarize(center = mean(bp_cum))
  
  ylim <- resultfile  %>%
    filter(p_wald == min(p_wald)) %>%
    mutate(ylim = abs(floor(log10(p_wald))) + 2) %>%
    pull(ylim)
  
  # Set p-values
  # First is bonferroni based on number of SNPs
  # Second is based on permutations test, top 0.0001% of top 1% of permutation p-values
  p_bonf = -log10(.05/3021696)
  p_perm_0.0001 = -log10(as.numeric(perminfo[1,2]))
  
  # Plot in ggplot
  manhplot <- ggplot(resultfile, aes(x = bp_cum, y = -log10(p_wald),
                                    color = as_factor(chr),
                                    alpha = -log10(p_wald)/max(-log10(p_wald)))) +
    geom_point(size = 1, shape = 16) +
    geom_hline(yintercept = p_perm_0.0001, colour = "red", linetype = 2) +
    geom_hline(yintercept = p_bonf, colour = "blue", linetype = 2) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(1, ylim+1)) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"),
                                    unique(length(axis_set$chr)))) +
    scale_size_continuous(range = c(0.5, 3)) +
    labs(x = NULL,
         y = expression(-log[10](p)),
         title = phenotype) +
    annotate("text", x = min(resultfile$bp_cum), y = ylim +1, hjust = 0, vjust = 1.5,
             label = paste0("Permutation threshold: ", signif(10^(-p_perm_0.0001), 4)),
             color = "red", size = 3.5) +
    annotate("text", x = min(resultfile$bp_cum), y = ylim +0.2, hjust = 0, vjust = 1.5,
             label = paste0("Bonferroni threshold: ", signif(10^(-p_bonf), 4)),
             color = "blue", size = 3.5) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5))
  
  #add plot to plot list
  plot_list[[i]] <- manhplot
}

# Prepare plot design with arrays of 6 figures.
group_into_six <- function(number) {
  groups <- list()
  current_group <- c()
  for (i in 1:number) {
    current_group <- c(current_group, i)
    if (length(current_group) == 6) {
      groups <- append(groups, list(current_group))
      current_group <- c()
    }
  }
  if (length(current_group) > 0) {
    groups <- append(groups, list(current_group))
  }
  return(groups)
}

# Use this function to index the plot_list
myindex <- group_into_six(length(plot_list))

for(i in 1:length(myindex)){
  #Specify the layout of the plot
  plotdesign <-
    "12
   34
   56"
  
  # Specify the plot itself
  finalplot <-
    plot_list[[myindex[[i]][1]]] + plot_list[[myindex[[i]][2]]] + 
    plot_list[[myindex[[i]][3]]] + plot_list[[myindex[[i]][4]]] + 
    plot_list[[myindex[[i]][5]]] + plot_list[[myindex[[i]][6]]] +
    plot_layout(design = plotdesign)
  
  # Generate the plot
  png(filename = paste("manhattan-plot-ULMM-noparents-normalized-LDpruned",i,"-sierras.png", sep = ""),
      res = 400, units = "in", height = 9, width = 13)
  print(finalplot)
  dev.off()
}
################################################################################



# 3. Focusing on specific traits
################################################################################
# This is the list of important traits for differentiating species
important_traits <- c("modified_hue","avg_tube_height_at_base","avg_tube_height_at_inflection",
                      "avg_longstamen_tube_ratio","chroma","avg_tube_height_at_mouth",
                      "avg_floral_tube_length", "avg_mouth_constriction_angle","avg_pistil_length",
                      "avg_longstamen_pistil_ratio","avg_tube_width_at_base","avg_tube_width_at_mouth", 
                      "avg_mouth_height_width_ratio")

# crease to bottom height, hue, mouth angle, platform angle, stamen length ratio,
# tube height at inflection, and upper petal angle all have "legit" GWAS peaks
# Of that list, hue and tube height at inflection are the only "important" traits
# These must be obtained from supplemental data repo
filelist = c("filtered_ulmm.NOPARENTS_normalized_LDpruned.modified_hue.assoc.txt",
             "filtered_ulmm.NOPARENTS_normalized_LDpruned.avg_tube_height_at_inflection.assoc.txt",
             "filtered_ulmm.NOPARENTS_normalized_LDpruned.avg_stamen_length_ratio.assoc.txt")

# Path to directory with permutation results:
# These must be obtained from supplemental data repo
permfiles = c("modified_hue_significance_quantiles.txt",
              "avg_tube_height_at_inflection_significance_quantiles.txt",
              "avg_stamen_length_ratio_significance_quantiles.txt")

chr_colors <- c("1" = "grey", "0" = "black")


# Loop to generate generalized manhattan plots for each trait
for (i in 1:length(filelist)){
  
  # Specify file
  infile <- filelist[i]
  
  # Read in GWAS results
  resultfile <- read_delim(infile, delim = "\t") %>%
    na.exclude() %>%
    mutate(chr = str_extract_all(rs, "(?<=:)(.*?)(?=:)")) %>%
    mutate(chr = as.numeric(str_remove_all(chr, "scaffold_"))) %>%
    mutate(ps = as.numeric(str_extract_all(rs, "(?<=:)[^:]+$")))
  
  # Specify the phenotype being examined
  phenotype <- str_extract(infile, "(?<=LDpruned\\.).*?(?=\\.assoc)") %>%
    str_remove(., "avg_") %>%
    str_replace_all(., "_", " ") %>%
    str_to_title()
  
  # Read in the permutation file
  perminfo <- read_delim(permfiles[i], delim = "\t")
  
  # Make cumulative count for chr lengths and add color
  chr_lengths <- resultfile %>%
    group_by(chr) %>%
    summarise(chr_length = max(ps)) %>%
    arrange(chr) %>%
    mutate(offset = cumsum(lag(chr_length, default = 0)),
           chr_color = factor(row_number() %% 2))
  
  # Make labels for chromosomes
  chr_labels <- chr_lengths %>%
    mutate(mid = offset + chr_length / 2)
  
  # Refresh with cumulative chr lengths
  plot_resultfile <- resultfile %>%
    left_join(chr_lengths %>% select(chr, offset, chr_color), by = "chr") %>%
    mutate(cum_pos = ps + offset)
  
  # Set y limits
  ylim <- resultfile  %>%
    filter(p_wald == min(p_wald)) %>%
    mutate(ylim = abs(floor(log10(p_wald))) + 2) %>%
    pull(ylim)
  
  # Set p-values -- permutation
  p_perm_0.0001 = -log10(as.numeric(perminfo[1,2]))
  
  # # If phenotype is hue, add line at f3'5'H
  # if(phenotype == "Modified Hue"){
  #   genecoords = (40697697+as.numeric(chr_lengths[2,3]))/1e6
  # }else if(phenotype == "Tube Height At Inflection"){
  #   genecoords = (38440924+as.numeric(chr_lengths[7,3]))/1e6
  # }
  
  # Plot in ggplot
 manhplot <- ggplot(plot_resultfile, aes(x = cum_pos/1e6, y = -log10(p_wald),
                                     color = chr_color,
                                     alpha = -log10(p_wald)/max(-log10(p_wald)))) +
    geom_hline(yintercept = p_perm_0.0001, colour = "red", linetype = 2) +
    geom_point(size = 0.3, shape = 16) +
    # geom_vline(xintercept = genecoords, colour = "orange", linetype = 2) +
    scale_x_continuous(breaks = chr_labels$mid/1e6, labels = chr_labels$chr,
                       expand = c(0.001,0,0.001,0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(1, ylim)) +
    scale_color_manual(values = chr_colors, guide = "none") +
    labs(x = NULL,
         y = expression(-log[10](p)),
         title = NULL) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 0.95),
          axis.text.y = element_text(size = 7),
          axis.title = element_text(size = 7),
          title = element_text(size = 8))
  
  # Save to file
  ggsave(manhplot, filename = paste("PLOT", str_replace_all(phenotype, " ", "-"), "singlemanhattan.tiff", sep = "-"),
         width = 6.5, height = 2, units = "in", dpi = 900)
}
#####
