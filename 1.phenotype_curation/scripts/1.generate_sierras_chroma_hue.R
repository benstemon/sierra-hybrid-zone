#notes:
#1. BRIGHTNESS is AUC from 400-700 nm

#first separate spectrum among four equal segments
#blue -- 400-475
#green -- 475-550
#yellow -- 550-625
#red -- 625-700

#relative brightness of each segment determined by
#brightness in segment divided by total brightness


#convert to LM contrast and MS contrast:
#red - green
#yellow - blue


#plotting MS vs LM, CHROMA (saturation) = distance from origin
#sqrt(LM^2 + Ms^2)

#HUE (type of color) = angle clockwise from the y-axis
#red colors are low hue, while yellow, green, blue, and purple 
#have progressively higher values


#Begin computation
setwd("1.phenotype_curation/")
library(pavo)
library(tidyverse)

#load spec data
# This is available in supplemental data repo
rawdata <- getspec("~/sierras_spec_data/",
                   ext = "txt")


#Generate chroma brightness hue for each FLOWER
################################################################################
#clean up the raw data
#remove decoloration individuals (if present)
#transform data so that lowest value is 0 and highest is 1
#rename rownames to be set to wl, and remove wl column
cleandata <- procspec(rawdata, fixneg = "zero", opt = c("min", "max")) %>%
  column_to_rownames("wl")
  

#correct the column names, with 3-digit individual IDS
newcolnames <- colnames(cleandata) %>%
  str_split(., "-") %>%
  map(., ~{
    if (endsWith(.[[2]], "A") || endsWith(.[[2]], "B")) {
      width <- 4
    } else {
      width <- 3
    }
    .[[2]] <- str_pad(.[[2]], width = width, side = "left", pad = "0")
    return(.)
  }) %>%
  map(., ~paste(.x, collapse = "-")) %>%
  unlist()

#just a little sanity check to ensure this worked correctly
for(i in 1:length(colnames(cleandata))){
  print(paste(colnames(cleandata)[i], newcolnames[i]))
}

#replace column names
colnames(cleandata) <- newcolnames

write.csv(cleandata, "data/sierras_rawspec_cleaned_allflowers.csv")





#for loop to generate chroma, brightness, hue for each sample
for (i in 1:ncol(cleandata)){
  #blue, green, yellow, red brightness -- total values
  b_raw <- sum(cleandata[as.character(400:475),i])
  g_raw <- sum(cleandata[as.character(475:550),i])
  y_raw <- sum(cleandata[as.character(550:625),i])
  r_raw <- sum(cleandata[as.character(625:700),i])
  
  #total brightness
  brightness <- sum(b_raw, g_raw, y_raw, r_raw)
  
  #relative brightness
  #doing this way because otherwise doesn't sum to 1
  B <- b_raw/brightness
  G <- g_raw/brightness
  Y <- y_raw/brightness
  R <- r_raw/brightness
  
  #calculate chroma
  LM = R - G
  MS = Y - B
  chroma <- sqrt(LM^2 + MS^2)
  
  #calculate hue
  hue_radians <- (sign(MS)*acos(LM/chroma)) %% (2*pi)
  hue_degrees <- ((sign(MS)*acos(LM/chroma)) %% (2*pi)) * (180/pi)
  
  #If first iteration, make the outlist
  if(i == 1){
    outlist <- list()
  }
  
  #Write to the outlist
  outlist[[colnames(cleandata[i])]] <- data.frame("blue_brightness" = B,
                                                  "green_brightness" = G,
                                                  "yellow_brightness" = Y,
                                                  "red_brightness" = R,
                                                  "total_brightness" = brightness,
                                                  "chroma" = chroma,
                                                  "hue_radians" = hue_radians,
                                                  "hue_degrees" = hue_degrees)
}

#Turn the list into a data frame which we can then export
out_df <- do.call(rbind, outlist)

#add newID column, which is just the rowname
out_df$newID <- rownames(out_df)

#add pop column and reorganize
out_df <- out_df %>%
  rowwise() %>%
  mutate(pop = str_split(newID, "-")[[1]][1]) %>%
  select(newID, pop, everything()) %>%
  mutate(modified_hue = ifelse(hue_degrees > 270, (hue_degrees + 90) - 360, hue_degrees + 90))

write.csv(out_df, file = "data/ALLSAMPLES_sierras_color_chroma_hue_allflowers.csv", row.names = F)
################################################################################





#Generate chroma brightness hue for AVERAGED flowers
################################################################################
#find names of flowers that need averaged
no_decolor <- rawdata %>%
  select(-contains("old"))

a_names <- names(no_decolor)[grep("A", names(no_decolor))]
b_names <- names(no_decolor)[grep("B", names(no_decolor))]

#make new df with only flowers to be averaged (also sort and add wl column to front)
avg_flowers <- no_decolor %>%
  select(contains(sort(c(a_names, b_names)))) %>%
  mutate(wl = no_decolor$wl, .before = names(.)[1])
  
  
#verify that the lengths of the dfs are correct
length(a_names) + length(b_names) +1 == ncol(avg_flowers)

# IF THIS IS NOT CORRECT, THE NAMING SCHEME IS WRONG:
#####

# find the names that "should" exist, but don't (no "B" flower for an "A" flower, vice versa)
shouldb_names <- str_replace(a_names, "A", "B")
shoulda_names <- str_replace(b_names, "B", "A")

# Manually change file names that fall into these categories until everything is fixed
a_names[!(a_names %in% shoulda_names)]
b_names[!(b_names %in% shouldb_names)]
#####

#also find flowers that do not need averaged (and do not keep the wl column)
noavg_flowers <- no_decolor %>%
  select(-contains(names(avg_flowers))) %>%
  select(-contains("wl"))



#pavo command to find averaged values:
avg_flowers <- aggspec(avg_flowers, by = 2, FUN = mean)

#combine all values into new data set, clean up names, and normalize
# Changing to minmax and fixneg has no effect on hue calculations
# But does change brightness and chroma calculations
cleandata_average <- cbind(avg_flowers, noavg_flowers) %>%
  rename_all(~gsub("A", "", .)) %>%
  procspec(., opt = c("max")) %>%
  column_to_rownames("wl")


#correct the column names, with 3-digit individual IDS
newcolnames <- colnames(cleandata_average) %>%
  str_split(., "-") %>%
  map(., ~{
    if (endsWith(.[[2]], "A") || endsWith(.[[2]], "B")) {
      width <- 4
    } else {
      width <- 3
    }
    .[[2]] <- str_pad(.[[2]], width = width, side = "left", pad = "0")
    return(.)
  }) %>%
  map(., ~paste(.x, collapse = "-")) %>%
  unlist()

#just a little sanity check to ensure this worked correctly
for(i in 1:length(colnames(cleandata_average))){
  print(paste(colnames(cleandata_average)[i], newcolnames[i]))
}


# #replace the old column names and write to file
# colnames(cleandata_average) <- newcolnames
# write.csv(cleandata_average, "sierras_rawspec_cleaned_averageflowers2.csv")
# NOTE -- missing ALL PHENOTYPE data for 133-013. Unfortunate.


#GENERATE CHROMA, BRIGHTNESS, HUE
#for loop to generate chroma, brightness, hue for each sample
for (i in 1:ncol(cleandata_average)){
  #blue, green, yellow, red brightness -- total values
  b_raw <- sum(cleandata_average[as.character(400:475),i])
  g_raw <- sum(cleandata_average[as.character(475:550),i])
  y_raw <- sum(cleandata_average[as.character(550:625),i])
  r_raw <- sum(cleandata_average[as.character(625:700),i])
  
  #total brightness
  brightness <- sum(b_raw, g_raw, y_raw, r_raw)
  
  #relative brightness
  #doing this way because otherwise doesn't sum to 1
  B <- b_raw/brightness
  G <- g_raw/brightness
  Y <- y_raw/brightness
  R <- r_raw/brightness
  
  #calculate chroma
  LM = R - G
  MS = Y - B
  chroma <- sqrt(LM^2 + MS^2)
  
  #calculate hue
  hue_radians <- (sign(MS)*acos(LM/chroma)) %% (2*pi)
  hue_degrees <- ((sign(MS)*acos(LM/chroma)) %% (2*pi)) * (180/pi)
  
  #If first iteration, make the outlist
  if(i == 1){
    outlist <- list()
  }
  
  #Write to the outlist
  outlist[[colnames(cleandata_average[i])]] <- data.frame("blue_brightness" = B,
                                                          "green_brightness" = G,
                                                          "yellow_brightness" = Y,
                                                          "red_brightness" = R,
                                                          "total_brightness" = brightness,
                                                          "chroma" = chroma,
                                                          "hue_radians" = hue_radians,
                                                          "hue_degrees" = hue_degrees)
}


#Turn the list into a data frame which we can then export
# Add NEWID column and pop column
out_df <- do.call(rbind, outlist) %>%
  mutate(newID = rownames(.)) %>%
  rowwise() %>%
  mutate(pop = str_split(newID, "-")[[1]][1]) %>%
  select(newID, pop, everything()) %>%
  arrange(newID) %>%
  mutate(modified_hue = ifelse(hue_degrees > 270, (hue_degrees + 90) - 360, hue_degrees + 90))

write.csv(out_df, file = "data/ALLSAMPLES_sierras_color_chroma_hue_averageflowers.csv", row.names = F)


################################################################################

#make rgb
rgbinfile <- cleandata_average %>%
  mutate(wl = avg_flowers$wl) %>%
  relocate(wl)
rgbvalues <- as.data.frame(spec2rgb(as.rspec(rgbinfile, lim = c(300, 700)))) %>%
  rename(rgb = 1) %>%
  mutate(newID = rownames(.))

# Because we can easily distinguish everything with a single axis, don't need 2D
# Rewrite the output matrix
finaldf <-full_join(x = out_df, y = rgbvalues)
write.csv(finaldf, file = "data/ALLSAMPLES_sierras_color_chroma_hue_averageflowers.csv", row.names = F)



