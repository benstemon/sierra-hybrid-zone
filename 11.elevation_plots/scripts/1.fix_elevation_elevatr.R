library(tidyverse)
library(elevatr)
library(sf)

setwd("~/11.elevation_plots/")

# Read in current lat/lon coordinates
# And rearrange so lat is first, lon is second, and rename
setwd("~/Desktop/fix_elevation/")
latlon <- read_delim("data/sierras-pops-elevation.txt") %>%
  rename(x = lon, y = lat, oldelevation = elevation) %>%
  relocate(c(x, y)) %>%
  as.data.frame(.)

# Estimate new elevation
t <- get_elev_point(latlon, prj = "EPSG:4326", units = "feet")

# Write to new file
newlatlon <- t %>%
  rename(newelevation = elevation) %>%
  select(c(accession, newelevation, oldelevation)) %>%
  st_drop_geometry(.)

write_delim(newlatlon, file = "data/sierras-pops-NEWelevation.txt", delim = "\t")
