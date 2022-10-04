### checking values in figures


# packages

library(tidyverse)
library(tidylog)
library("easystats")
library(scales)

## directory for figures
out.dir <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V2/Figures/"

## data

load(file = "output_data/sync/04a_all_sync_for_figures.RData") #allsync

# plotting over distance----------------------------------------------------------------

head(allsync)
unique(allsync$Region)

syncVals <- allsync %>%
  group_by(Region, Trait, Connectivity) %>%
  summarise(minSync = min(na.omit(Sync)),
            maxSync = max(na.omit(Sync)),
            meanSync = mean(na.omit(Sync)),
            medianSync = median(na.omit(Sync)))

syncVals
## isolate high values bteween basins
## could be means are higher between basins as less decay, max is higher within basins with more decay?



