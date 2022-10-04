### matrix regression

library(tidyverse)
library(tidylog)
library(ecodist)

# install.packages("remotes")
# remotes::install_github("reumandc/mms")

# install.packages("PopGenReport")
library(PopGenReport)


# Upload and format synchrony data ----------------------------------------

## functional synchrony

load(file = "output_data/sync/03_sync_traits_CWM_CWV_distances.RData")
head(syncDF)
unique(syncDF$Trait)

funcsync <- syncDF %>%
  rename(Sync = synchrony) %>%
  select(-X) %>%
  filter(Trait == "Tp_pref")

head(funcsync)

## environment synchrony
load(file="output_data/sync/03_sync_temp_distances.RData")
tempsync <- syncDF %>%
  # mutate(SyncType = "TSync")  %>%
  pivot_wider(names_from = Metric, values_from = synchrony) %>%
  filter( !Site_ID1 == Site_ID2)## remove pairs comrised of same sites
head(tempsync)

## remove for now  - fix later if needed
# tempsync <- na.omit(tempsync)

### join functinal synchrony

allsync <- full_join(funcsync, tempsync, by = c("Pair","Site_ID1", "Site_ID2", "Region", "Connectivity", "Euclid_Dist_Meters",
                                               "Similarity", "MeanLat", "MeanLon", "MaxDist")) %>%
  mutate(DistKM = Euclid_Dist_Meters/1000)



# Format sites ------------------------------------------------------------

## some pairs are the same sites but other way around
head(allsync)
length(unique(allsync$Pair))

  
  sites <- allsync %>%
    ungroup() %>%
    # filter(Region == paste(region[r])) %>%
    # dplyr::select(Site_ID1, Site_ID2,Pair, Region) %>% 
    separate(Pair, into =c("Site_ID2a", "Site_ID1a"), remove = F) %>%
    mutate(rev =  paste0(Site_ID1a, ".", Site_ID2a, sep = "")) # #calculate reverse entries
  
  sites$no <- seq_along(sites$Pair) #number the entries
  sites$whererev <- match(sites$rev, sites$Pair) #identify where reversed entries occur
  sites$whererev[sites$whererev>sites$no] <- NA #remove the first of each duplicated pair 
  sites$Pair[!is.na(sites$whererev)] <- sites$rev[!is.na(sites$whererev)] #replace duplicates
  
  all_sites <- unique(sites$Pair) ## pairs with no reverse duplicate
  length(all_sites) ## 208066 
  dim(sites)
  
## subset original df 
  names(sites)
  allsyncx  <- sites %>%
    ungroup() %>%
    dplyr::select(-Site_ID2a, -Site_ID1a, -rev,-no,-whererev) %>%
    distinct(Pair, .keep_all = TRUE)

save(allsyncx, file = "output_data/sync/06_temp_pref_env_dist_no_dupl_pairs.RData")


# Create matrices ---------------------------------------------------------


head(allsyncx)
## temp pref
bio_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "Europe") %>%
  dplyr::select(Site_ID1, Site_ID2,  Sync) %>%
  pivot_wider(names_from = Site_ID2, values_from = Sync) %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows
# nrow(bio_wide)


## get mirror of lower/upper triangle
bio_wide[upper.tri(bio_wide)] <- t(bio_wide)[upper.tri(bio_wide)]
View(bio_wide)##  check visually 
bio_wide <- as.matrix(bio_wide)
bio_wide
write.csv(bio_wide, "output_data/temp_pref_sync_matrix_eu.csv")

## temp sync
temp_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "Europe") %>%
  dplyr::select(Site_ID1, Site_ID2,  annual_avg) %>%
  pivot_wider(names_from = Site_ID2, values_from = annual_avg)  %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows

## get mirror of lower/upper triangle
temp_wide[upper.tri(temp_wide)] <- t(temp_wide)[upper.tri(temp_wide)]
as.data.frame(temp_wide) ##  check visually 
temp_wide <- as.matrix(temp_wide)
View(temp_wide)
## convert NAs to 1

temp_wide[is.na(temp_wide)] <- 1

write.csv(temp_wide, "output_data/temp_sync_matrix_eu.csv")

## diversity
div_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "Europe") %>%
  dplyr::select(Site_ID1, Site_ID2,  diversity) %>%
  pivot_wider(names_from = Site_ID2, values_from = diversity) %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows

## get mirror of lower/upper triangle
div_wide[upper.tri(div_wide)] <- t(div_wide)[upper.tri(div_wide)]
as.data.frame(div_wide) ##  check visually site-site pairs aren't 1?
div_wide <- as.matrix(div_wide)
View(div_wide)


write.csv(div_wide, "output_data/temp_diversity_matrix_eu.csv")


## distance
dist_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "Europe") %>%
  dplyr::select(Site_ID1, Site_ID2,  diversity) %>%
  pivot_wider(names_from = Site_ID2, values_from = diversity)   %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows

## get mirror of lower/upper triangle
dist_wide[upper.tri(dist_wide)] <- t(dist_wide)[upper.tri(dist_wide)]
as.data.frame(dist_wide) ##  check visually site-site pairs aren't 1?
dist_wide <- as.matrix(dist_wide)
View(dist_wide)

write.csv(dist_wide, "output_data/temp_diversity_matrix_eu.csv")

class(dist_wide)

# Matrix regression -------------------------------------------------------

data(landgen)
library(raster)
fric.raster <- readRDS(system.file("extdata","fric.raster.rdata", package="PopGenReport"))
glc <- genleastcost(landgen, fric.raster, "D", NN=4, path="leastcost")
class(glc$gen.mat)
glc$gen.mat
lgrMMRR(glc$gen.mat, glc$cost.mats, glc$eucl.mat, nperm=999)

lgrMMRR(bio_wide, temp_wide, dist_wide, nperm=9)
?lgrMMRR
# matrix regression. so far does not include interactions

MRM(bio_wide ~ dist_wide, nperm = 10)
?MRM
