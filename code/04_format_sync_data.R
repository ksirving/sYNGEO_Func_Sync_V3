## format for figures and stats

# packages

library(tidyverse)
library(tidylog)
library("easystats")
library(scales)
getwd()
## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V2/Figures/"


# Temp pref synchrony -----------------------------------------------------


load(file = "output_data/sync/03_sync_traits_CWM_CWV_distances.RData")

funcsync <- syncDF %>%
  rename(Sync = synchrony) %>%
  dplyr::select(-X) 

head(funcsync)

## checking missing pairs script 5a

singSpeciesorig <- funcsync %>%
  filter(Pair %in% otherSites) ## all sync NA values
dim(singSpeciesorig) ##  1763

# Temp synchrony ----------------------------------------------------------

load(file="output_data/sync/03_sync_temp_distances.RData")

tempsync <- syncDF %>%
  # mutate(SyncType = "TSync")  %>%
  pivot_wider(names_from = Metric, values_from = synchrony) %>%
  filter( !Site_ID1 == Site_ID2)## remove pairs comrised of same sites
head(tempsync)

## checking missing pairs script 5a

singSpeciesorig <- tempsync %>%
  filter(Pair %in% otherSites) 
dim(singSpeciesorig) ##  1763

# Join --------------------------------------------------------------------


allsync <- full_join(funcsync, tempsync, by = c("Pair","Site_ID1", "Site_ID2", "Region", "Connectivity", "Euclid_Dist_Meters",
                                                "Similarity", "MeanLat", "MeanLon", "MaxDist")) %>%
  mutate(DistKM = Euclid_Dist_Meters/1000)

save(allsync, file = "output_data/sync/04_temp_pref_env_dist_ALL.RData")
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

## checking missing pairs script 5a

singSpeciesorig <- sites %>%
  filter(Pair %in% otherSites) ## all the duplicates



## subset original df 
names(sites)
allsyncx  <- sites %>%
  ungroup() %>%
  dplyr::select(-Site_ID2a, -Site_ID1a, -rev,-no,-whererev) %>%
  distinct(Pair, .keep_all = TRUE)

save(allsyncx, file = "output_data/sync/04_temp_pref_env_dist_no_dupl_pairs.RData")


