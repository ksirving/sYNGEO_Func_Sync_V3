## calculate differences in synchrony between all sync and LOO sync

library(tidyr)
library(tidyverse)

## large data so takes a while....

# Single traits -----------------------------------------------------------

## upload
load("output_data/sync/03_sync_data_funcgroup_traitgroup_similarity_euclidean_dist_interpolated.RData")
ssAll <- syncDF
load("output_data/sync/03_sync_data_funcgroup_traitgroup_similarity_euclidean_dist_interpolated_LOO.RData")
ssLOO <- syncDF
rm(syncDF)

head(ssAll)
head(ssLOO)

## workflow: get all sync values, calculate difference, add to LOO df,

## get only traits, pairs and sync from both dfs

ssAll <- ssAll %>% 
  dplyr::select(Trait, TraitGroup, Correlation, Pair) 


ssLOO <- ssLOO %>% 
  dplyr::select(Trait, TraitGroup, Correlation, Pair, YearRemoved) %>%
  rename(CorrelationLOO = Correlation)

head(ssAll)
head(ssLOO)


## join together

ssCor <- left_join(ssAll, ssLOO, by = c("Trait", "TraitGroup", "Pair", "Region"))
head(ssCor)
names(ssCor)
str(ssCor)

head(unique(ssCor$Pair))
tail(unique(ssCor$Pair))

## calculate difference, raw difference and percentage

ssCor_Diff <- ssCor %>%
  # group_by("Trait", "TraitGroup", "Pair", "Region", "YearRemoved") %>%
  mutate(CorrelationLOO = as.numeric(as.character(CorrelationLOO))) %>%
  mutate(SyncDiff = CorrelationLOO - Correlation, SyncDiffPerc = (SyncDiff/Correlation)*100)


head(ssCor_Diff)

## look at ranges

range(na.omit(ssCor_Diff$SyncDiff))

range(na.omit(ssCor_Diff$SyncDiffPerc)) ## silly numbers, stick with absolute difference

### add distances etc back to LOO df

load("output_data/sync/03_sync_data_funcgroup_traitgroup_similarity_euclidean_dist_interpolated_LOO.RData")
ssLOO <- syncDF

ssLOO <- ssLOO %>% 
  rename(CorrelationLOO = Correlation) %>%
  mutate(CorrelationLOO = as.numeric(as.character(CorrelationLOO)))
head(ssLOO)

ssLOO_join <- left_join(ssLOO, ssCor_Diff, by = c("Region", "Trait", "TraitGroup", "CorrelationLOO", "Pair", "YearRemoved"))
head(ssLOO_join)  

ssLOO_join <- ssLOO_join %>% 
  rename(CorrelationALL = Correlation)

save(ssLOO_join, file = "output_data/sync/04_sync_single_traits_LOO.RData")

# Ordination -----------------------------------------------------------

## upload
load("output_data/sync/03_sync_data_ordination_traitgroup_similarity_euclidean_dist_interpolated.RData")
ssAll <- syncDF
load("output_data/sync/03_sync_data_funcgroup_traitgroup_similarity_euclidean_dist_interpolated_ordination_LOO.RData")
ssLOO <- syncDF
rm(syncDF)

head(ssAll)
head(ssLOO)

## workflow: get all sync values, calculate difference, add to LOO df,

## get only traits, pairs and sync from both dfs

ssAll <- ssAll %>% 
  dplyr::select(Axis, TraitGroup, Correlation, Pair) 


ssLOO <- ssLOO %>% 
  dplyr::select(Axis, TraitGroup, Correlation, Pair, YearRemoved) %>%
  rename(CorrelationLOO = Correlation)

head(ssAll)
head(ssLOO)
unique(ssCor$YearRemoved)

## join together

ssCor <- left_join(ssAll, ssLOO, by = c("Axis", "TraitGroup", "Pair", "Region"))
head(ssCor)
names(ssCor)
str(ssCor)

head(unique(ssCor$Pair))
tail(unique(ssCor$Pair))

## calculate difference, raw difference and percentage

ssCor_Diff <- ssCor %>%
  # group_by("Trait", "TraitGroup", "Pair", "Region", "YearRemoved") %>%
  mutate(CorrelationLOO = as.numeric(as.character(CorrelationLOO))) %>%
  mutate(SyncDiff = CorrelationLOO - Correlation, SyncDiffPerc = (SyncDiff/Correlation)*100)


head(ssCor_Diff)

## look at ranges

range(na.omit(ssCor_Diff$SyncDiff))

range(na.omit(ssCor_Diff$SyncDiffPerc))

## add differences back to LOO df

load("output_data/sync/03_sync_data_funcgroup_traitgroup_similarity_euclidean_dist_interpolated_ordination_LOO.RData")
ssLOO <- syncDF

ssLOO <- ssLOO %>% 
  rename(CorrelationLOO = Correlation) %>%
  mutate(CorrelationLOO = as.numeric(as.character(CorrelationLOO)))
head(ssLOO)

ssLOO_join <- left_join(ssLOO, ssCor_Diff, by = c("Region", "Axis", "TraitGroup", "CorrelationLOO", "Pair", "YearRemoved"))
head(ssLOO_join)  

ssLOO_join <- ssLOO_join %>% 
  rename(CorrelationALL = Correlation)

save(ssLOO_join, file = "output_data/sync/04_sync_ordination_LOO.RData")


# Temperature -------------------------------------------------------------

## upload LOO dfs
load(file= "output_data/sync/03_sync_data_temperature_LOO.RData")
ssLOO <- syncDF

## upload all sync

ssAll <- read.csv("output_data/sync/02_temperature_between_all_sites_biogeographic_regions.csv")

head(ssAll)
head(ssLOO)
tail(ssAll$Pair)
## workflow: get all sync values, calculate difference, add to LOO df,

## get only env vars, pairs and sync from both dfs

ssAll <- ssAll %>% 
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  select(-X) %>%
  dplyr::select(env_var, Correlation, Pair, Region) 


ssLOO <- ssLOO %>% 
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  dplyr::select(env_var, Correlation, Pair, YearRemoved, Region) %>%
  rename(CorrelationLOO = Correlation)

head(ssAll)
head(ssLOO)


## join together

ssCor <- left_join(ssAll, ssLOO, by = c("env_var", "Pair", "Region"))
head(ssCor)
names(ssCor)
str(ssCor)
unique(ssCor$YearRemoved)

head(unique(ssCor$Pair))
tail(unique(ssCor$Pair))

## calculate difference, raw difference and percentage

ssCor_Diff <- ssCor %>%
  # group_by("Trait", "TraitGroup", "Pair", "Region", "YearRemoved") %>%
  mutate(CorrelationLOO = as.numeric(as.character(CorrelationLOO))) %>%
  mutate(SyncDiff = CorrelationLOO - Correlation, SyncDiffPerc = (SyncDiff/Correlation)*100)


head(ssCor_Diff)

## look at ranges

range(na.omit(ssCor_Diff$SyncDiff))

range(na.omit(ssCor_Diff$SyncDiffPerc))

## add differences back to LOO df

load(file= "output_data/sync/03_sync_data_temperature_LOO.RData")
ssLOO <- syncDF

ssLOO <- ssLOO %>% 
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  rename(CorrelationLOO = Correlation) %>%
  mutate(CorrelationLOO = as.numeric(as.character(CorrelationLOO)))
head(ssLOO)

ssLOO_join <- left_join(ssLOO, ssCor_Diff, by = c("Region", "env_var", "CorrelationLOO", "Pair", "YearRemoved"))
head(ssLOO_join)  

ssLOO_join <- ssLOO_join %>% 
  rename(CorrelationALL = Correlation)

save(ssLOO_join, file = "output_data/sync/04_sync_temp_LOO.RData")


# Flow --------------------------------------------------------------------

## upload LOO dfs
load(file= "output_data/sync/03_sync_data_flow_LOO.RData")
ssLOO <- syncDF

## upload all sync

ssAll <- read.csv("output_data/sync/02_flow_between_all_sites_biogeographic_regions.csv")

head(ssAll)
head(ssLOO)

## workflow: get all sync values, calculate difference, add to LOO df,

## get only env vars, pairs and sync from both dfs

ssAll <- ssAll %>% 
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  select(-X) %>%
  dplyr::select(env_var, Correlation, Pair, Region) 


ssLOO <- ssLOO %>% 
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  dplyr::select(env_var, Correlation, Pair, YearRemoved, Region) %>%
  rename(CorrelationLOO = Correlation)

head(ssAll)
head(ssLOO)
unique(ssCor$YearRemoved)
# tail(ssAll$Pair)

## join together

ssCor <- left_join(ssAll, ssLOO, by = c("env_var", "Pair", "Region"))
head(ssCor)
names(ssCor)
str(ssCor)

head(unique(ssCor$Pair))
tail(unique(ssCor$Pair))

## calculate difference, raw difference and percentage

ssCor_Diff <- ssCor %>%
  # group_by("Trait", "TraitGroup", "Pair", "Region", "YearRemoved") %>%
  mutate(CorrelationLOO = as.numeric(as.character(CorrelationLOO))) %>%
  mutate(SyncDiff = CorrelationLOO - Correlation, SyncDiffPerc = (SyncDiff/Correlation)*100)


head(ssCor_Diff)

## look at ranges

range(na.omit(ssCor_Diff$SyncDiff))

range(na.omit(ssCor_Diff$SyncDiffPerc))

## add differences back to LOO df

load(file= "output_data/sync/03_sync_data_flow_LOO.RData")
ssLOO <- syncDF

ssLOO <- ssLOO %>% 
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  rename(CorrelationLOO = Correlation) %>%
  mutate(CorrelationLOO = as.numeric(as.character(CorrelationLOO)))
head(ssLOO)

ssLOO_join <- left_join(ssLOO, ssCor_Diff, by = c("Region", "env_var", "CorrelationLOO", "Pair", "YearRemoved"))
head(ssLOO_join)  

ssLOO_join <- ssLOO_join %>% 
  rename(CorrelationALL = Correlation)

save(ssLOO_join, file = "output_data/sync/04_sync_flow_LOO.RData")
