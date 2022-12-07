### table of synchrony

## aim: get synchrony at specific distances e.g. 1km, 10km, 100km
## per trait, per site?

# packages
library(tidylog)
library(tidyverse)

## function to find roots
load(file="functions/root_interpolation_function.Rdata")


# Calculate sync at distances - Single traits -----------------------------------------------------------

## data 

load(file = "output_data/sync/03_sync_data_funcgroup_traitgroup_similarity_euclidean_dist_interpolated.RData") ## syncDF
head(SyncBasin)

SyncBasin <- syncDF %>%
  select(-X) %>%
  mutate(Euclid_Dist_KM = Euclid_Dist_Meters/1000) %>%
  mutate(Euclid_Dist_KM_round = round(Euclid_Dist_KM, digits = 0))


## get mean values per site - lots of NaNs produced, fix!!!!
MeanSync <- SyncBasin %>%
  group_by(Trait, Site_ID2, Connectivity, TraitGroup, Region) %>%
  summarise(Mean_Cor = mean(na.omit(Correlation)), 
            Mean_Dist = mean(na.omit(Euclid_Dist_KM)), 
            Mean_Lat =mean(na.omit(MeanLat)))


head(MeanSync)
max(MeanSync$Mean_Dist)

## all traits

## calculates synchrony at specific distances, i.e., 1km, 10km, 100km & 1000km

names(MeanSync)

## lots of NaNs and -inf produced, fix!!!!
SyncDists <- MeanSync %>%
  group_by(Trait, Region, Connectivity, TraitGroup) %>%
  summarise(Sync_1_Min = min(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 1)),
         Sync_1_Mean = mean(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 1)),
         Sync_1_Max = max(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 1)),
         Sync_10_Min = min(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 10)),
         Sync_10_Mean = mean(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 10)),
         Sync_10_Max = max(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 10)),
         Sync_100_Min = min(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 100)),
         Sync_100_Mean = mean(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 100)),
         Sync_100_Max = max(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 100)),
         Sync_1000_Min = min(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 1000)),
         Sync_1000_Mean = mean(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 1000)),
         Sync_1000_Max = max(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 1000)),
         Sync_2000_Min = min(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 2000)),
         Sync_2000_Mean = mean(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 2000)),
         Sync_2000_Max = max(RootLinearInterpolant(na.omit(Mean_Cor),na.omit(Mean_Dist), 2000)))


### -inf and NaN values occur where that distance doesn't have a sync value 
### i.e., 1km for unconnected sites, or 1000km for connected (same basin) sites                  


# Figure - Single traits --------------------------------------------------

head(SyncDists)

SyncDistsLong <- SyncDists %>%
  pivot_longer(Sync_1_Min:Sync_2000_Max, names_to = "Distances", values_to = "Sync") %>%
  separate(Distances, into = c("s", "km", "SumStat")) %>%
  select(-s) %>%
  mutate(Connectivity= as.factor(Connectivity))

head(SyncDistsLong)

ggplot(SyncDistsLong, aes(y=Sync, x = km, col = Trait))+
  geom_point() +
  facet_grid(rows = vars(Connectivity), cols = vars(Region))



ggplot(SyncDistsLong, aes(y=Sync, x = km, col = SumStat))+
  geom_point(aes(shape=Connectivity)) +
  facet_wrap(~TraitGroup)


ggplot(SyncDistsLong, aes(y=Sync, x = km, col = SumStat))+
  geom_point(aes(shape=Connectivity, size = 1.2)) +
  facet_grid(rows = vars(Connectivity), cols = vars(TraitGroup))
