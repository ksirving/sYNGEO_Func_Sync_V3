## add dummy variabkle within/between basins

library(sp)
library(raster)
library(ggplot2)
library(dplyr)
library(tidyverse)

getwd()

## directory for figures
out.dir <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V2/Figures/"

## upload fish abundance and site data
originaldata <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")
head(originaldata)

## take only sites and basins
sites <- originaldata %>%
  dplyr::select(SiteID, HydroBasin) %>%
  distinct()

## all site synchrony

sync <- read.csv("output_data/sync/02_between_all_sites_single_traits_CWM_CWV.csv")
head(sync)

## join first site basin
site1 <- left_join(sync, sites, by = c("Site_ID1" = "SiteID"))
head(site1)

## change names so we know it's for site 1
site1 <- site1 %>%
  select(-X) %>%
  rename(HydroBasin1 = HydroBasin)

## join site 2

site2 <- left_join(site1, sites, by = c("Site_ID2" = "SiteID"))
head(site2)

## change names so we know it's for site 2
site2 <- site2 %>%
  # select(-X) %>%
  rename(HydroBasin2 = HydroBasin)

head(site2)


## create dummy variable. If hydrobasin 1 matches hydrobasin 2 = within site (1), if not then between sites (0)

con <- site2 %>%
  mutate(Connectivity = ifelse(HydroBasin1 == HydroBasin2, 1, 0))

head(con)

write.csv(con, "output_data/sync/03_func_sync_connectivity.csv")

rm(con)
# calculate distance ------------------------------------------------------


sync <- read.csv("output_data/sync/03_func_sync_connectivity.csv")
head(sync)

## define pairs
syncsites <- sync %>%
  dplyr::select(Site_ID1, Site_ID2, Pair, Connectivity) %>%
  # filter(TraitGroup == "FoodAquisition") %>%
  # tidyr::pivot_wider(names_from = Trait, values_from = Correlation) %>%
  # mutate(Pair = paste(Site_ID1, ".", Site_ID2, sep="")) %>%
  mutate(Euclid_Dist_Meters = 0, Similarity = 0, MeanLat = 0, MeanLon = 0) %>%
  # dplyr::select(-FeedingGroup, ) %>%
  distinct()

head(syncsites)

pairs <- unique(syncsites$Pair)
tail(pairs)
length(pairs) ##  207299

## site coords

fish_ab <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")
head(fish_ab)

## get coords
SiteCoords <- fish_ab %>%
  dplyr::select(SiteID, Latitude, Longitude) %>%
  distinct()

## loop over pairs - takes a long ass time, go do something else for a bit...

for(p in 1:length(pairs)) {

  ## get pair from sync data

  pairx <- syncsites %>%
    filter(Pair == pairs[p])
  # pairx
  ## define sites
  S1 <- pairx$Site_ID1
  S2 <- pairx$Site_ID2

  # S1

  ## get coords for each site
  CoordsS1 <- SiteCoords %>%
    filter(SiteID == S1) %>%
    dplyr::select(Longitude, Latitude, SiteID)

  CoordsS2 <- SiteCoords %>%
    filter(SiteID == S2) %>%
    dplyr::select(Longitude, Latitude, SiteID)

  sp::coordinates(CoordsS1) <- c("Longitude", "Latitude")
  sp::coordinates(CoordsS2) <- c("Longitude", "Latitude")

  #Make a distance matrix
  dst <- pointDistance(CoordsS1,CoordsS2, lonlat=TRUE)
  # str(dst)
  # get mean latitude/longitude
  MeanLat <- (CoordsS1$Latitude+CoordsS2$Latitude)/2
  MeanLon <- (CoordsS1$Longitude+CoordsS2$Longitude)/2

  ## add to dataframe
  syncsites[p,5] <- dst
  syncsites[p,7] <- MeanLat
  syncsites[p,8] <- MeanLon
 # head(syncsites)

}

head(syncsites)

save(syncsites, file = "output_data/sync/03_all_pair_distances.RData")



# Combine with single trait DF --------------------------------------------

sync <- read.csv("output_data/sync/02_between_all_sites_single_traits_CWM_CWV.csv")
head(sync)

load(file = "output_data/sync/03_all_pair_distances.RData") ## syncsites

## make Pair column
# sync <- sync %>%
#   unite(Pair, Site_ID1:Site_ID2,sep = ".", remove=F)

## take only distance columns
sync_sub <- syncsites %>%
  dplyr::select(Pair:MeanLon)

## join
all_sync <- left_join(sync, sync_sub, by = "Pair")
head(all_sync)
## convert to similarities

syncDF <- all_sync %>%
  group_by(Region, Trait) %>%
  mutate(MaxDist = max(Euclid_Dist_Meters)) %>%
  mutate(Similarity = 1-(Euclid_Dist_Meters/MaxDist))

head(syncDF)


save(syncDF, file = "output_data/sync/03_sync_traits_CWM_CWV_distances.RData")


# Join with temp and flow -------------------------------------------------

sync <- read.csv( "output_data/sync/02_between_all_sites_temp_synchrony.csv")
head(sync)

load(file = "output_data/sync/03_all_pair_distances.RData") ## syncsites

## make Pair column
sync <- sync %>%
  select(-X) %>%
  unite(Pair, Site_ID1:Site_ID2,sep = ".", remove=F)

## take only distance columns
sync_sub <- syncsites %>%
  dplyr::select(Pair:MeanLon)

## join
all_sync <- left_join(sync, sync_sub, by = "Pair")
head(all_sync)
## convert to similarities

syncDF <- all_sync %>%
  group_by(Region, Metric) %>%
  mutate(MaxDist = max(Euclid_Dist_Meters)) %>%
  mutate(Similarity = 1-(Euclid_Dist_Meters/MaxDist))

head(syncDF)


save(syncDF, file = "output_data/sync/03_sync_temp_distances.RData")

##### Flow

sync <- read.csv( "output_data/sync/02_between_all_sites_flow_synchrony.csv")
head(sync)

load(file = "output_data/sync/03_all_pair_distances.RData") ## syncsites

## make Pair column
sync <- sync %>%
  select(-X) %>%
  unite(Pair, Site_ID1:Site_ID2,sep = ".", remove=F)

## take only distance columns
sync_sub <- syncsites %>%
  dplyr::select(Pair:MeanLon)

## join
all_sync <- left_join(sync, sync_sub, by = "Pair")
head(all_sync)
## convert to similarities

syncDF <- all_sync %>%
  group_by(Region, Metric) %>%
  mutate(MaxDist = max(Euclid_Dist_Meters)) %>%
  mutate(Similarity = 1-(Euclid_Dist_Meters/MaxDist))

head(syncDF)


save(syncDF, file = "output_data/sync/03_sync_flow_distances.RData")


# watercourse V euclid distance -------------------------------------------

library(scales)
## euclid distance
load(file = "output_data/sync/03_all_pair_distances.RData") ## syncsites
head(syncsites)

## format euclid - remove sites not in same basin
syncsites <- syncsites %>%
  mutate(DistMetersEuclid = Euclid_Dist_Meters) %>%
  filter(Connectivity == 1) %>%
  select(Pair, DistMetersEuclid)

## water course distance
watersites <- read.csv2("input_data/sites/Sites_DistancesRiverATLAS_280621.csv")
head(watersites)

## format water course data

watersites <- watersites %>%
  mutate(Pair = paste(SiteID_Orig, SiteID_Dest, sep =".")) %>%
  mutate(DistMetersWater = as.numeric(TotLong_Meters)) %>%
  select(Pair, DistMetersWater, BioRealm, HydroBasin, Country)

head(watersites)
head(syncsites)

length(syncsites$Pair) # 8459
length(watersites$Pair) # 17221
sum(watersites$Pair %in% syncsites$Pair) ## ~200 site pairs missing - could be related to filtering of species/traits/sites

all_sites <- left_join(syncsites, watersites, by = "Pair")
head(all_sites)

str(all_sites)

cor(all_sites$DistMetersEuclid, all_sites$DistMetersWater, use = "complete.obs") ## 0.93

## make data long for plot

# all_sites_long <- all_sites %>%
#   pivot_longer(DistMetersEuclid:DistMetersWater, names_to = "Type", values_to = "Meters")
# 
# head(all_sites_long)

t1 <- ggplot(all_sites, aes(x=DistMetersWater/1000, y = DistMetersEuclid/1000)) +
  geom_point() +
  scale_y_log10(name="Log Eucliean Distance (km)", labels = comma) +
  scale_x_log10(name="Log Water Course Distance (km)", labels = comma) 

file.name1 <- paste0(out.dir, "watercourse_v_euclid.jpg")
ggsave(t1, filename=file.name1, dpi=300, height=5, width=6)
t1

ggplot(all_sites, aes(x=DistMetersWater/1000, y = DistMetersEuclid/1000)) +
  geom_point() +
  scale_y_continuous(name="Eucliean Distance (km)", labels = comma, limits = c(0, 3000)) +
  scale_x_continuous(name="Water Course Distance (km)", labels = comma, limits = c(0, 3000))

# Checking weird sites ----------------------------------------------------

## problem is some euclid distance is longer than water course - water course should always be longer

## get ratio of watercourse/euclid distance

head(all_sites)

ratios <- all_sites %>%
  mutate(Ratio = DistMetersEuclid/DistMetersWater)

head(ratios)

ratios_big <- ratios %>%
  filter(Ratio > 1) %>%
  separate(Pair, into = c("Site_ID1", "Site_ID2"), remove = F)

dim(ratios_big) ## 242

write.csv(ratios_big, "output_data/03_watercourse_smaller_than_euclid_dist.csv")

head(ratios_big)
## add coords

head(SiteCoords)

pairs <- ratios_big$Pair
pairs



p=1
  
  ## get pair from sync data
  
  pairx <- ratios_big %>%
    filter(Pair == pairs[p])
  pairx
  ## define sites
  S1 <- pairx$Site_ID1
  S2 <- pairx$Site_ID2
  
  S1
  
  ## get coords for each site
  CoordsS1 <- SiteCoords %>%
    filter(SiteID == S1) %>%
    dplyr::select(Longitude, Latitude, SiteID)
  
  CoordsS1
  
  CoordsS2 <- SiteCoords %>%
    filter(SiteID == S2) %>%
    dplyr::select(Longitude, Latitude, SiteID)
  
  CoordsS2
  
  sp::coordinates(CoordsS1) <- c("Longitude", "Latitude")
  sp::coordinates(CoordsS2) <- c("Longitude", "Latitude")
  
  #Make a distance matrix
  dst <- pointDistance(CoordsS1,CoordsS2, lonlat=TRUE)
  dst
  # str(dst)
  # get mean latitude/longitude
  MeanLat <- (CoordsS1$Latitude+CoordsS2$Latitude)/2
  MeanLon <- (CoordsS1$Longitude+CoordsS2$Longitude)/2
  
  ## add to dataframe
  syncsites[p,8] <- dst
  syncsites[p,10] <- MeanLat
  syncsites[p,11] <- MeanLon
  head(syncsites)
  



# Add water course distance to synchrony ----------------------------------

head(watersites)

## convert to km

watersites <- watersites %>%
  mutate(WCDistkm = DistMetersWater/1000)
  
## upload 
load(file = "output_data/sync/03_sync_traits_CWM_CWV_distances.RData")

head(syncDF)

unique(syncDF$Trait)
syncDF <- left_join(syncDF, watersites, by = "Pair") 

syncDF <- syncDF %>%
  filter(Connectivity == 1, Trait %in% c("AVG_MXL", "AVG_FECUND", "Tp_pref", "Q_pref")) %>%
  mutate(Euclid_Dist_KM = Euclid_Dist_Meters/1000 ) %>%
  pivot_longer(c(Euclid_Dist_KM, WCDistkm), names_to = "Dist_Type", values_to = "KM")


sm1a <- ggplot(syncDF, aes(x=KM, y=synchrony, color = Dist_Type)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Trait) +
  scale_color_discrete(name = "Distance Type", labels = c("Euclidean", "Water Course")) +
  scale_y_continuous(name="Synchrony") +
  scale_x_log10(name="Log Distance (km)", labels = comma) 
sm1a


file.name1 <- paste0(out.dir, "Euclid_v_waterCourse_distance_4_traits.jpg")
ggsave(sm1a, filename=file.name1, dpi=300, height=5, width=6)








