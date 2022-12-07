### checking values in figures
getwd()

# packages

library(tidyverse)
library(tidylog)
library("easystats")
library(scales)

## directory for figures
out.dir <- "/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V3/Figures/"
# getwd()
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




# Exploration of sites, species & values ----------------------------------

## workflow
## get av and max number of species per time series
## get av relative abundance per time series
## compare no species ~ rel abundance
## compare diversity value with single species sites (or low number of species)
## compare no species ~ diversity & distance
## calculate diversity as coef of variation


load(file="output_data/sync/01_fish_abundances.RData")
names(fish_ab_rel_int)

## number of species and av rel abundance per TS
fish_ab_rel_int_sp <- fish_ab_rel_int %>%
  ungroup() %>%
  group_by(site_year) %>% 
  mutate(NumberSpeciesSiteYear = length(unique(Species))) %>% ## number of species per site year
  ungroup() %>%
  group_by(SiteID) %>%
  mutate(RelAbundTS = mean(RelAbundanceSiteYear)) %>%## average relative abundance over TS
  ungroup() %>%
  group_by(SiteID) %>%
  mutate(NumberSpeciesTS = mean(NumberSpeciesSiteYear)) %>% #%>% ## number of species per site
  dplyr::select(SiteID, RelAbundTS, NumberSpeciesTS) %>%
  distinct()

cor(fish_ab_rel_int_sp$TotalAbundanceSite, fish_ab_rel_int_sp$NumberSpeciesTS) ## 0.63



## difference in no of species/rel abundance between site pairs
?expand.grid
## get site pairs list
cc <- expand.grid(colnames(fishRel), colnames(fishRel), KEEP.OUT.ATTRS = FALSE)
head(cc)
## define sites
sites <- fish_ab_rel_int_sp$SiteID
s=2
sites

pairsx <- NULL

for(s in 1:length(sites)) {
   ## define first site in pair
  site <- sites[s]
  ## get values
  vals1 <- filter(fish_ab_rel_int_sp, SiteID == site)
  vals1
  ## add values to site1 
  pairs <- filter(cc, Var1 == site) %>%
    mutate(RelAbundTS1 = vals1$RelAbundTS,
           NumberSpeciesTS1 = vals1$NumberSpeciesTS) %>%
    rename(SiteID = Var2)
  ## define paired sites
  site2 <- pairs$SiteID
  ## get values for paired sites
  vals2 <- filter(fish_ab_rel_int_sp, SiteID %in% site2)
  
  ## join values for all pairs and get differences
  pairs2 <- full_join(pairs, vals2, by = "SiteID") %>%
    mutate(RelAbundDiff = RelAbundTS1-RelAbundTS,
           NumSpeciesDiff = NumberSpeciesTS1 - NumberSpeciesTS) %>%
    distinct()
 
  head(pairsx)
  
  pairsx <- bind_rows(pairsx, pairs2)
  
}

head(pairsx)
## remove reverse duplicates
sites <- pairsx %>%
  ungroup() %>%
  mutate(Pair = paste0(Var1, ".", SiteID)) %>%
  # dplyr::select(Site_ID1, Site_ID2,Pair, Region) %>% 
  separate(Pair, into =c("Site_ID2a", "Site_ID1a"), remove = F) %>%
  mutate(rev =  paste0(Site_ID1a, ".", Site_ID2a, sep = "")) # #calculate reverse entries

head(sites)

sites$no <- seq_along(sites$Pair) #number the entries
sites$whererev <- match(sites$rev, sites$Pair) #identify where reversed entries occur
sites$whererev[sites$whererev>sites$no] <- NA #remove the first of each duplicated pair 
sites$Pair[!is.na(sites$whererev)] <- sites$rev[!is.na(sites$whererev)] #replace duplicates

all_sites <- unique(sites$Pair) ## pairs with no reverse duplicate
length(all_sites) ## 294528
dim(sites)

## subset original df 
names(sites)
allvals  <- sites %>%
  ungroup() %>%
  dplyr::select(-Site_ID2a, -Site_ID1a, -rev,-no,-whererev) %>%
  distinct(Pair, .keep_all = TRUE)

head(allvals)

save(allvals, file = "output_data/sync/05a_diff_species_rel_abundance.RData")

cor(allvals$RelAbundDiff, allvals$NumSpeciesDiff) ##  -0.52

ggplot(allvals, aes(x= RelAbundDiff, y = NumSpeciesDiff)) +
  geom_point()

#### compare with diversity etc

load(file = "output_data/sync/04_temp_pref_env_dist_no_dupl_pairs.RData")



# allsyncx <- na.omit(allsyncx)

head(allsyncx)

### join with difference

allvals <- allvals %>%
  select(Pair, RelAbundDiff, NumSpeciesDiff)

df <- inner_join(allsyncx, allvals, by = "Pair")
head(df)
names(df)

cordf <- df %>%
  select(distance, diversity, MeanLat, MeanLon, annual_avg, RelAbundDiff, NumSpeciesDiff)

cor(cordf)
write.csv(cordf, "output_data/05a_correlations_Species_rel_abund.csv")

cd <- cor(cordf)
write.csv(cd, "output_data/05a_correlations_Species_rel_abund_mat.csv")



### sites with only one species

## check sites with 1 species
## filter on 1 species
singSpecies <- fish_ab_rel_int %>%
  ungroup() %>%
  group_by(site_year) %>% 
  mutate(NumberSpeciesSiteYear = length(unique(Species))) %>%
  filter(NumberSpeciesSiteYear %in% c(1))

singSpeciessites <- unique(singSpecies$SiteID)
singSpeciessites

## get site pairs list
ss <- expand.grid(singSpeciessites, singSpeciessites, KEEP.OUT.ATTRS = FALSE)

ss$Pair <- paste0(ss$Var1, ".", ss$Var2)
head(ss)


ss_sites <- unique(ss$Pair)
ss_sites
length(unique(ss_sites)) ## 2401

singSpeciesdf <- df %>%
  filter(Pair %in% ss_sites)

sum(!ss_sites %in% singSpeciesdf$Pair) ## 638
sum(ss_sites %in% df$Pair) ## 638
ind <- !ss_sites %in% singSpeciesdf$Pair
otherSites <- ss_sites[ind]

singSpeciesdfmis <- df %>%
  filter(Pair %in% otherSites)
singSpeciesdfmis ## not in dataset???
## 1176 are duplicates, 638 are 1 species - where are the remaining 587? they are all sync NAs or duplicates
head(singSpeciesdf) ## nothing here, are these the NAs??


## check in unforamtted df

load(file = "output_data/sync/04_temp_pref_env_dist_ALL.RData")

singSpeciesorig <- allsync %>%
  filter(Pair %in% ss_sites)

singSpeciesorig <- allsync %>%
  filter(Pair %in% otherSites)

## some NAs in sync. find and remove - check in sync code later!!!!!
ind <- which(is.na(df))
test <- df[ind,]
# allsyncx <- allsyncx[-ind,]
sum(is.na(test))

naPairs <- unique(df$Pair)

# check cwm ~ cwv 

weiVals <- read.csv("output_data/01_trt_single_traits_interpolated_cwm_cmv.csv")

head(weiVals)
cor(weiVals$CWM, weiVals$CWV) ## 0.49

## check one species sites in cmv

weiValsss <- weiVals %>%
  filter(SiteID %in% singSpeciessites)
  
length(unique(weiValsss$SiteID))

## check sites with 2 species
## filter on 2 species
singSpecies <- fish_ab_rel_int %>%
  ungroup() %>%
  group_by(site_year) %>% 
  mutate(NumberSpeciesSiteYear = length(unique(Species))) %>%
  filter(NumberSpeciesSiteYear %in% c(2))

singSpeciessites <- unique(singSpecies$SiteID)
singSpeciessites ## 54

## get site pairs list
ss <- expand.grid(singSpeciessites, singSpeciessites, KEEP.OUT.ATTRS = FALSE)

ss$Pair <- paste0(ss$Var1, ".", ss$Var2)
head(ss)


ss_sites <- unique(ss$Pair)
ss_sites
length(unique(ss_sites)) ## 2916

singSpeciesdf <- df %>%
  filter(Pair %in% ss_sites)
singSpeciesdf 
mean(singSpeciesdf$diversity)

singSpeciesdfmis <- df %>%
  filter(Pair %in% otherSites)
names(fish_ab_rel_int)
whichSpecies <- fish_ab_rel_int  %>%
  filter(SiteID %in% singSpeciessites) %>%
  select(Species, site_year, SiteID, Abundance, RelAbundanceSiteYear, Country) 

## some NAs in country, fill and re join for boxplot
nasF <- whichSpecies %>%
  filter(Country == "FRA") #%>%
  # mutate(Country = "FRA")
  
fra <- unique(nasF$SiteID)
fra
nasS <- whichSpecies %>%
  filter(Country == "SWE") #%>%
  # mutate(Country = "SWE")

swe <- unique(nasS$SiteID)
swe
whichSpeciesF <- whichSpecies %>%
  filter(SiteID %in% fra) %>%
  mutate(Country = "FRA")

whichSpeciesS <- whichSpecies %>%
  filter(SiteID %in% swe) %>%
  mutate(Country = "SWE")

whichSpecies <- bind_rows(whichSpeciesS, whichSpeciesF)

length(unique(whichSpecies$SiteID)) ## 54

## boxplot of sites with 2 species rel abundance
B1 <- ggplot(whichSpecies, aes(x=Species, y = RelAbundanceSiteYear)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Country)

B1

file.name1 <- paste0(out.dir, "Boxplot_rel_abund_2_species.jpg")
ggsave(B1, filename=file.name1, dpi=300, height=5, width=6)


## boxplot of sites with 2 species, diversity & distance

singSpeciesdfLong <- singSpeciesdf %>%
  pivot_longer(c(Sync:diversity), names_to="Variable", values_to = "Values")

B2 <- ggplot(singSpeciesdfLong, aes(x=Variable, y = Values)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

B2

file.name1 <- paste0(out.dir, "Values_2_species.jpg")
ggsave(B2, filename=file.name1, dpi=300, height=5, width=6)

