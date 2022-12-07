### interpolated abundances on trait ordination

library(FD)
library(vegan)
library(ape)
library(plyr)
library(reshape2)
library(tidyr)
library(dplyr)
library(ade4)
library(cluster)
library(ggplot2)
library(Amelia)
library("StatMatch")

setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V2")
getwd()

## directory for figures
out.dir <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V2/Figures/"

# Traits ------------------------------------------------------------------
  
## upload traits
trt1 <- read.csv("input_data/Bio/matSpecies_imputedCloseTaxa.csv")

head(trt1)
dim(trt1) ## 272, 21 traits

## remove CTMax and species for morpho also
## count NAs in each row - nas are missing traits
trt <- trt1 %>%
  select(-CTmax, -Species_used_for_morpho, -imput) %>%
  mutate(number_nas = rowSums(is.na(trt1)))

rare_species <- trt %>%
  filter(number_nas >= 2) ## number of missing traits

RSp <- rare_species$Species
RSp

head(trt)

str(cor_trt)

cor_trt <- trt %>%
  select(-Species, - number_nas, -AVG_RGUILD, -AVG_Troph)

cor_mat <- cor(cor_trt, use = "complete.obs")

?cor

write.csv(cor_mat, "output_data/01b_trait_correlation.csv")


# log traits
trt$AVG_MXL<-(trt$AVG_MXL)/10
trt$AVG_MXL<-log(trt$AVG_MXL+1)
trt$AVG_LMAT<-log(trt$AVG_LMAT+1)
trt$AVG_AGEMAT<-log(trt$AVG_AGEMAT+1)
trt$AVG_LONGY<-log(trt$AVG_LONGY+1)
trt$AVG_FECUND<-log(trt$AVG_FECUND+1)
trt$AVG_EGGSIZE<-log(trt$AVG_EGGSIZE+1)
trt$Q_pref <- log(trt$Q_pref+1)
trt$JlHd <- log(trt$JlHd+1)
trt$BlBd <- log(trt$BlBd+1)
trt$CFdCPd <- log(trt$CFdCPd+1)
trt$Length <- log(trt$Length+1)

## match lengths
trt$AVG_MXL
trt$Length
## most slightly out
sum(is.na(trt$AVG_MXL)) ## 1
sum(is.na(trt$Length)) ## 26

cor(trt$AVG_MXL, trt$Length, use = "complete.obs") ## 0.99

## keep both and remove length later (no Nas)

## rank reproductive guild

trt <- trt %>%
  mutate(AVG_RGUILD_ORD = NA) %>%
  mutate(AVG_RGUILD_ORD = replace(AVG_RGUILD_ORD, AVG_RGUILD == "Bearer", 5)) %>%
  mutate(AVG_RGUILD_ORD = replace(AVG_RGUILD_ORD, AVG_RGUILD == "G_NS", 4)) %>%
  mutate(AVG_RGUILD_ORD = replace(AVG_RGUILD_ORD, AVG_RGUILD == "G_SC", 3)) %>%
  mutate(AVG_RGUILD_ORD = replace(AVG_RGUILD_ORD, AVG_RGUILD == "NG_BH", 2)) %>%
  mutate(AVG_RGUILD_ORD = replace(AVG_RGUILD_ORD, AVG_RGUILD == "NG_OS", 1)) %>%
  select(-AVG_RGUILD)

head(trt)
str(trt)

trt$AVG_RGUILD_ORD =  as.ordered(trt$AVG_RGUILD_ORD)

## trait correlation
str(trt_cor)
trt_cor <- trt %>%
  select(-AVG_Troph, -AVG_RGUILD_ORD, -Species)

head(trt_cor)
# ?cor
cor_mat <- cor(trt_cor,  use = "pairwise.complete.obs")

## sepatate into effect/response/both and categories

# EdHd
# MoBd
# JlHd
# EhBd
# BlBd
## food aquisition - effect

# HdBd
# PFiBd
# PFlBl
# Length
# AVG_MXL
## swimming ability - both

# AVG_LMAT
# AVG_AGEMAT
# AVG_LONGY
# AVG_FECUND
# AVG_EGGSIZE
# AVG_RGUILD
## reproduction - response

# tp_pref
# Q_pref
# env preferences - response

# AVG_Troph
#  trophic level - grouping

head(trt)
dim(trt)

write.csv(trt, "output_data/01_traits_corrected.csv")


# Abundances --------------------------------------------------------------

fish_ab <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")

head(fish_ab)
str(fish_ab)
## remove basins - Sweden 2080030650 and 2080031160
## keep only origin Ohio and LTRM in mississippi
## change to relative abundance
## missing trait values - remove fish with less than 2 traits (check that it's less than 5%)
## remove basins - Sweden 2080030650 and 2080031160
## keep only origin Ohio and LTRM in mississippi
basins_remove <- c(2080031160, 2080030650)
origin_remove <- c("Tennessee", "ltr")

fish_ab <- fish_ab %>%
  filter(!HydroBasin %in% basins_remove, !ORIGIN %in% origin_remove) 

sites <- fish_ab %>%
  select(SiteID, Latitude, Longitude, BioRealm,  HydroBasin, Country) %>%
  distinct()

## count how many years per site
tally <- fish_ab %>%
  group_by(SiteID) %>%
  select(Year, SiteID) %>%
  distinct() %>%
  tally()

head(tally) ## some have less than 8 years - remove

remove_sites <- tally %>%
  filter(n < 8) 

remove_sites <- remove_sites$SiteID
remove_sites ## sites to be removed

## remove from main DF
fish_ab  <- fish_ab %>%
  filter(!SiteID %in% remove_sites)
write_csv(fish_ab, here("input_data", "Bio", "fish_ab.csv"))

## look at abundance per site over all years
## define years
years <- c("2004" ,"2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013")

## define sites
sites <- unique(fish_ab$SiteID)
s=21
s
## define df 
mediansx <- NULL
names(fish_ab)

for(s in 1:length(sites)) {
  
  ### change all NAs to 0
  
  FishInt <- fish_ab %>%
    filter(SiteID == sites[s]) %>%
    select(-c(TimeSeriesID,SurveyID, Quarter, UnitAbundance,SourceID, Protocol:ORIGIN )) %>% 
    pivot_wider(names_from = "Year", values_from = "Abundance") %>%
    replace(is.na(.), 0)
  
  FishInt
  ### change all NAs to 0
  
  site_years <- sort(names(FishInt[, 3:length(colnames(FishInt))]))
  site_years
  
  medians <- FishInt %>%
    rowwise() %>% 
    mutate(med = median(c_across(where(is.numeric)), na.rm = TRUE))
  
  
  if (sum(years %in% site_years == F)==1) {
    ## if 1 year is missing then interpolate
    missing_yearx <- which(years %in% site_years == F)
    missing_year <- years[missing_yearx]
    
    colnames(medians)[colnames(medians) == "med"] <- paste(missing_year)
    
  } else if (sum(years %in% site_years == F)==2) {
    ## if 2 year is missing then interpolate
    missing_yearx <- which(years %in% site_years == F)
    missing_year <- years[missing_yearx]
    
    medians$med2 <- medians$med
    
    colnames(medians)[colnames(medians) == "med"] <- paste(missing_year)[1]
    
    colnames(medians)[colnames(medians) == "med2"] <- paste(missing_year)[2]
    
  } else {
    medians <-  medians %>% select(-med)
    
  }
  
  mediansx <- bind_rows(mediansx, medians)
  
}

head(mediansx)
mediansx <- mediansx %>%
  pivot_longer(cols = `2009`:`2008`, names_to = "Year", values_to = "Abundance") %>%
  unite("site_year", c(SiteID, Year), sep="_", remove = F)
dim(mediansx)

## add back to main DF
head(fish_ab)
dim(fish_ab)
fish_ab_sub <- fish_ab %>%
  unite("site_year", c(SiteID, Year), sep="_", remove = F) %>%
  select(-Species, -Abundance, -Year) ## remove original abundances

## join new abundances
fish_ab_int <- full_join(mediansx, fish_ab_sub, by = "site_year")
fish_ab_int <- fish_ab_int %>%
  rename(SiteID = SiteID.x) %>%
  select(-SiteID.y)

names(fish_ab_int)
head(fish_ab_int)

## change to relative abundance - site/year
fish_ab_rel_int <- fish_ab_int %>%
  group_by(SiteID, Year) %>%
  mutate(TotalAbundanceSiteYear = sum(Abundance)) %>% ## total species abundance at each site/year
  ungroup %>%
  group_by(SiteID) %>%
  mutate(TotalAbundanceSite = sum(Abundance)) %>% ## total species abundance at each site
  ungroup() %>%
  group_by(Species, SiteID) %>%
  mutate(TotalAbundanceSiteSpecies = sum(Abundance)) ## abundance of each species at each site


fish_ab_rel_int <- fish_ab_rel_int %>%
  ungroup() %>%
  mutate(RelAbundanceSiteYear = (Abundance/TotalAbundanceSiteYear)*100) %>%
  mutate(RelAbundanceSite = (TotalAbundanceSiteSpecies/TotalAbundanceSite)*100)

## remove rare species (species with 2 or less traits) from main df

fish_ab_rel_int <- fish_ab_rel_int %>%
  filter(!Species %in% RSp)

RSp
## filter trt to same species as fish df

fish_sp <- unique(fish_ab_rel_int$Species)
fish_sp

trt <- trt %>%
  filter(Species %in% fish_sp)
str(trt)

## match species 
# check the matchin of spp
setdiff(trt$Species, fish_ab_rel_int$Species)
setdiff(fish_ab_rel_int$Species, trt$Species) # all matched!!!

# make sure basin is a factor
fish_ab_rel_int$HydroBasin<-as.factor(fish_ab_rel_int$HydroBasin)

### define traits groups
resource <- trt %>%
  select(Species, EdHd:BlBd)

dispersal <- trt %>%
  select(Species, AVG_MXL, HdBd:CFdCPd)

reproduction <- trt %>%
  select(Species, AVG_LMAT:AVG_EGGSIZE, AVG_RGUILD_ORD)

envPref <- trt %>%
  select(Species, Tp_pref:Q_pref)


### format fish abundances
fish_ab2 <- fish_ab_rel_int

fish_ab2$site_seas_year<-paste(fish_ab2$SiteID, fish_ab2$Month,fish_ab2$Year, sep="_")

# convert to wide format

fish_mat2<-dcast(fish_ab2, site_seas_year  ~ Species, value.var="RelAbundanceSiteYear", fun.aggregate = sum)
names(fish_ab2)
# add  columns of year, site and seasons to the fish_mat2 matrix using "colsplit"
# Some sites have NAs on season; season values with NAs are pasted but turned into character#
fish_mat3<-cbind(colsplit(fish_mat2$site_seas_year,"_", c("site","month", "year")),fish_mat2)

# add column site_year (if duplicates are here means two seasons are covered; then we can aggregate these within a site ##
fish_mat3$site_year<-paste(fish_mat3$site, fish_mat3$year, sep="_")

# find sites with two seasons
which((duplicated(fish_mat3$site_year)==TRUE)) ### no sites sampled twice in a year

#assign row names to the final fish abund matrix
row.names(fish_mat3)<-fish_mat3$site_year
sum(is.na(fish_mat3)) ## 1262
write_csv(fish_mat3, here("input_data", "Bio", "fish_mat3.csv"))


### merge with traits - food

## functcomp requires spp and site names to be row names and column names etc, and the df should not contain other info ##
resource<- resource[order(resource$Species),] # sort species names in the trait df (they should match the fish matrix)
row.names(resource)<-resource$Species
resource$Species<-NULL
# om_resource$number_nas<-NULL

head(resource) ## traits per species
head(fish_mat3) ### abundance of species per site year
names(fish_mat3)
dim(fish_mat3)

# # create a "clean" df (called fish for traits "fish_fortr") with fish abundance for the functcomp command (weighting traits by spp relative abund)
fish_fortrt<-fish_mat3[,c(5:202)]
row.names(fish_fortrt)<-fish_mat3$site_year
dim(fish_fortrt)
dim(resource)
## computes the functional composition of functional traits, by community weighted mean
res_matrix<-functcomp(resource, as.matrix(fish_fortrt), CWM.type = "dom")  
head(res_matrix)

res_matrix <- res_matrix %>%
  mutate(SiteYear = rownames(res_matrix)) %>%
  pivot_longer(EdHd:BlBd, names_to="Trait", values_to= "Values") %>%
  mutate(TraitGroup = "FoodAquisition", TraitType = "Effect") 


## functcomp requires spp and site names to be row names and column names etc, and the df should not contain other info ##
dispersal<- dispersal[order(dispersal$Species),] # sort species names in the trait df (they should match the fish matrix)
row.names(dispersal)<-dispersal$Species
dispersal$Species<-NULL
# om_dispersal$number_nas<-NULL

head(dispersal) ## traits per species
head(fish_mat3) ### abundance of species per site year
names(fish_mat3)
dim(fish_mat3)

# # create a "clean" df (called fish for traits "fish_fortr") with fish abundance for the functcomp command (weighting traits by spp relative abund)
fish_fortrt<-fish_mat3[,c(5:202)]
row.names(fish_fortrt)<-fish_mat3$site_year
dim(fish_fortrt)
dim(dispersal)
## computes the functional composition of functional traits, by community weighted mean
disp_matrix<-functcomp(dispersal, as.matrix(fish_fortrt), CWM.type = "dom")  
head(disp_matrix)
disp_matrix <- disp_matrix %>%
  mutate(SiteYear = rownames(disp_matrix)) %>%
  pivot_longer(cols = c(AVG_MXL, HdBd:CFdCPd), names_to="Trait", values_to= "Values") %>%
  mutate(TraitGroup = "Dispersal", TraitType = "Both")

## functcomp requires spp and site names to be row names and column names etc, and the df should not contain other info ##
reproduction<- reproduction[order(reproduction$Species),] # sort species names in the trait df (they should match the fish matrix)
row.names(reproduction)<-reproduction$Species
reproduction$Species<-NULL
reproduction$AVG_RGUILD_ORD<-NULL #### need to fix this
# om_reproduction$number_nas<-NULL

head(reproduction) ## traits per species
str(reproduction)
head(fish_mat3) ### abundance of species per site year
names(fish_mat3)
dim(fish_mat3)
str(fish_mat3)

# # create a "clean" df (called fish for traits "fish_fortr") with fish abundance for the functcomp command (weighting traits by spp relative abund)
fish_fortrt<-fish_mat3[,c(5:202)]
row.names(fish_fortrt)<-fish_mat3$site_year
dim(fish_fortrt)
dim(reproduction)
## computes the functional composition of functional traits, by community weighted mean
repro_matrix<-functcomp(reproduction, as.matrix(fish_fortrt), CWM.type = "dom")  

repro_matrix <- repro_matrix %>%
  mutate(SiteYear = rownames(repro_matrix)) %>%
  # mutate(AVG_RGUILD_ORD = as.numeric(AVG_RGUILD_ORD)) %>%
  pivot_longer(cols= c(AVG_LMAT:AVG_EGGSIZE), names_to="Trait", values_to= "Values") %>%
  mutate(TraitGroup = "Reproduction", TraitType = "Response")

## functcomp requires spp and site names to be row names and column names etc, and the df should not contain other info ##
envPref<- envPref[order(envPref$Species),] # sort species names in the trait df (they should match the fish matrix)
row.names(envPref)<-envPref$Species
envPref$Species<-NULL
# om_envPref$number_nas<-NULL

head(envPref) ## traits per species
head(fish_mat3) ### abundance of species per site year
names(fish_mat3)
dim(fish_mat3)

# # create a "clean" df (called fish for traits "fish_fortr") with fish abundance for the functcomp command (weighting traits by spp relative abund)
fish_fortrt<-fish_mat3[,c(5:202)]
row.names(fish_fortrt)<-fish_mat3$site_year
dim(fish_fortrt)
dim(envPref)
## computes the functional composition of functional traits, by community weighted mean
env_matrix<-functcomp(envPref, as.matrix(fish_fortrt), CWM.type = "dom")  
env_matrix <- env_matrix %>%
  mutate(SiteYear = rownames(env_matrix)) %>%
  pivot_longer(cols= c(Tp_pref:Q_pref), names_to="Trait", values_to= "Values") %>%
  mutate(TraitGroup = "EnvPreferences", TraitType = "Response")

# combine all
traitmatrix <- rbind(res_matrix, disp_matrix, repro_matrix, env_matrix)

#### format
site_year_basin<-fish_mat3[,c(1,3)]
head(site_year_basin)
fish_ab_t <- fish_ab2 %>%
  select(HydroBasin, SiteID, Country)
fish_ab_t <- na.omit(fish_ab_t)
# add the basin id 
site_year_basin$HydroBasin<- fish_ab_t$HydroBasin[match(site_year_basin$site, fish_ab_t$SiteID)]
# add the origin
site_year_basin$Country<-fish_ab_t$Country[match(site_year_basin$site, fish_ab_t$SiteID)]
# add the year
# site_year_basin<-cbind(site_year_basin, colsplit(site_year_basin$site_year, "_", c("SiteID", "Year")))
names(site_year_basin)
head(site_year_basin)
site_year_basin <- site_year_basin %>%
  mutate(SiteYear = rownames(site_year_basin))

## add sites 

traitmatrix<-merge(traitmatrix, site_year_basin, by="SiteYear")

# traitmatrix <- merge(traitmatrix, FGroup, by = "Species")
head(traitmatrix)
dim(traitmatrix)

## save
write.csv(traitmatrix, "output_data/01_trt_single_traits_interpolated.csv")


