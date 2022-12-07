## matrix regression

library(tidyverse)
library(tidylog)
library(ecodist)

# install.packages("remotes")
# remotes::install_github("reumandc/mms")

# install.packages("PopGenReport")
library(PopGenReport)
?lgrMMRR


# Single Traits -----------------------------------------------------------

## format matrix of synchrony, temp, flow & dummy variable
## regression

load(file = "output_data/sync/03_sync_data_funcgroup_traitgroup_similarity_euclidean_dist_interpolated.RData") ## syncDF
head(syncDF)

SyncBasin <- syncDF %>%
  select(-X) %>%
  mutate(Euclid_Dist_KM = Euclid_Dist_Meters/1000)

head(SyncBasin)

# test <- SyncBasin %>%
#   filter(Pair == "S10203.S10089")


# Temp and Flow -----------------------------------------------------------

TempSync <- read.csv("output_data/sync/02_temperature_between_all_sites_biogeographic_regions.csv")
head(TempSync)

unique(TempSync$env_var)

TempSync <- rename(TempSync, Pair = X, TempCor = Correlation)
TempSync <- TempSync %>%
  filter(env_var == "clim_max_raw") %>%
  select(TempCor, Site_ID1, Site_ID2) %>%
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  select(-Site_ID2, - Site_ID1)

FlowSync <- read.csv("output_data/sync/02_flow_between_all_sites_biogeographic_regions.csv")
head(FlowSync)

FlowSync <- rename(FlowSync, Pair = X, FlowCor = Correlation)
FlowSync <- FlowSync %>%
  filter(env_var == "qmax_raw") %>%
  select(FlowCor, Site_ID1, Site_ID2) %>%
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  select(-Site_ID2, - Site_ID1)

tail(TempSync$Pair)
tail(FlowSync$Pair)


## make sync values wider and format for join
## some site pairs don't match here, flow sync missing sone pairs - fix!!!!

SyncBasinT <- left_join(SyncBasin, TempSync, by = "Pair")
SyncBasinTF <- left_join(SyncBasinT, FlowSync, by = "Pair")


rm(TempSync, syncDF,SyncBasin, SyncBasinT,FlowSync)

head(SyncBasinTF)
names(SyncBasinTF)

# test <- TempSync %>%
#   filter(Pair == "S10203.S10089")


# Regression --------------------------------------------------------------

# matrix regression. so far does not include interactions

## filter to one region - USA

SyncDataUSA <- SyncBasinTF %>%
  filter(Region == "USA") %>%
  select(Site_ID1, Site_ID2, Trait, TraitGroup, Correlation, Connectivity, Euclid_Dist_KM, TempCor, FlowCor) %>%
  # filter(Trait %in% c("AVG_MXL", "AVG_FECUND")) %>%
  distinct()

## define traits
traits <- unique(SyncDataUSA$Trait)
traits

SyncDataUSA <- SyncDataUSA %>%
  pivot_wider(names_from = Trait, values_from = Correlation)

names(SyncDataUSA)
## run all traits in one function

# matrix regression. so far does not include interactions

MRM(dist(SyncDataUSA[,25]) ~ dist(TempCor) * dist(Connectivity),  data=SyncDataUSA, nperm = 1000)

mrmTrait<-lapply(1:length(traits), function(i)
{
  

  MRM(dist(SyncDataUSA[,8+i]) ~ dist(TempCor)+ dist(FlowCor) * dist(Connectivity),  data=SyncDataUSA, nperm = 1000)
  
  
})

save(mrmTrait, file="models/07_MRM_single_traits_USA.RData")
mrmTrait

## filter to one region - AUS
unique(SyncBasinTF$Region)

SyncDataAUS <- SyncBasinTF %>%
  filter(Region == "Oceania") %>%
  select(Site_ID1, Site_ID2, Trait, TraitGroup, Correlation, Connectivity, Euclid_Dist_KM, TempCor, FlowCor) %>%
  # filter(Trait %in% c("AVG_MXL", "AVG_FECUND")) %>%
  distinct()

## define traits
traits <- unique(SyncDataAUS$Trait)
traits

SyncDataAUS <- SyncDataAUS %>%
  pivot_wider(names_from = Trait, values_from = Correlation)

## run all traits in one function

# matrix regression. so far does not include interactions

mrmTrait<-lapply(1:length(traits), function(i)
{
  
  
  MRM(dist(SyncDataAUS[,8+i]) ~ dist(TempCor)+ dist(FlowCor) * dist(Connectivity),  data=SyncDataAUS, nperm = 1000)
  
  
})

save(mrmTrait, file="models/07_MRM_single_traits_Oceania.RData")
mrmTrait

## filter to one region - AUS
unique(SyncBasinTF$Region)

SyncDataEU <- SyncBasinTF %>%
  filter(Region == "Europe") %>%
  select(Site_ID1, Site_ID2, Trait, TraitGroup, Correlation, Connectivity, Euclid_Dist_KM, TempCor, FlowCor) %>%
  # filter(Trait %in% c("AVG_MXL", "AVG_FECUND")) %>%
  distinct()

## define traits
traits <- unique(SyncDataEU$Trait)
traits

SyncDataEU <- SyncDataEU %>%
  pivot_wider(names_from = Trait, values_from = Correlation)

## run all traits in one function

# matrix regression. so far does not include interactions

mrmTrait<-lapply(1:length(traits), function(i)
{
  
  
  MRM(dist(SyncDataEU[,8+i]) ~ dist(TempCor)+ dist(FlowCor) * dist(Connectivity),  data=SyncDataEU, nperm = 1000)
  
  
})

save(mrmTrait, file="models/07_MRM_single_traits_Europe.RData")
mrmTrait

### testing

# USAMod1 <- MRM(dist(SyncDataUSA[,8+i]) ~ dist(TempCor)+ dist(FlowCor) * dist(Connectivity),  data=SyncDataUSA, nperm = 1000)
# USAMod1
# USAMod2 <- MRM(dist(AVG_FECUND) ~ dist(TempCor)+ dist(FlowCor) * dist(Connectivity),  data=SyncDataUSA, nperm = 1000)
# USAMod2
# 
# USAMod3 <- MRM(dist(AVG_MXL) ~ dist(TempCor)+ dist(FlowCor) * dist(Connectivity),  data=SyncDataUSA, nperm = 1000)
# USAMod4
# USAMod2 <- MRM(dist(AVG_FECUND) ~ dist(TempCor)+ dist(FlowCor) * dist(Connectivity),  data=SyncDataUSA, nperm = 1000)
# USAMod2



# Model selection --------------------------------------------------------

### not working!!!!

mxl_dist <- as.matrix(dist(SyncDataUSA$AVG_MXL))
temp_dist <- as.matrix(dist(SyncDataUSA$TempCor))
flow_dist <- as.matrix(dist(SyncDataUSA$FlowCor))
basin_dist <- as.matrix(dist(SyncDataUSA$Connectivity))

class(mxl_dist)
## Function to rank models using a leave-n-out cross validation (LNOCV) procedure for matrix regression models

mmsrank<-function(mats,model.names=NA,n,maxruns,rank.mod=F)
{
  #error checking
  errcheck_mats("mmsrank",mats)
  errcheck_n("mmsrank",dim(mats[[1]])[1],n,maxruns)
  
  #if the user does not provide a list of models names, make one with 
  #all names
  if(length(model.names)==1 && is.na(model.names)==T){
    model.names<-makenames(length(mats))
  } else
  {
    errcheck_pred("mmsrank",model.names,length(mats))
  }
  
  return(mmsrank_int(mats,model.names,n,maxruns,rank.mod))
}

mats <- list(mxl_dist, temp_dist) 
length(mats)

mmsrank(mats, n=500, maxruns = 10)


# alternative packages ----------------------------------------------------

?lgrMMRR

SyncDataUSA <- SyncBasinTF %>%
  filter(Region == "USA") %>%
  select(Site_ID1, Site_ID2, Trait, TraitGroup, Correlation, Connectivity, Euclid_Dist_KM, TempCor, FlowCor) %>%
  # filter(Trait %in% c("AVG_MXL", "AVG_FECUND")) %>%
  distinct()

## define traits
traits <- unique(SyncDataUSA$Trait)
traits

SyncDataUSA <- SyncDataUSA %>%
  pivot_wider(names_from = Trait, values_from = Correlation)

## remove NAs based on sync column for test of function

head(SyncDataUSA)
## subset to test 
SyncDataUSATest
SyncDataUSATest <- SyncDataUSA[1:1000,]

names(SyncDataUSA)
SyncDataUSA
SyncDataUSA[15:20]


syncmat <- na.omit(dist(SyncDataUSA[,25]))
tempmat <- na.omit(dist(SyncDataUSA$TempCor))
distmat <- na.omit(dist(SyncDataUSA$Euclid_Dist_KM))
conmat <- na.omit(dist(SyncDataUSA$Connectivity)
syncmat
MRM(syncmat ~ tempmat + distmat * conmat, nperm = 1000)

# Make a list of the explanatory (X) matrices.
# Names are optional.  Order doesn't matter.
# Can include more than two matrices, if desired.
Xmats <- list(geography=distmat,ecology=tempmat)
Xmats

# Run MMRR function using genMat as the response variable and Xmats as the explanatory variables.
# nperm does not need to be specified, default is nperm=999)
MMRR(syncmat,Xmats,nperm=999)

# These data should generate results of approximately:
# Coefficient of geography = 0.778 (p = 0.001)
# Coefficient of ecology = 0.167 (p = 0.063)
# Model r.squared = 0.727 (p = 0.001)
# Note that significance values may change slightly due to the permutation procedure.


