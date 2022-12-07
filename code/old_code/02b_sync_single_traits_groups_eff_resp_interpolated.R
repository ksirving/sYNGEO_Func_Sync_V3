### synchrony on single traits

## packages
library(tidyverse)
library(reshape2)
library(tidyr)
library(dplyr)
library(synchrony)
library(codyn)
# library(rfishbase)
library(munfold)
library(data.table)
library(gdata)
library(here)

getwd()
## upload fish abundance and site data
originaldata <- read.csv(here("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv"))
head(originaldata)

# test <- originaldata %>% filter(SiteID == "S6654")
# test

## upload and format community weighted mean traits - all groups

trait_matrix <- read.csv(here("output_data/01_trt_single_traits_interpolated.csv"))

## combine all groups

all_groups <- trait_matrix %>%
  select(-X) %>%
  rename(SiteID = site)

all_groups <- all_groups %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA"))

# within Basin synchrony -----------------------------------------------


synchrony_axis = NULL
basin = 1
f=1
t=1
ax=1
tg
## loop over trait group


tg <- unique(all_groups$TraitGroup)

## loop over trait group
for(t in 1:length(tg)) {
  
  tgData <- all_groups %>%
    filter(TraitGroup == tg[t])
  head(tgData)
  
    ## define basins
    basinsID<-unique(tgData$HydroBasin) # 40 basins

  
  ### loop over basins
  for (basin in 1:length(basinsID)) {
    basindata<-tgData[tgData$HydroBasin==basinsID[basin],]
    # head(basindata)
    basindata <- basindata[order(basindata$SiteID),]
    
    ### loop over axis
    Ntraits<-unique(basindata$Trait)
    
    for (ax in 1: length(Ntraits)) {
      Ntraits[ax]
      trait_data<-basindata[basindata$Trait==unique(basindata$Trait)[ax],]
      # sum(is.na(trait_data))
      dim(trait_matrix)
      years <- unique(sort(trait_data$year)) ## define years for columns
      years
      # make df wide
      trait_data  <- trait_data %>% 
        dplyr::select(-c(SiteYear) ) %>%
        spread(SiteID, Values) #%>% ## some NAs appear here, don't have all trait scores for all site years
      head(trait_data)
      # flip data
      trait_data <- t(trait_data)[-c(1:7),]
      
      # define rows and column names and convert to numbers
      sitesIDs<-rownames(trait_data)
      colnames(trait_data) <- years
      trait_data<-apply(trait_data, 2, as.numeric)
      rownames(trait_data)<-sitesIDs
      
      
      ### synchrony
      correlation<-cor(t(trait_data), use = "pairwise.complete.obs")
      head(correlation)
      # write.csv(correlation, paste("output_data/cor_matrix/02_site_sync_", basinsID[basin], 
      #                              "_", Ntraits[ax], ".csv", sep=""))
      
      vector_data_correl<- unmatrix(correlation,byrow=F)
      lower_triangle<-lower.tri(correlation)
      vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
      correl_result<-vector_data_correl[vector_data_triangle]
      site_ID1<-sapply(strsplit(names(correl_result),":"),'[',1)
      site_ID2<-sapply(strsplit(names(correl_result),":"),'[',2)
      # correl_result
      names(basindata)
      synchrony_axis<-rbind(synchrony_axis, cbind(drop.levels(basinsID[basin]), 
                                                  drop.levels(Ntraits[ax]),
                                                  basindata$Country[1], basindata$TraitGroup[1],
                                                  correl_result,site_ID1,site_ID2))
      
      
    }
    
    
  }
  
  
}


head(synchrony_axis)
warnings()
synchrony_axis<-data.frame(synchrony_axis)
colnames(synchrony_axis)<-c("basin_ID", "Trait","Country", "TraitGroup", "Correlation","Site_ID1","Site_ID2")
nlevels(factor(synchrony_axis$basin_ID)) # 43
nlevels(factor(synchrony_axis$Trait)) # 7
sum(is.na(synchrony_axis))

###save results
write.csv(synchrony_axis, "output_data/sync/02_funcgroup_traitgroup_withinsite_single_traits_interpolated.csv")


# between basin sites - within biogeogrpahic region -----------------------

## define biogeogrpahic region

head(all_groups)

regionsID<-unique(all_groups$BiogeoRegion) # 3 regions
regionsID
synchrony_axis = NULL
region = 3
f=1
t=1
ax=1

tg <- unique(all_groups$TraitGroup)

## loop over trait group
for(t in 1:length(tg)) {
  
  tgData <- all_groups %>%
    filter(TraitGroup == tg[t])
  head(tgData)
  
  
  ### loop over regions
  for (region in 1:length(regionsID)) {
    
    basindata<-tgData[tgData$BiogeoRegion==regionsID[region],]
    # head(basindata)
    basindata <- basindata[order(basindata$SiteID),]
    
    ### loop over axis
    Ntraits<-unique(basindata$Trait)
    Ntraits
    
    for (ax in 1: length(Ntraits)) {
      # Ntraits
      trait_data<-basindata[basindata$Trait==unique(basindata$Trait)[ax],]
      # sum(is.na(trait_data))
      years <- unique(sort(trait_data$year)) ## define years for columns
      # years
      # make df wide
      trait_data  <- trait_data %>% 
        dplyr::select(-c(SiteYear, Country, HydroBasin) ) %>%
        spread(SiteID, Values) #%>% ## some NAs appear here, don't have all trait scores for all site years
      # head(trait_data)
      # flip data
      trait_data <- t(trait_data)[-c(1:5),]
      
      # define rows and column names and convert to numbers
      
      sitesIDs<-rownames(trait_data)
      colnames(trait_data) <- years
      trait_data<-apply(trait_data, 2, as.numeric)
      rownames(trait_data)<-sitesIDs
      
      
      ### synchrony
      correlation<-cor(t(trait_data), use = "pairwise.complete.obs")
      head(correlation)
      # write.csv(correlation, paste("output_data/cor_matrix/02_site_sync_", regionsID[region], 
      #                              "_", Ntraits[ax], ".csv", sep=""))
      
      vector_data_correl<- unmatrix(correlation,byrow=F)
      lower_triangle<-lower.tri(correlation)
      vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
      correl_result<-vector_data_correl[vector_data_triangle]
      site_ID1<-sapply(strsplit(names(correl_result),":"),'[',1)
      site_ID2<-sapply(strsplit(names(correl_result),":"),'[',2)
      
      # correl_result
      names(basindata)
      head(basindata)
      synchrony_axis<-rbind(synchrony_axis, cbind(basindata$BiogeoRegion[1], 
                                                  drop.levels(Ntraits[ax]),
                                                  basindata$TraitGroup[1],
                                                  correl_result,site_ID1,site_ID2))
      
      
    }
    
    
  }
  
  
}


head(synchrony_axis)
warnings()
synchrony_axis<-data.frame(synchrony_axis)
colnames(synchrony_axis)<-c("Region", "Trait","TraitGroup", "Correlation","Site_ID1","Site_ID2")
nlevels(factor(synchrony_axis$Trait)) # 7
sum(is.na(synchrony_axis))

###save results
write.csv(synchrony_axis, "output_data/sync/02_funcgroup_traitgroup_between_all_sites_single_traits_biogeographic_regions_interpolated.csv")

# test <- synchrony_axis %>% filter(Pair == "S7990.S7991")
# test
# 


# Leaving one year out ----------------------------------------------------


## synchrony leaving one year out

### loop over regions

regionsID<-unique(all_groups$BiogeoRegion) # 3 regions
regionsID
head(all_groups)
region <- 2
tg <- unique(all_groups$TraitGroup)
tg
t=1
ax=1

### loop over basins
for (region in 1:length(regionsID)) {
  synchrony_axis = NULL
  
  basindata<-all_groups[all_groups$BiogeoRegion==regionsID[region],]
  # head(basindata)
  basindata <- basindata[order(basindata$SiteID),]
  
  ## loop over trait group
  for(t in 1:length(tg)) {
    
    tgData <- basindata %>%
      filter(TraitGroup == tg[t])
    
    Ntraits<-unique(tgData$Trait)
    
    
    ## loop over axis
    for (ax in 1: length(Ntraits)) {
      # Ntraits
      trait_data<-tgData[tgData$Trait==unique(tgData$Trait)[ax],]
      # sum(is.na(trait_data))
      years <- unique(sort(trait_data$year)) ## define years for columns
      # years
      
      # make df wide
      trait_data  <- trait_data %>% 
        dplyr::select(-c(SiteYear, Country, HydroBasin) ) %>%
        spread(SiteID, Values) #%>% ## some NAs appear here, don't have all trait scores for all site years
      
      # flip data
      trait_data <- t(trait_data)[-c(1:5),]
      
      # define rows and column names and convert to numbers
      sitesIDs<-rownames(trait_data)
      colnames(trait_data) <- years
      trait_data<-apply(trait_data, 2, as.numeric)
      rownames(trait_data)<-sitesIDs
      
      # change NAs to zero - can change later if needed
      trait_data[which(is.na(trait_data))] <- 0
      
      ### synchrony - leaving one out
      ## loop over years
      for (y in 1: length(years)) {
        
        year_removed <- paste(years[y])
        trait_data_reduced <- trait_data[,-y] 
        correlation<-cor(t(trait_data_reduced), use = "pairwise.complete.obs")
        vector_data_correl<- unmatrix(correlation,byrow=F)
        lower_triangle<-lower.tri(correlation)
        vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
        correl_result<-vector_data_correl[vector_data_triangle]
        site_ID1<-sapply(strsplit(names(correl_result),":"),'[',1)
        site_ID2<-sapply(strsplit(names(correl_result),":"),'[',2)
        correl_result
        synchrony_axis<-rbind(synchrony_axis, cbind(basindata$BiogeoRegion[1],  
                                                    drop.levels(Ntraits[ax]),
                                                    tgData$TraitGroup[t],
                                                    correl_result,site_ID1,site_ID2, year_removed))
        synchrony_axis
      }
      
    }
  }
  synchrony_axis<-data.frame(synchrony_axis)
  colnames(synchrony_axis)<-c("Region","Trait","TraitGroup", "Correlation","Site_ID1","Site_ID2", "YearRemoved")
  save(synchrony_axis, file=paste0("output_data/sync/02_", paste(regionsID[region]), "_single_traits_interpolated_site_sync_one_out.RData", sep=""))
  rm(synchrony_axis)
}

# Climate synchrony -------------------------------------------------------
load(file="input_data/Env/clim_data_melt_raw_new_sites.RData")
head(melt_clim_raw)

# S10203.S10089 e.g. missing pair, both sites missing 

sites <- all_groups %>%
  select(SiteID, HydroBasin, Country) %>%
  # rename(SiteID = site) %>%
  distinct()

sites <- sites %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA"))

head(sites)
dim(sites)

clim_sites <- full_join(melt_clim_raw, sites, by="SiteID")

clim_sites <- clim_sites %>%
  filter(!year == 2003)

test <- clim_sites %>%
  filter(SiteID == "S10203")

clim_sites <- na.omit(clim_sites)
dim(clim_sites)
length(unique(clim_sites$SiteID))

head(clim_sites)


regionsID<-unique(clim_sites$BiogeoRegion) # 3 regions
regionsID
synchrony_axis = NULL
region = 3
ax=1

### loop over regions
for (region in 1:length(regionsID)) {
  
  basindata<-clim_sites[clim_sites$BiogeoRegion==regionsID[region],]
  head(basindata)
  basindata <- basindata[order(basindata$SiteID),]
  
  ### loop over axis
  Ntraits<-unique(basindata$env_var)
  
  for (ax in 1: length(Ntraits)) {
    # Ntraits
    trait_data<-basindata[basindata$env_var==unique(basindata$env_var)[ax],]
    # sum(is.na(trait_data))
    years <- unique(sort(trait_data$year)) ## define years for columns
    years
    # make df wide
    trait_data  <- trait_data %>% 
      dplyr::select(-c( Country, HydroBasin) ) %>%
      spread(SiteID, Temp) #%>% ## some NAs appear here, don't have all trait scores for all site years
    head(trait_data)
    # flip data
    trait_data <- t(trait_data)[-c(1:3),]
    
    # define rows and column names and convert to numbers
    sitesIDs<-rownames(trait_data)
    colnames(trait_data) <- years
    trait_data<-apply(trait_data, 2, as.numeric)
    rownames(trait_data)<-sitesIDs
    
    
    ### synchrony
    correlation<-cor(t(trait_data), use = "pairwise.complete.obs")
    head(correlation)
    # write.csv(correlation, paste("output_data/cor_matrix/02_site_sync_", regionsID[region], 
    #                              "_", Ntraits[ax], ".csv", sep=""))
    
    vector_data_correl<- unmatrix(correlation,byrow=F)
    lower_triangle<-lower.tri(correlation)
    vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
    correl_result<-vector_data_correl[vector_data_triangle]
    site_ID1<-sapply(strsplit(names(correl_result),":"),'[',1)
    site_ID2<-sapply(strsplit(names(correl_result),":"),'[',2)
    # correl_result
    # names(basindata)
    synchrony_axis<-rbind(synchrony_axis, cbind(drop.levels(regionsID[region]), 
                                                as.character(unique(basindata$env_var)[ax]),
                                                
                                                correl_result,site_ID1,site_ID2))
    
    
  }
  
  
}



head(synchrony_axis)
warnings()
synchrony_axis<-data.frame(synchrony_axis)
colnames(synchrony_axis)<-c("Region", "env_var", "Correlation","Site_ID1","Site_ID2")
nlevels(factor(synchrony_axis$basin_ID)) # 43
nlevels(factor(synchrony_axis$Trait)) # 7
sum(is.na(synchrony_axis))

###save results
write.csv(synchrony_axis, "output_data/sync/02_temperature_between_all_sites_biogeographic_regions.csv")


# Flow synchrony ----------------------------------------------------------

load(file="input_data/Env/flow_data_melt_raw_new_sites.RData")
head(melt_flow_raw)

sites <- all_groups %>%
  select(SiteID, HydroBasin, Country) %>%
  # rename(SiteID = site) %>%
  distinct()

sites <- sites %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA"))

head(sites)
dim(sites)

flow_sites <- full_join(melt_flow_raw, sites, by="SiteID")

flow_sites <- flow_sites %>%
  filter(!year == 2003)
flow_sites <- na.omit(flow_sites)
dim(flow_sites)
length(unique(flow_sites$SiteID)) ## 3 sites missing???

head(flow_sites)


regionsID<-unique(flow_sites$BiogeoRegion) # 3 regions
regionsID
synchrony_axis = NULL
region = 1


### loop over regions
for (region in 1:length(regionsID)) {
  
  basindata<-flow_sites[flow_sites$BiogeoRegion==regionsID[region],]
  head(basindata)
  basindata <- basindata[order(basindata$SiteID),]
  
  ### loop over axis
  Ntraits<-unique(basindata$env_var)
  
  for (ax in 1: length(Ntraits)) {
    # Ntraits
    trait_data<-basindata[basindata$env_var==unique(basindata$env_var)[ax],]
    # sum(is.na(trait_data))
    years <- unique(sort(trait_data$year)) ## define years for columns
    years
    # make df wide
    trait_data  <- trait_data %>% 
      dplyr::select(-c( Country, HydroBasin) ) %>%
      spread(SiteID, Flow) #%>% ## some NAs appear here, don't have all trait scores for all site years
    head(trait_data)
    # flip data
    trait_data <- t(trait_data)[-c(1:3),]
    
    # define rows and column names and convert to numbers
    sitesIDs<-rownames(trait_data)
    colnames(trait_data) <- years
    trait_data<-apply(trait_data, 2, as.numeric)
    rownames(trait_data)<-sitesIDs
    
    
    ### synchrony
    correlation<-cor(t(trait_data), use = "pairwise.complete.obs")
    # head(correlation)
    # write.csv(correlation, paste("output_data/cor_matrix/02_site_sync_", regionsID[region], 
    #                              "_", Ntraits[ax], ".csv", sep=""))
    
    vector_data_correl<- unmatrix(correlation,byrow=F)
    lower_triangle<-lower.tri(correlation)
    vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
    correl_result<-vector_data_correl[vector_data_triangle]
    site_ID1<-sapply(strsplit(names(correl_result),":"),'[',1)
    site_ID2<-sapply(strsplit(names(correl_result),":"),'[',2)
    # correl_result
    # names(basindata)
    synchrony_axis<-rbind(synchrony_axis, cbind(drop.levels(regionsID[region]), 
                                                as.character(unique(basindata$env_var)[ax]),
                                                
                                                correl_result,site_ID1,site_ID2))
    
    
  }
  
  
}



head(synchrony_axis)
warnings()
synchrony_axis<-data.frame(synchrony_axis)
colnames(synchrony_axis)<-c("Region", "env_var", "Correlation","Site_ID1","Site_ID2")
nlevels(factor(synchrony_axis$basin_ID)) # 43
nlevels(factor(synchrony_axis$Trait)) # 7
sum(is.na(synchrony_axis))

###save results
write.csv(synchrony_axis, "output_data/sync/02_flow_between_all_sites_biogeographic_regions.csv")


# Leave one out - Environment ---------------------------------------------

# Climate synchrony -------------------------------------------------------

load(file="input_data/Env/clim_data_melt_raw_new_sites.RData")
head(melt_clim_raw)

sites <- all_groups %>%
  select(SiteID, HydroBasin, Country) %>%
  # rename(SiteID = site) %>%
  distinct()

sites <- sites %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA"))

head(sites)
dim(sites)



clim_sites <- full_join(melt_clim_raw, sites, by="SiteID")

clim_sites <- clim_sites %>%
  filter(!year == 2003)

## NAs are the basins and orgins etc taken out of species df
clim_sites <- na.omit(clim_sites)

dim(clim_sites)
length(unique(clim_sites$SiteID)) ## 767

head(clim_sites)


regionsID<-unique(clim_sites$BiogeoRegion) # 3 regions
regionsID
synchrony_axis = NULL
region = 1
ax = 1

### loop over basins
for (region in 1:length(regionsID)) {
  synchrony_axis = NULL
  
  basindata<-clim_sites[clim_sites$BiogeoRegion==regionsID[region],]
  # head(basindata)
  basindata <- basindata[order(basindata$SiteID),]
  
    
    Ntraits<-unique(basindata$env_var)
    
    
    ## loop over axis
    for (ax in 1: length(Ntraits)) {
      # Ntraits
      trait_data<-basindata[basindata$env_var==unique(basindata$env_var)[ax],]
      # sum(is.na(trait_data))
      years <- unique(sort(trait_data$year)) ## define years for columns
      # years
      
      # make df wide
      trait_data  <- trait_data %>% 
        dplyr::select(-c( Country, HydroBasin) ) %>%
        spread(SiteID, Temp) #%>% ## some NAs appear here, don't have all trait scores for all site years
      head(trait_data)
      # flip data
      trait_data <- t(trait_data)[-c(1:3),]
      
      # define rows and column names and convert to numbers
      sitesIDs<-rownames(trait_data)
      colnames(trait_data) <- years
      trait_data<-apply(trait_data, 2, as.numeric)
      rownames(trait_data)<-sitesIDs
      
      # change NAs to zero - can change later if needed
      # trait_data[which(is.na(trait_data))] <- 0
      
      ### synchrony - leaving one out
      ## loop over years
      for (y in 1: length(years)) {
        
        year_removed <- paste(years[y])
        trait_data_reduced <- trait_data[,-y] 
        correlation<-cor(t(trait_data_reduced), use = "pairwise.complete.obs")
        vector_data_correl<- unmatrix(correlation,byrow=F)
        lower_triangle<-lower.tri(correlation)
        vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
        correl_result<-vector_data_correl[vector_data_triangle]
        site_ID1<-sapply(strsplit(names(correl_result),":"),'[',1)
        site_ID2<-sapply(strsplit(names(correl_result),":"),'[',2)
        correl_result
        synchrony_axis<-rbind(synchrony_axis, cbind(basindata$BiogeoRegion[1],  
                                                    drop.levels(Ntraits[ax]),
                                                    correl_result,site_ID1,site_ID2, year_removed))
        synchrony_axis
      }
      
    
  }
  synchrony_axis<-data.frame(synchrony_axis)
  colnames(synchrony_axis)<-c("Region","env_var", "Correlation","Site_ID1","Site_ID2", "YearRemoved")
  save(synchrony_axis, file=paste0("output_data/sync/02_", paste(regionsID[region]), "_temperature_between_all_sites_biogeographic_regions_interpolated_sync_one_out.RData", sep=""))
  rm(synchrony_axis)
}


# Flow synchrony ----------------------------------------------------------


load(file="input_data/Env/flow_data_melt_raw_new_sites.RData")
head(melt_flow_raw)

sites <- all_groups %>%
  dplyr::select(SiteID, HydroBasin, Country) %>%
  # rename(SiteID = site) %>%
  distinct()

sites <- sites %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA"))

head(sites)
dim(sites)

flow_sites <- full_join(melt_flow_raw, sites, by="SiteID")

clim_sx <- unique(clim_sites$SiteID)
length(clim_sx)



flow_sites <- flow_sites %>%
  filter(!year == 2003)
flow_sites <- na.omit(flow_sites)
dim(flow_sites)
length(unique(flow_sites$SiteID)) ## 4 sites missing???

head(flow_sites)


## missing sites - "S6801" "S804"  "S1117" "S1209"

?which

regionsID<-unique(flow_sites$BiogeoRegion) # 3 regions
regionsID
# synchrony_axis = NULL
region = 1


for (region in 1:length(regionsID)) {
  synchrony_axis = NULL
  
  basindata<-flow_sites[flow_sites$BiogeoRegion==regionsID[region],]
  # head(basindata)
  basindata <- basindata[order(basindata$SiteID),]
  
  
  Ntraits<-unique(basindata$env_var)
  
  
  ## loop over axis
  for (ax in 1: length(Ntraits)) {
    # Ntraits
    trait_data<-basindata[basindata$env_var==unique(basindata$env_var)[ax],]
    # sum(is.na(trait_data))
    years <- unique(sort(trait_data$year)) ## define years for columns
    # years
    
    # make df wide
    trait_data  <- trait_data %>% 
      dplyr::select(-c( Country, HydroBasin) ) %>%
      spread(SiteID, Flow) #%>% ## some NAs appear here, don't have all trait scores for all site years
    head(trait_data)
    # flip data
    trait_data <- t(trait_data)[-c(1:3),]
    
    # define rows and column names and convert to numbers
    sitesIDs<-rownames(trait_data)
    colnames(trait_data) <- years
    trait_data<-apply(trait_data, 2, as.numeric)
    rownames(trait_data)<-sitesIDs
    
    # change NAs to zero - can change later if needed
    # trait_data[which(is.na(trait_data))] <- 0
    
    ### synchrony - leaving one out
    ## loop over years
    for (y in 1: length(years)) {
      
      year_removed <- paste(years[y])
      trait_data_reduced <- trait_data[,-y] 
      correlation<-cor(t(trait_data_reduced), use = "pairwise.complete.obs")
      vector_data_correl<- unmatrix(correlation,byrow=F)
      lower_triangle<-lower.tri(correlation)
      vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
      correl_result<-vector_data_correl[vector_data_triangle]
      site_ID1<-sapply(strsplit(names(correl_result),":"),'[',1)
      site_ID2<-sapply(strsplit(names(correl_result),":"),'[',2)
      correl_result
      synchrony_axis<-rbind(synchrony_axis, cbind(basindata$BiogeoRegion[1],  
                                                  drop.levels(Ntraits[ax]),
                                                  correl_result,site_ID1,site_ID2, year_removed))
      synchrony_axis
    }
    
    
  }
  synchrony_axis<-data.frame(synchrony_axis)
  colnames(synchrony_axis)<-c("Region","env_var", "Correlation","Site_ID1","Site_ID2", "YearRemoved")
  save(synchrony_axis, file=paste0("output_data/sync/02_", paste(regionsID[region]), "_flow_between_all_sites_biogeographic_regions_interpolated_sync_one_out.RData", sep=""))
  rm(synchrony_axis)
}

load(file = "output_data/sync/02_Europe_flow_between_all_sites_biogeographic_regions_interpolated_sync_one_out.RData")


nlevels(factor(synchrony_axis$basin_ID)) # 43
nlevels(factor(synchrony_axis$env_var)) # 7
sum(is.na(synchrony_axis))



