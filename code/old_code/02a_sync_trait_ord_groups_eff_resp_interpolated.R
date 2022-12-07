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

getwd()

## upload fish abundance and site data
originaldata <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")
head(originaldata)

## upload community traits - all groups

trait_matrix <- read.csv("output_data/01_trait_ordination_interpolated.csv")

## format
all_groups <- trait_matrix %>%
  select(-X) %>%
  rename(SiteID = site)

head(all_groups)

## add bioregion
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
  
  ## filter to one trait group
  tgData <- all_groups %>%
    filter(TraitGroup == tg[t])
  head(tgData)

    ## define basins
    basinsID<-unique(tgData$HydroBasin) 
    basinsID

  
  ### loop over basins
  for (basin in 1:length(basinsID)) {
    
    basindata<-tgData[tgData$HydroBasin==basinsID[basin],]
    # head(basindata)
    basindata <- basindata[order(basindata$SiteID),]
    
    ### loop over axis
    Ntraits<-unique(basindata$Axis)
    
    for (ax in 1: length(Ntraits)) {
      Ntraits[ax]
      trait_data<-basindata[basindata$Axis==unique(basindata$Axis)[ax],]
      # sum(is.na(trait_data))
      # dim(trait_matrix)
      years <- unique(sort(trait_data$year)) ## define years for columns
      
      # make df wide
      trait_data  <- trait_data %>% 
        dplyr::select(-c(SiteYear) ) %>%
        spread(SiteID, Scores) #%>% ## some NAs appear here, don't have all trait scores for all site years
      # head(trait_data)
      
      
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
colnames(synchrony_axis)<-c("basin_ID", "Axis","Country", "TraitGroup", "Correlation","Site_ID1","Site_ID2")
nlevels(factor(synchrony_axis$basin_ID)) # 43
nlevels(factor(synchrony_axis$Axis)) # 7
sum(is.na(synchrony_axis))

###save results
write.csv(synchrony_axis, "output_data/sync/02_funcgroup_traitgroup_withinsite_ordination_interpolated.csv")


# between basin sites - within biogeogrpahic region -----------------------

## define biogeogrpahic region

regionsID<-unique(all_groups$BiogeoRegion) # 3 regions
regionsID
synchrony_axis = NULL
region = 3
f=1
t=1
ax=1

tg <- unique(all_groups$TraitGroup)
t
## loop over trait group
for(t in 1:length(tg)) {
  
  tgData <- all_groups %>%
    filter(TraitGroup == tg[t])
  head(tgData)
  
  
  ### loop over regions
  for (region in 1:length(regionsID)) {
    
    basindata<-tgData[tgData$BiogeoRegion==regionsID[region],]
    head(basindata)
    basindata <- basindata[order(basindata$SiteID),]
    
    ### loop over axis
    Ntraits<-unique(basindata$Axis)
    Ntraits
    
    for (ax in 1: length(Ntraits)) {
      # Ntraits
      trait_data<-basindata[basindata$Axis==unique(basindata$Axis)[ax],]
      # sum(is.na(trait_data))
      years <- unique(sort(trait_data$year)) ## define years for columns
      # years
      # make df wide
      trait_data  <- trait_data %>% 
        dplyr::select(-c(SiteYear, Country, HydroBasin) ) %>%
        spread(SiteID, Scores) #%>% ## some NAs appear here, don't have all trait scores for all site years
      head(trait_data)
      # flip data
      trait_data <- t(trait_data)[-c(1:5),]
      
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
colnames(synchrony_axis)<-c("Region", "Axis","TraitGroup", "Correlation","Site_ID1","Site_ID2")
nlevels(factor(synchrony_axis$Axis)) # 7
unique(synchrony_axis$Axis)
sum(is.na(synchrony_axis))

###save results
write.csv(synchrony_axis, "output_data/sync/02_funcgroup_traitgroup_between_all_sites_ordination_biogeographic_regions_interpolated.csv")


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
    unique(basindata$TraitGroup)
    ## loop over trait group
    for(t in 1:length(tg)) {
      
      tgData <- basindata %>%
        filter(TraitGroup == tg[t])
  
      Ntraits<-unique(tgData$Axis)
      unique(tgData$TraitGroup)
  
  ## loop over axis
    for (ax in 1: length(Ntraits)) {
      # Ntraits
      trait_data<-tgData[tgData$Axis==unique(tgData$Axis)[ax],]
      # sum(is.na(trait_data))
      years <- unique(sort(trait_data$year)) ## define years for columns
      # years
      
      # make df wide
      trait_data  <- trait_data %>% 
        dplyr::select(-c(SiteYear, Country, HydroBasin) ) %>%
        spread(SiteID, Scores) #%>% ## some NAs appear here, don't have all trait scores for all site years
   
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
  colnames(synchrony_axis)<-c("Region", "Axis","TraitGroup", "Correlation","Site_ID1","Site_ID2", "YearRemoved")
  save(synchrony_axis, file=paste0("output_data/sync/02_", paste(regionsID[region]), "_ordination_interpolated_site_sync_one_out.RData", sep=""))
  rm(synchrony_axis)
}



load(file= "output_data/sync/02_Europe_ordination_interpolated_site_sync_one_out.RData")
head(synchrony_axis)
unique(synchrony_axis$TraitGroup)
