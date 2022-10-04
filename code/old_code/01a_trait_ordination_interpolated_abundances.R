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


setwd("/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync")
getwd()

## directory for figures
out.dir <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync/"

# Traits ------------------------------------------------------------------


## upload traits
trt1 <- read.csv("input_data/Bio/matSpecies_imputedCloseTaxa.csv")

head(trt1)
dim(trt1) ## 272, 21 traits

## check missing traits

mis <- missmap(trt1) 


## remove CTMax and species for morpho also
## count NAs in each row - nas are missing traits
trt <- trt1 %>%
  select(-CTmax, -Species_used_for_morpho, -imput) %>%
  mutate(number_nas = rowSums(is.na(trt1)))

rare_species <- trt %>%
  filter(number_nas >= 2) ## number of missing traits

RSp <- rare_species$Species
RSp


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
  select(-AVG_Troph, -AVG_RGUILD_ORD)

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

# Interpolate abundances --------------------------------------------------

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
dim(trt)
str(trt)

## match species 
# check the matchin of spp
setdiff(trt$Species, fish_ab_rel_int$Species)
setdiff(fish_ab_rel_int$Species, trt$Species) # all matched!!!

# make sure basin is a factor
fish_ab_rel_int$HydroBasin<-as.factor(fish_ab_rel_int$HydroBasin)


# Trait ordination --------------------------------------------------------


resource <- trt %>%
  select(Species, EdHd:BlBd)

write.csv(resource, "output_data/01_resource_traits_species_for_dist_mat.csv")

dispersal <- trt %>%
  select(Species, AVG_MXL, HdBd:CFdCPd)

write.csv(dispersal, "output_data/01_dispersal_traits_species_for_dist_mat.csv")

reproduction <- trt %>%
  select(Species, AVG_LMAT:AVG_EGGSIZE, AVG_RGUILD_ORD)

write.csv(reproduction, "output_data/01_reproduction_traits_species_for_dist_mat.csv")

envPref <- trt %>%
  select(Species, Tp_pref:Q_pref)

write.csv(envPref, "output_data/01_envPreferences_traits_species_for_dist_mat.csv")

## resource
rownames(resource) <- resource$Species
resource <- resource[,-1] ## remove species column

sum(is.na(resource))
resource
### gower distance matrix
distance.matrix <- vegdist(resource, method="gower", na.rm = T)

# distance.matrix = gower.dist(resource) ## another option

## PCoA
mds.stuff <- wcmdscale(distance.matrix, eig=TRUE, add = "cailliez") ## negative eigenvalues correction

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.var.per

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.values
mds.data <- data.frame(Species=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data

r1 <- ggplot(data=mds.data, aes(x=X, y=Y, label=Species)) +
  geom_point() +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("Food Aquisition Traits: PCoA with Gower")

r1


file.name1 <- paste0(out.dir, "Figures/Resource_trait_ordination.jpg")
ggsave(r1, filename=file.name1, dpi=300, height=5, width=6)

res_ord <- mds.data

## add arrows to plot

## test plotting with plot function
plot(mds.stuff, resource, na.rm = T, add = T)
plot(envfit(mds.stuff, resource, na.rm = T, add = T))

scrs <- res_ord
scrs

## use envfit to get trait loadings
set.seed(123)
vf <- envfit(mds.stuff, resource, perm = 999, na.rm = T)
vf

## get vectors of traits, add trait names
spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Traits = rownames(spp.scrs))

## add easy names
spp.scrs <- spp.scrs %>%
  mutate(TraitName = case_when(Traits == "EdHd" ~ "Eye Size",
                               Traits == "MoBd" ~ "Oral Gape Position",
                               Traits == "JlHd" ~ "Maxillary Length",
                               Traits == "EhBd" ~ "Eye Position",
                               Traits == "BlBd" ~ "Body Elongation"))

spp.scrs
## plot with arrows
p <- ggplot(scrs) +
  geom_point(mapping = aes(x = X, y = Y)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "blue") +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("Food Aquisition Traits: PCoA with Gower") +
  geom_text(data = spp.scrs, aes(x = Dim1, y = Dim2, label = TraitName),
            size = 3, colour = "blue")
p
file.name1
file.name1 <- paste0(out.dir, "Figures/Resource_trait_ordination_arrows.jpg")
ggsave(p, filename=file.name1, dpi=300, height=5, width=6)

## weak relationship = short arrows
?envfit

### reproduction

rownames(reproduction) <- reproduction$Species
reproduction <- reproduction[,-c(1)] ## remove species and X column

reproduction$AVG_RGUILD_ORD <- as.numeric(reproduction$AVG_RGUILD_ORD) ## need to change this

### gower distance matrix
distance.matrix <- vegdist(reproduction, method="gower", na.rm = T) 


# distance.matrix = gower.dist(resource) ## another option

## PCoA
mds.stuff <- wcmdscale(distance.matrix, eig=TRUE, add = "cailliez") ## negative eigenvalues correction

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.var.per

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.values
mds.data <- data.frame(Species=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data

r1 <- ggplot(data=mds.data, aes(x=X, y=Y, label=Species)) +
  geom_point() +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("Reproduction Traits: PCoA with Gower")

r1


file.name1 <- paste0(out.dir, "Figures/Reproduction_trait_ordination.jpg")
ggsave(r1, filename=file.name1, dpi=300, height=5, width=6)

repro_ord <- mds.data

## add arrows to plot

## test plotting with plot function
plot(mds.stuff, reproduction, na.rm = T, add = T)
plot(envfit(mds.stuff, reproduction, na.rm = T, add = T))

scrs <- repro_ord
scrs

## use envfit to get trait loadings
set.seed(123)
vf <- envfit(mds.stuff, reproduction, perm = 999, na.rm = T)
vf

## get vectors of traits, add trait names
spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Traits = rownames(spp.scrs))
spp.scrs$Traits
# ## add easy names
spp.scrs <- spp.scrs %>%
  mutate(TraitName = case_when(Traits == "AVG_LMAT" ~ "Length at Maturity",
                               Traits == "AVG_AGEMAT" ~ "Age at Maturity",
                               Traits == "AVG_LONGY" ~ "Longevity",
                               Traits == "AVG_FECUND" ~ "Fecund",
                               Traits == "AVG_EGGSIZE" ~ "Eggsize",
                               Traits == "AVG_RGUILD_ORD" ~ "Guild"))

spp.scrs
## plot with arrows
p <- ggplot(scrs) +
  geom_point(mapping = aes(x = X, y = Y)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "blue") +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("Reproduction Traits: PCoA with Gower") +
  geom_text(data = spp.scrs, aes(x = Dim1, y = Dim2, label = TraitName),
            size = 3, colour = "blue")
p

file.name1 <- paste0(out.dir, "Figures/Reproduction_trait_ordination_arrows.jpg")
ggsave(p, filename=file.name1, dpi=300, height=5, width=6)

### dispersal

rownames(dispersal) <- dispersal$Species
dispersal <- dispersal[,-c(1)] ## remove species and X column


### gower distance matrix
distance.matrix <- vegdist(dispersal, method="gower", na.rm = T) 


# distance.matrix = gower.dist(resource) ## another option

## PCoA
mds.stuff <- wcmdscale(distance.matrix, eig=TRUE, add = "cailliez") ## negative eigenvalues correction

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.var.per

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.values
mds.data <- data.frame(Species=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data

r1 <- ggplot(data=mds.data, aes(x=X, y=Y, label=Species)) +
  geom_point() +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("Dispersal Traits: PCoA with Gower")

r1


file.name1 <- paste0(out.dir, "Figures/Dispersal_trait_ordination.jpg")
ggsave(r1, filename=file.name1, dpi=300, height=5, width=6)

disp_ord <- mds.data

## test plotting with plot function
plot(mds.stuff, dispersal, na.rm = T, add = T)
plot(envfit(mds.stuff, dispersal, na.rm = T, add = T))

scrs <- disp_ord
scrs

## use envfit to get trait loadings
set.seed(123)
vf <- envfit(mds.stuff, dispersal, perm = 999, na.rm = T)
vf

## get vectors of traits, add trait names
spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Traits = rownames(spp.scrs))
spp.scrs$Traits
# ## add easy names
spp.scrs <- spp.scrs %>%
  mutate(TraitName = case_when(Traits == "AVG_MXL" ~ "Max. Length",
                               Traits == "HdBd" ~ "Lateral Shape",
                               Traits == "PFiBd" ~ "Pectoral Fin Position",
                               Traits == "PFlBl" ~ "Pectoral Fin Size",
                               Traits == "CFdCPd" ~ "Throttling"))

spp.scrs
## plot with arrows
p <- ggplot(scrs) +
  geom_point(mapping = aes(x = X, y = Y)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "blue") +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("Dispersal Traits: PCoA with Gower") +
  geom_text(data = spp.scrs, aes(x = Dim1, y = Dim2, label = TraitName),
            size = 3, colour = "blue")
p

file.name1 <- paste0(out.dir, "Figures/Dispersal_trait_ordination_arrows.jpg")
ggsave(p, filename=file.name1, dpi=300, height=5, width=6)


## temp/flow preferences

rownames(envPref) <- envPref$Species
envPref <- envPref[,-c(1)] ## remove species and X column
head(envPref)
### gower distance matrix
distance.matrix <- vegdist(envPref, method="gower", na.rm = T) 


# distance.matrix = gower.dist(resource) ## another option

## PCoA
mds.stuff <- wcmdscale(distance.matrix, eig=TRUE, add = "cailliez") ## negative eigenvalues correction

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.var.per

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.values
mds.data <- data.frame(Species=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data

r1 <- ggplot(data=mds.data, aes(x=X, y=Y, label=Species)) +
  geom_point() +
  theme_bw() +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("Env Preferences Traits: PCoA with Gower")

r1


file.name1 <- paste0(out.dir, "Figures/Env_prefs_trait_ordination.jpg")
ggsave(r1, filename=file.name1, dpi=300, height=5, width=6)

env_ord <- mds.data

## test plotting with plot function
plot(mds.stuff, envPref, na.rm = T, add = T)
plot(envfit(mds.stuff, envPref, na.rm = T, add = T))

scrs <- env_ord
scrs

## use envfit to get trait loadings
set.seed(123)
vf <- envfit(mds.stuff, envPref, perm = 999, na.rm = T)
vf

## get vectors of traits, add trait names
spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Traits = rownames(spp.scrs))
spp.scrs$Traits
# ## add easy names
spp.scrs <- spp.scrs %>%
  mutate(TraitName = case_when(Traits == "Tp_pref" ~ "Temperature",
                               Traits == "Q_pref" ~ "Flow"))

spp.scrs
## plot with arrows
p <- ggplot(scrs) +
  geom_point(mapping = aes(x = X, y = Y)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "blue") +
  xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
  ggtitle("Env Preferences Traits: PCoA with Gower") +
  geom_text(data = spp.scrs, aes(x = Dim1, y = Dim2, label = TraitName),
            size = 3, colour = "blue")
p

file.name1 <- paste0(out.dir, "Figures/Env_prefs_trait_ordination_arrows.jpg")
ggsave(p, filename=file.name1, dpi=300, height=5, width=6)



# Community weighted mean -------------------------------------------------

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


### merge with trait scores - food resources

## functcomp requires spp and site names to be row names and column names etc, and the df should not contain other info ##
res_ord <- res_ord[order(res_ord$Species),] # sort species names in the trait df (they should match the fish matrix)
row.names(res_ord)<-res_ord$Species
res_ord$Species<-NULL
# om_resource$number_nas<-NULL

head(res_ord) ## traits per species
dim(res_ord)
head(fish_mat3) ### abundance of species per site year
names(fish_mat3)
dim(fish_mat3)

# # create a "clean" df (called fish for traits "fish_fortr") with fish abundance for the functcomp command (weighting traits by spp relative abund)
fish_fortrt<-fish_mat3[,c(5:202)]
row.names(fish_fortrt)<-fish_mat3$site_year
dim(fish_fortrt)
dim(resource)
## computes the functional composition of functional traits, by community weighted mean
res_matrix<-functcomp(res_ord, as.matrix(fish_fortrt), CWM.type = "dom")  
head(res_matrix)

res_matrix <- res_matrix %>%
  mutate(SiteYear = rownames(res_matrix)) %>%
  pivot_longer(X:Y, names_to="Axis", values_to= "Scores") %>%
  mutate(TraitGroup = "FoodAquisition", TraitType = "Effect") 


## dispersal
head(disp_ord)
## functcomp requires spp and site names to be row names and column names etc, and the df should not contain other info ##
disp_ord <- disp_ord[order(disp_ord$Species),] # sort species names in the trait df (they should match the fish matrix)
row.names(disp_ord)<-disp_ord$Species
disp_ord$Species<-NULL

head(fish_mat3) ### abundance of species per site year
names(fish_mat3)
dim(fish_mat3)

# # create a "clean" df (called fish for traits "fish_fortr") with fish abundance for the functcomp command (weighting traits by spp relative abund)
fish_fortrt<-fish_mat3[,c(5:202)]
row.names(fish_fortrt)<-fish_mat3$site_year
dim(fish_fortrt)
dim(resource)
## computes the functional composition of functional traits, by community weighted mean
disp_matrix<-functcomp(disp_ord, as.matrix(fish_fortrt), CWM.type = "dom")  
head(disp_matrix)

disp_matrix <- disp_matrix %>%
  mutate(SiteYear = rownames(disp_matrix)) %>%
  pivot_longer(X:Y, names_to="Axis", values_to= "Scores") %>%
  mutate(TraitGroup = "Dispersal", TraitType = "Both") 


## functcomp requires spp and site names to be row names and column names etc, and the df should not contain other info ##
repro_ord<- repro_ord[order(repro_ord$Species),] # sort species names in the trait df (they should match the fish matrix)
row.names(repro_ord)<-repro_ord$Species
repro_ord$Species<-NULL
repro_ord$AVG_RGUILD_ORD<-NULL #### need to fix this
# om_reproduction$number_nas<-NULL

head(repro_ord) ## traits per species
str(repro_ord)
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
repro_matrix<-functcomp(repro_ord, as.matrix(fish_fortrt), CWM.type = "dom")  

repro_matrix <- repro_matrix %>%
  mutate(SiteYear = rownames(repro_matrix)) %>%
  # mutate(AVG_RGUILD_ORD = as.numeric(AVG_RGUILD_ORD)) %>%
  pivot_longer(X:Y, names_to="Axis", values_to= "Scores") %>%
  mutate(TraitGroup = "Reproduction", TraitType = "Response")

## functcomp requires spp and site names to be row names and column names etc, and the df should not contain other info ##
env_ord<- env_ord[order(env_ord$Species),] # sort species names in the trait df (they should match the fish matrix)
row.names(env_ord)<-env_ord$Species
env_ord$Species<-NULL
# om_envPref$number_nas<-NULL

head(env_ord) ## traits per species
head(fish_mat3) ### abundance of species per site year
names(fish_mat3)
dim(fish_mat3)

# # create a "clean" df (called fish for traits "fish_fortr") with fish abundance for the functcomp command (weighting traits by spp relative abund)
fish_fortrt<-fish_mat3[,c(5:202)]
row.names(fish_fortrt)<-fish_mat3$site_year
dim(fish_fortrt)
dim(env_ord)
## computes the functional composition of functional traits, by community weighted mean
env_matrix<-functcomp(env_ord, as.matrix(fish_fortrt), CWM.type = "dom")  

env_matrix <- env_matrix %>%
  mutate(SiteYear = rownames(env_matrix)) %>%
  pivot_longer(X:Y, names_to="Axis", values_to= "Scores") %>%
  mutate(TraitGroup = "EnvPreferences", TraitType = "Response")

# combine all
traitmatrix <- rbind(res_matrix, disp_matrix, repro_matrix, env_matrix)

head(traitmatrix)

#### format for synchrony
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
write.csv(traitmatrix, "output_data/01_trait_ordination_interpolated.csv")

head(traitmatrix)

sum(is.na(traitmatrix)) ## 80 NAs


