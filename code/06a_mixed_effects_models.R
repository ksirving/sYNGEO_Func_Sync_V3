### mixed effects models

# install.packages("haven")
library(haven)

library(tidyverse)
library(tidylog)
library(ecodist)

library(performance)
library(lme4)

library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(sjlabelled)

library(effects)
library(sjstats) #use for r2 functions
library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
library("easystats") ## multicolinearality
library(lmerMultiMember)
library("scales")

## doers this work with membership
# fit.nlmer = nlmer(y ~ power.f(time, k, a, b) ~ k|id, start=list(nlpars=c(k=2, a=1, b=1)), data=D)
# 
# https://rdrr.io/cran/lme4/man/nlmer.html

# library(devtools)
# install_github("jvparidon/lmerMultiMember")
# install.packages("performance")

# install.packages("remotes")
# remotes::install_github("reumandc/mms")
getwd()
# install.packages("PopGenReport")
# library(PopGenReport)

load(file = "output_data/sync/04_temp_pref_env_dist_no_dupl_pairs.RData")
head(allsyncx)

range(allsyncx$diversity2)
range(allsyncx$distance)
## standardize variables

allsyncx <- allsyncx %>%
  rename(overlap = diversity2)
  # group_by(Region) %>%
  # mutate(distance1 = rescale(distance))# %>%
  # mutate(ValuesSt = (Values - mean(na.omit(Values)))/sd(na.omit(Values))) %>% ## makes negative numbers
  # pivot_wider(names_from=Metric, values_from = Values)

## some NAs in sync. find and remove - from sites with 1 spoecies, cannot compute synchrony
ind <- which(is.na(allsyncx$Sync))
allsyncx[ind,]
allsyncx <- allsyncx[-ind,]
sum(is.na(allsyncx))

allsyncx <- na.omit(allsyncx)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V3/Figures/"


# CORRELATIONS ------------------------------------------------------------

allsyncxCor <- allsyncx %>%
  select(distance, overlap, annual_avg, DistKM)

syncCor <- cor(allsyncxCor)


write.csv(syncCor, "output_data/06_cor_dist_overlap_temp.csv")
# Modes: spatial decay ----------------------------------------------------

names(allsyncx)

head(allsyncx)
dim(allsyncx)

hist(allsyncx$diversity)
hist(allsyncx$overlap)
hist(allsyncx$distance)

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 

## distance as log and sqrt

allsyncx <- allsyncx %>%
  mutate(DistKMLog = log(DistKM+1),
         DistKMsqrt = sqrt(DistKM+1))

range(allsyncx$DistKM)

## rescale diversity measures between 0-1
# allsyncx$diversity <- rescale(allsyncx$diversity)
# allsyncx$diversity2 <- rescale(allsyncx$diversity2)
# allsyncx$distance <- rescale(allsyncx$distance)

# ## make longer to get site names for membership model
allsyncx <- allsyncx %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName") 

Wa <- lmerMultiMember::weights_from_vector(allsyncx$Region)
Wj <- Matrix::fac2sparse(allsyncx$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

## mixed model with distance and connectivity
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8784019/#:~:text=As%20linear%20mixed%2Deffects%20models,associated%20with%20a%20random%20effect.
## paper suggesting that fewer than 5 levels in Random effect is ok when only looking at fixed effects, but not good when analysing random effects

mem_mixed0 <- lmerMultiMember::lmer(log(Sync) ~ DistKM  *Connectivity
                                      + (1 | Region ) + ## add predictors here to get random effect per region
                                      + (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

spatPlot2a <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("DistKM","Connectivity"),
                                 axis.title = c("Euclidean Distance KM", "Thermal Trait Synchrony"), pred.type="re", 
                                 ci.lvl=NA, show.data  = T, dot.size = 0.1)

spatPlot2a 



summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) 
check_singularity(mem_mixed0) ## False

# lattice::qqmath(mem_mixed1)
# plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)

### plots 
class(mem_mixed0) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_spatial_decay.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed0, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

# ?set_theme
mem_mixed0
spatPlot2a <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("DistKM","Connectivity"),
                                 axis.title = c("Euclidean Distance KM", "Thermal Trait Synchrony"), pred.type="re", 
                                 ci.lvl=NA, show.data  = T, dot.size = 0.1)

spatPlot2a
#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "spatial_decay_mem_mod_exp.jpg")
ggsave(spatPlot2a, filename=file.name1, dpi=300, height=8, width=10)


### mixed model only - no membership


modelMix<-lmer(formula = Sync ~ log(Euclid_Dist_Meters)  *Connectivity +
                (1|SiteName),
              data    = allsyncx) 

modelMix<-lmer(formula = log(Sync) ~ DistKM  *Connectivity +
                 (1|SiteName),
               data    = allsyncx) 

summary(modelMix)
anova(modelMix, ddf = "Satterthwaite")
r2_nakagawa(modelMix) 
check_singularity(mem_mixed0) ## False

# lattice::qqmath(mem_mixed1)
# plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)

### plots 
class(mem_mixed0) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_spatial_decay.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed0, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size


spatPlot3 <- sjPlot::plot_model(modelMix, type="pred", terms=c("DistKM","Connectivity"),
                                axis.title = c("Euclidean Distance (km)", "Thermal Trait Synchrony"), pred.type="re", 
                                ci.lvl=NA)

spatPlot3
#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "spatial_decay_mixed_mod_exp.jpg")
ggsave(spatPlot3, filename=file.name1, dpi=300, height=8, width=10)


# Membership Model: Region as fixed ---------------------------------------

names(allsyncx)
head(allsyncx)
dim(allsyncx)

hist(allsyncx$diversity)
hist(sqrt(allsyncx$diversity2))
hist(allsyncx$distance)

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 

## recode region
allsyncx$Region <- as.factor(allsyncx$Region)
allsyncx$Region <- recode_factor(allsyncx$Region,  "USA" = "USA", "Oceania" = "Oceania", "Europe" = "Europe") 
unique(allsyncx$Region)
## rescale diversity measures between 0-1
allsyncx$diversity <- rescale(allsyncx$diversity)
allsyncx$diversity2 <- rescale(allsyncx$diversity2)
allsyncx$distance <- rescale(allsyncx$distance)

# ## make longer to get site names for membership model
allsyncx <- allsyncx %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName")

Wa <- lmerMultiMember::weights_from_vector(allsyncx$Region)
Wj <- Matrix::fac2sparse(allsyncx$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

## model with new diversity measure

mem_mixed0 <- lmerMultiMember::lmer(Sync ~  (overlap*annual_avg + distance*annual_avg)*Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

## add in interactions between other variables
#(annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) 
check_singularity(mem_mixed0) ## False

# lattice::qqmath(mem_mixed1)
# plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)

### plots 
class(mem_mixed0) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_overlap_region_removed.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed0, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

## plot
set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

theme_set(theme_sjplot())

quantile(allsyncx$annual_avg, probs = c(0.05,0.5,0.95))
round(quantile(allsyncx$overlap, probs = c(0.05,0.5,0.95)), digits = 3)

tempDiv <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap", "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDiv
tempDist <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance", "annual_avg [0.68, 0.92, 0.998]"),
                              axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                              legend.title = "Environmental synchrony")
tempDist
tempCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg", "Connectivity"),
                               axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                               legend.title = "Connectivity")

divCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap", "Connectivity"),
                              axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                              legend.title = "Connectivity")
divCon

distCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance", "Connectivity"),
                             axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                             legend.title = "Connectivity")
distCon

tempDistCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("distance","Connectivity", "annual_avg [0.68, 0.92, 0.998]" ),
                               axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                               legend.title = "Connectivity")
tempDistCon

tempDivCon <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("overlap","Connectivity", "annual_avg [0.68, 0.92, 0.998]" ),
                                  axis.title = c("Trait Overlap", "Thermal Trait Synchrony"), 
                                  legend.title = "Connectivity")

tempDivCon2 <- sjPlot::plot_model(mem_mixed0, type="pred", terms= c("annual_avg","Connectivity", "overlap [0, 0.035, 0.795]" ),
                                 axis.title = c("Environmental synchrony", "Thermal Trait Synchrony"), 
                                 legend.title = "Connectivity")
tempDivCon2
?plot_model
## significant interactions
CompPlotSig <- plot_grid(list(divCon, tempCon, tempDivCon, tempDivCon2)) ## grid figures

##insignificant interactions
CompPlotInSig <- plot_grid(list(tempDiv, tempDist, distCon, tempDistCon)) ## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "comp_mod_sig_interactions_region_removed.jpg")
ggsave(CompPlotSig, filename=file.name1, dpi=300, height=12, width=18)

file.name1 <- paste0(out.dir, "comp_mod_insig_interactions_region_removed.jpg")
ggsave(CompPlotInSig, filename=file.name1, dpi=300, height=10, width=15)

?plot_model
tempDist <- compPlot1[[2]]
tempCon <- compPlot1[[3]]
tempReg <- compPlot1[[4]]
DivRegion <- compPlot1[[5]]
DistRegion<- compPlot1[[6]]
ConReg <- compPlot1[[7]]
tempDivReg <- compPlot1[[8]]
tempDistReg <- compPlot1[[9]]
tempConReg <- compPlot1[[10]]



## non significant interactions
CompPlotnon <- plot_grid(list(tempCon, tempReg, ConReg, tempConReg))## grid figures

#terms=c("Euclid_Dist_Meters", "Connectivity")
file.name1 <- paste0(out.dir, "comp_mod_interactions_non_significant_region_removed.jpg")
ggsave(CompPlotnon, filename=file.name1, dpi=300, height=8, width=10)



# compPlot2 <- sjPlot::plot_model(mem_mixed0, type="int", terms=c("annual_avg","distance", "Connectivity","Region"),
#                                 pred.type="re", 
#                                 ci.lvl=NA)
# #axis.title = c("Temp sync", "Thermal Trait Synchrony"),
# compPlot2
# #terms=c("Euclid_Dist_Meters", "Connectivity")
# file.name1 <- paste0(out.dir, "distance_composition_model_by_region.jpg")
# ggsave(compPlot2, filename=file.name1, dpi=300, height=8, width=10)


# conPlot3 <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("annual_avg","overlap", "distance","Connectivity"),
#                                 pred.type="re", 
#                                 ci.lvl=NA)
# 
# conPlot3
# #terms=c("Euclid_Dist_Meters", "Connectivity")
# file.name1 <- paste0(out.dir, "composition_model_by_connectivity.jpg")
# ggsave(conPlot3, filename=file.name1, dpi=300, height=8, width=10)
# 
# 
# RegPlot3 <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("annual_avg","overlap", "Region","Connectivity"),
#                                 pred.type="re", 
#                                ci.lvl=NA)
# 
# RegPlot3
# #terms=c("Euclid_Dist_Meters", "Connectivity")
# file.name1 <- paste0(out.dir, "composition_model_by_connectivityRegion.jpg")
# ggsave(RegPlot3, filename=file.name1, dpi=300, height=8, width=10)



# Membership mixed model --------------------------------------------------

names(allsyncx)
head(allsyncx)
dim(allsyncx)

# hist(allsyncx$Sync)
# hist(allsyncx$distance)
hist(log(allsyncx$diversity))
hist(log(allsyncx$diversity2))
# hist(allsyncx$annual_avg)

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 

## rescale diversity measures between 0-1
# allsyncx$diversity <- rescale(allsyncx$diversity)
# allsyncx$diversity2 <- rescale(allsyncx$diversity2)
# allsyncx$distance <- rescale(allsyncx$distance)

# ## make longer to get site names for membership model
allsyncx <- allsyncx %>%
  pivot_longer(Site_ID1:Site_ID2, names_to = "SiteNumber", values_to = "SiteName")

Wa <- lmerMultiMember::weights_from_vector(allsyncx$Region)
Wj <- Matrix::fac2sparse(allsyncx$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)


## original model removing 3-way

mem_mixed <- lmerMultiMember::lmer(Sync ~ (annual_avg*log(diversity2) + distance*log(diversity2)) * Connectivity
                                   + (1 | Region) + ## add predictors here to get random effect per region
                                     (1 | RegionXSiteName), 
                                   memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                   REML = T,
                                   data = allsyncx)

summary(mem_mixed, ddf = "Satterthwaite")
anova(mem_mixed, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed) 
check_singularity(mem_mixed) ## False

# lattice::qqmath(mem_mixed)
# plot(mem_mixed, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed)

### plots 
class(mem_mixed) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_original.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab


## model with new diversity measure

mem_mixed0 <- lmerMultiMember::lmer(Sync ~  annual
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)
# (annual_avg+diversity2 + distance+Connectivity)^2
# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed0, ddf = "Satterthwaite")
anova(mem_mixed0, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed0) 
check_singularity(mem_mixed0) ## False

# lattice::qqmath(mem_mixed1)
# plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)

### plots 
class(mem_mixed0) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed0, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_diversity2.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed0, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

##plot
summary(mem_mixed0)


lmmod <-  lm(Sync ~ 1 + as.factor(Region), data=allsyncx)
summary(lmmod)

tempPlot <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("annual_avg","Region", "Connectivity"),
                               axis.title = c("Temperature Synchrony", "Thermal Trait Synchrony"), pred.type="re", 
                               ci.lvl=NA)




divPlot <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("diversity2","Region", "Connectivity"),
                              axis.title = c("Trait Diversity", "Thermal Trait Synchrony"), 
                              pred.type="re", ci.lvl=NA)
divPlot

distPlot <- sjPlot::plot_model(mem_mixed0, type="pred", terms=c("distance","Region", "Connectivity"),
                               axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                               pred.type="re", ci.lvl=NA)
distPlot




effects_diversity <- effects::effect(term= "diversity2:Connectivity", mod= mem_mixed0)
summary(effects_diversity) #output of what the values are


x_div <- as.data.frame(effects_diversity) 
names(x_div)[3] <- "Sync"
names(x_div)[1] <- "Value"
x_div$Variable <- "diversity2"

effects_annualAvg <- effects::effect(term= "annual_avg:Connectivity", mod= mem_mixed0)
summary(effects_annualAvg) #output of what the values are
effects_annualAvg
x_temp <- as.data.frame(effects_annualAvg)
names(x_temp)[3] <- "Sync"
names(x_temp)[1] <- "Value"
x_temp$Variable <- "annual_avg"
x_temp

effects_distance <- effects::effect(term= "distance:Connectivity", mod= mem_mixed0)
summary(effects_distance) #output of what the values are
effects_distance
x_dist <- as.data.frame(effects_distance)
names(x_dist)[3] <- "Sync"
names(x_dist)[1] <- "Value"
x_dist$Variable <- "distance"
x_dist

x_all <- bind_rows(x_div, x_temp, x_dist)
x_all
## add column for variable name
## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin") 

x_all$Connectivity <- as.factor(x_all$Connectivity)
x_all$Connectivity <- recode_factor(x_all$Connectivity, "1" = "Within Basin", "0" = "Between Basin")

head(allsyncx)

allsyncx1 <- allsyncx %>%
  pivot_longer(c(diversity2, annual_avg, distance), names_to = "Variable", values_to = "Value")

allsyncx1

div1<-ggplot() +
  # geom_point(data = allsyncx1, aes(x=Value, y=Sync, colour = Variable), size = 0.5) +
  geom_line(data=x_all, aes(x=Value, y=Sync, colour = Variable)) +
  geom_ribbon(data= x_all, aes(x=Value, ymin=lower, ymax=upper, colour = Variable), alpha= 0.3) +
  facet_wrap(~Connectivity) +
  # # scale_color_hue(labels = c("Within Basin", "Between Basin")) +
  # scale_color_discrete(labels = c("Within Basin", "Between Basin")) +
  labs(x="Variable", y="Trait Synchrony")

div1

file.name1 <- paste0(out.dir, "random_intercept_model.jpg")
ggsave(div1, filename=file.name1, dpi=300, height=7, width=8)

## combine plots
library("cowplot")
library(ggpubr)
# install.packages("ggpubr")

tempSl <- plot_grid(tempPlot, divPlot, distPlot, div1, ### this not working anymore??
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
tempSl
file.name1 <- paste0(out.dir, "interactions_random_intercept.jpg")
ggsave(tempSl, filename=file.name1, dpi=300, height=8, width=10)


# stepwise model runs -----------------------------------------------------


mem_mixed2 <- lmerMultiMember::lmer(Sync ~  annual_avg
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## positive temp
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  diversity2
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## positive diversity - overlap
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  distance
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## positive distance
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  annual_avg + distance
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## both positive
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) ## lower r2

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  annual_avg + diversity2
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## both negative
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) ## highest r2 - 0.45

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  annual_avg + Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## both positive
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) ## lower r2

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  distance + diversity2 
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## both positive
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  distance + diversity2 +annual_avg
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## neg temp, pos dist & div
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  distance + diversity2 +annual_avg + Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## neg temp & div, pos dist & con
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 


mem_mixed2 <- lmerMultiMember::lmer(Sync ~  (annual_avg*diversity2 + distance*diversity2) * Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## neg dist and div, neg connectivity & temp
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  (annual_avg+diversity2 + distance+Connectivity)^2
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## neg dist, temp and div, pos connectivity
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  (annual_avg+distance+Connectivity)^2
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## neg dist & temp , pos connectivity &dist
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  (annual_avg+distance+diversity2)*Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## neg div, temp , pos connectivity &dist
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  annual_avg+distance+Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## pos dist, temp, neg conn
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 



## (diversity2 +distance+annual_avg)*Connectivity



# random slope tutorial ---------------------------------------------------


## tutorial
popular2data <- read_sav(file ="https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%202/popularity/SPSS/popular2.sav?raw=true")

popular2data <- select(popular2data, pupil, class, extrav, sex, texp, popular)

ggplot(data      = popular2data,
       aes(x     = extrav,
           y     = popular,
           col   = class,
           group = class))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(100))+
  geom_smooth(method = lm,
              se     = FALSE,
              linewidth   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "Popularity vs. Extraversion",
       subtitle = "add colours for different classes and regression lines")

interceptonlymodel <- lmer(formula = popular ~ 1 + (1|class),
                           data    = popular2data) #to run the model

summary(interceptonlymodel) #to get paramater estimates.

library(sjstats)
performance::icc(interceptonlymodel)

model1 <- lmer(formula = popular ~ 1 + sex + extrav + (1|class), 
               data    = popular2data)
summary(model1)

model2 <- lmer(popular ~ 1 + sex + extrav + texp + (1 | class), data=popular2data)
summary(model2)

model3 <- lmer(formula = popular ~ 1 + sex + extrav + texp + (1 + sex + extrav | class),
               data    = popular2data)

summary(model3)
ranova(model3)

model4 <- lmer(formula = popular ~ 1 + sex + extrav + texp + (1 + extrav |class), 
               data    = popular2data)
summary(model4)

model5<-lmer(formula = popular ~ 1 + sex + extrav + texp+ extrav:texp + (1 + extrav | class), 
             data    = popular2data)
summary(model5)

ggplot(data = popular2data,
       aes(x = extrav, 
           y = popular, 
           col = as.factor(texp)))+
  viridis::scale_color_viridis(discrete = TRUE)+
  geom_point(size     = .7,
             alpha    = .8, 
             position = "jitter")+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = 2,
              alpha  = .8)+
  theme_minimal()+
  labs(title    = "Linear Relationship for Different Years of Teacher Experience as Observed", 
       subtitle = "The linear relationship between the two is not the same for all classes", 
       col      = "Years of\nTeacher\nExperience")

# mixed model with random slopes ------------------------------------------
head(allsyncx)
names(allsyncx)
interceptonlymodel <- lmer(formula = Sync ~ 1 + (1|Region),
                           data    = allsyncx) #to run the model

summary(interceptonlymodel) #to get paramater estimates.


performance::icc(interceptonlymodel)

model1 <- lmer(formula = Sync ~ 1 + annual_avg + Connectivity + (1|Region),
               data    = allsyncx) #to run the model
summary(model1)

model2 <- lmer(formula = Sync ~ 1 + annual_avg + Connectivity + distance +diversity2 + (1|Region),
               data    = allsyncx) #to run the model
summary(model2)


model3 <- lmer(formula = Sync ~ 1 + annual_avg + Connectivity + distance +diversity2 + 
                 (1+annual_avg + Connectivity + distance +diversity2|Region),
                data    = allsyncx) # singular!!! 

summary(model3) ## fixed effects not significant
ranova(model3)

model4 <- lmer(formula = Sync ~ 1 + annual_avg + Connectivity + distance +diversity2 + 
                 (1+annual_avg + Connectivity |Region),
               data    = allsyncx) # singular!!!
summary(model4)
ranova(model4)



model5<-lmer(formula = Sync ~ 1 + annual_avg + Connectivity + distance +diversity2 + 
               (1+annual_avg + diversity2 |Region),
             data    = allsyncx) # singular!!!
summary(model5)
ranova(model5)


model6<-lmer(formula = Sync ~ 1 + annual_avg + Connectivity + distance +diversity2 + 
               (1+annual_avg |Region),
             data    = allsyncx) # not singular
summary(model6)
ranova(model6)

model7<-lmer(formula = Sync ~ (distance + annual_avg  + diversity2) *  Connectivity + 
               (1+annual_avg |Region),
             data    = allsyncx) # not singular!!!
summary(model7)
ranova(model7)


# Mixed model per region --------------------------------------------------

## mixed model with sitename as random effect
unique(allsyncx$Region)
## subset data
allsyncxEU <- allsyncx %>%
  filter(Region == "Europe")

allsyncxUS <- allsyncx %>%
  filter(Region == "USA")

allsyncxAU <- allsyncx %>%
  filter(Region == "Oceania")

## models per region

modelEU<-lmer(formula = Sync ~ (distance + annual_avg  + diversity2) *  Connectivity + 
               (1|SiteName),
             data    = allsyncxEU) 
summary(modelEU)
anova(modelEU)
r2_nakagawa(modelEU) ## 0.38

modelUS<-lmer(formula = Sync ~ (distance + annual_avg  + diversity2) *  Connectivity + 
                (1|SiteName),
              data    = allsyncxUS) 
summary(modelUS)
anova(modelUS)
r2_nakagawa(modelUS) ## 0.067

modelAU<-lmer(formula = Sync ~ (distance + annual_avg  + diversity2) *  Connectivity + 
                (1|SiteName),
              data    = allsyncxAU) 
# summary(modelAU)
# anova(modelAU)
r2_nakagawa(modelAU) ## 0.09

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.2,  #To change axis title size
          axis.textsize.x = 1,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

estsEU <- sjPlot::plot_model(modelEU, 
                           show.values=TRUE, show.p=TRUE,
                           title="")

estsEU
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_sep_eu_mod.jpg")
ggsave(estsEU, filename=file.name1, dpi=300, height=8, width=10)

estsUS <- sjPlot::plot_model(modelUS, 
                           show.values=TRUE, show.p=TRUE,
                           title="")

estsUS
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_sep_USA_mod.jpg")
ggsave(estsUS, filename=file.name1, dpi=300, height=8, width=10)

estsAU <- sjPlot::plot_model(modelAU, 
                           show.values=TRUE, show.p=TRUE,
                           title="")

estsAU
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_sep_AU_mod.jpg")
ggsave(estsAU, filename=file.name1, dpi=300, height=8, width=10)

## combine plots
# allplot <- plot_grid(estsEU + theme(axis.text.y = element_blank(),
#                                   axis.ticks.y = element_blank(),
#                                   axis.title.y = element_blank()), 
#                      estsUS + theme(axis.text.y = element_blank(),
#                                     axis.ticks.y = element_blank(),
#                                     axis.title.y = element_blank()),
#                       estsAU + theme(axis.text.y = element_blank(),
#                                      axis.ticks.y = element_blank(),
#                                      axis.title.y = element_blank()),
#                      labels = c("EU", "US", "AU"),
#                      # align = "v",
#                     nrow = 1)
# 
# allplot
# 
# file.name1 <- paste0(out.dir, "sync_sep_models_allregions.jpg")
# ggsave(allplot, filename=file.name1, dpi=300, height=5, width=6)

# estsTab <- sjPlot::tab_model(modelEU, 
#                              show.re.var= TRUE, 
#                              pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
#                              dv.labels= "Drivers of Thermal Synchrony")
# 
# estsTab



## interaction plots

effects_diversity <- effects::effect(term= "diversity2:Connectivity", mod= modelEU)
summary(effects_diversity) #output of what the values are


x_div <- as.data.frame(effects_diversity) 
names(x_div)[3] <- "Sync"
names(x_div)[1] <- "Value"
x_div$Variable <- "diversity2"

effects_annualAvg <- effects::effect(term= "annual_avg:Connectivity", mod= modelEU)
summary(effects_annualAvg) #output of what the values are

x_temp <- as.data.frame(effects_annualAvg)
names(x_temp)[3] <- "Sync"
names(x_temp)[1] <- "Value"
x_temp$Variable <- "annual_avg"
x_temp

effects_distance <- effects::effect(term= "distance:Connectivity", mod= modelEU)
summary(effects_distance) #output of what the values are

x_dist <- as.data.frame(effects_distance)
names(x_dist)[3] <- "Sync"
names(x_dist)[1] <- "Value"
x_dist$Variable <- "distance"
x_dist

x_all <- bind_rows(x_div, x_temp, x_dist)
x_all
## add column for variable name
## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "1" = "Within Basin", "0" = "Between Basin")

x_all$Connectivity <- as.factor(x_all$Connectivity)
x_all$Connectivity <- recode_factor(x_all$Connectivity, "1" = "Within Basin", "0" = "Between Basin")

head(allsyncxEU)

allsyncx1EU <- allsyncxEU %>%
  pivot_longer(c(diversity2, annual_avg, distance), names_to = "Variable", values_to = "Value")



EU1<-ggplot() +
  # geom_point(data = allsyncx1, aes(x=Value, y=Sync, colour = Variable), size = 0.5) +
  geom_line(data=x_all, aes(x=Value, y=Sync, colour = Variable)) +
  geom_ribbon(data= x_all, aes(x=Value, ymin=lower, ymax=upper, colour = Variable), alpha= 0.3) +
  # annotate("text", x = 0.9, y =1.1, label = "r2 =0.378") +
  # geom_text(data = ann_text,label = "Text") +
  facet_wrap(~Connectivity) +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6, 0.8,1)) +
  # # scale_color_hue(labels = c("Within Basin", "Between Basin")) +
  # scale_color_discrete(labels = c("Within Basin", "Between Basin")) +
  labs(x="Variable", y="Trait Synchrony")

EU1

file.name1 <- paste0(out.dir, "sync_sep_models_EU.jpg")
ggsave(EU1, filename=file.name1, dpi=300, height=5, width=6)

#####USA
effects_diversity <- effects::effect(term= "diversity2:Connectivity", mod= modelUS)
summary(effects_diversity) #output of what the values are


x_div <- as.data.frame(effects_diversity) 
names(x_div)[3] <- "Sync"
names(x_div)[1] <- "Value"
x_div$Variable <- "diversity2"

effects_annualAvg <- effects::effect(term= "annual_avg:Connectivity", mod= modelUS)
summary(effects_annualAvg) #output of what the values are

x_temp <- as.data.frame(effects_annualAvg)
names(x_temp)[3] <- "Sync"
names(x_temp)[1] <- "Value"
x_temp$Variable <- "annual_avg"
x_temp

effects_distance <- effects::effect(term= "distance:Connectivity", mod= modelUS)
summary(effects_distance) #output of what the values are

x_dist <- as.data.frame(effects_distance)
names(x_dist)[3] <- "Sync"
names(x_dist)[1] <- "Value"
x_dist$Variable <- "distance"
x_dist

x_all <- bind_rows(x_div, x_temp, x_dist)
x_all
## add column for variable name
## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "1" = "Within Basin", "0" = "Between Basin")

x_all$Connectivity <- as.factor(x_all$Connectivity)
x_all$Connectivity <- recode_factor(x_all$Connectivity, "1" = "Within Basin", "0" = "Between Basin")

head(allsyncxEU)

allsyncx1US <- allsyncxUS %>%
  pivot_longer(c(diversity2, annual_avg, distance), names_to = "Variable", values_to = "Value")


US1<-ggplot() +
  # geom_point(data = allsyncx1, aes(x=Value, y=Sync, colour = Variable), size = 0.5) +
  geom_line(data=x_all, aes(x=Value, y=Sync, colour = Variable)) +
  geom_ribbon(data= x_all, aes(x=Value, ymin=lower, ymax=upper, colour = Variable), alpha= 0.3) +
  # annotate("text", x = 0.9, y =1.1, label = "r2 =0.378") +
  # geom_text(data = ann_text,label = "Text") +
  facet_wrap(~Connectivity) +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6, 0.8,1)) +
  # # scale_color_hue(labels = c("Within Basin", "Between Basin")) +
  # scale_color_discrete(labels = c("Within Basin", "Between Basin")) +
  labs(x="Variable", y="Trait Synchrony")

US1

file.name1 <- paste0(out.dir, "sync_sep_models_US.jpg")
ggsave(US1, filename=file.name1, dpi=300, height=5, width=6)

#### australia

effects_diversity <- effects::effect(term= "diversity2:Connectivity", mod= modelAU)
summary(effects_diversity) #output of what the values are


x_div <- as.data.frame(effects_diversity) 
names(x_div)[3] <- "Sync"
names(x_div)[1] <- "Value"
x_div$Variable <- "diversity2"

effects_annualAvg <- effects::effect(term= "annual_avg:Connectivity", mod= modelAU)
summary(effects_annualAvg) #output of what the values are

x_temp <- as.data.frame(effects_annualAvg)
names(x_temp)[3] <- "Sync"
names(x_temp)[1] <- "Value"
x_temp$Variable <- "annual_avg"
x_temp

effects_distance <- effects::effect(term= "distance:Connectivity", mod= modelAU)
summary(effects_distance) #output of what the values are

x_dist <- as.data.frame(effects_distance)
names(x_dist)[3] <- "Sync"
names(x_dist)[1] <- "Value"
x_dist$Variable <- "distance"
x_dist

x_all <- bind_rows(x_div, x_temp, x_dist)
x_all
## add column for variable name
## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "1" = "Within Basin", "0" = "Between Basin")

x_all$Connectivity <- as.factor(x_all$Connectivity)
x_all$Connectivity <- recode_factor(x_all$Connectivity, "1" = "Within Basin", "0" = "Between Basin")

head(allsyncxEU)

# allsyncx1US <- allsyncxUS %>%
#   pivot_longer(c(diversity2, annual_avg, distance), names_to = "Variable", values_to = "Value")


AU1<-ggplot() +
  # geom_point(data = allsyncx1, aes(x=Value, y=Sync, colour = Variable), size = 0.5) +
  geom_line(data=x_all, aes(x=Value, y=Sync, colour = Variable)) +
  geom_ribbon(data= x_all, aes(x=Value, ymin=lower, ymax=upper, colour = Variable), alpha= 0.3) +
  # annotate("text", x = 0.9, y =1.1, label = "r2 =0.378") +
  # geom_text(data = ann_text,label = "Text") +
  facet_wrap(~Connectivity) +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6, 0.8,1)) +
  # # scale_color_hue(labels = c("Within Basin", "Between Basin")) +
  # scale_color_discrete(labels = c("Within Basin", "Between Basin")) +
  labs(x="Variable", y="Trait Synchrony")

AU1

file.name1 <- paste0(out.dir, "sync_sep_models_AU.jpg")
ggsave(AU1, filename=file.name1, dpi=300, height=5, width=6)


## combine plots

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.2,  #To change axis title size
          axis.textsize.x = 1.1,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1.1)  #To change y axis text size

allplot <- plot_grid(EU1, US1, AU1,
                     labels = c("EU", "US", "AU"),
                     ncol = 2, nrow = 2)

allplot

file.name1 <- paste0(out.dir, "sync_sep_models_allregions.jpg")
ggsave(allplot, filename=file.name1, dpi=300, height=6, width=9)

# Membership model with random slopes -------------------------------------

Wa <- lmerMultiMember::weights_from_vector(allsyncx$Region)
Wj <- Matrix::fac2sparse(allsyncx$SiteName)  # convert single membership vars to an indicator matrix with fac2sparse()
Waj <- interaction_weights(Wa, Wj)

## model with new diversity measure 

## random slope for temp

mem_mixed1 <- lmerMultiMember::lmer(Sync ~  (diversity2 +distance+annual_avg)*Connectivity
                                    + (1 + annual_avg| Region ) + ## add predictors here to get random effect per region
                                      (1 + annual_avg| RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed1, ddf = "Satterthwaite") 
# ranova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed1) ## 0.573

### plots 
class(mem_mixed1) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
estsTemp <- sjPlot::plot_model(mem_mixed1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

estsTemp
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_random_slope_temp.jpg")
ggsave(estsTemp, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed1, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.2,  #To change axis title size
          axis.textsize.x = 1,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 0.9)  #To change y axis text size

# ?set_theme

tempPlot <- sjPlot::plot_model(mem_mixed1, type="pred", terms=c("annual_avg","Region", "Connectivity"),
                   axis.title = c("Temperature Synchrony", "Thermal Trait Synchrony"), pred.type="re", 
                   ci.lvl=NA)

tempPlot


divPlot <- sjPlot::plot_model(mem_mixed1, type="pred", terms=c("diversity2","Region", "Connectivity"),
                               axis.title = c("Trait Diversity", "Thermal Trait Synchrony"), 
                              pred.type="re", ci.lvl=NA)
divPlot

distPlot <- sjPlot::plot_model(mem_mixed1, type="pred", terms=c("distance","Region", "Connectivity"),
                              axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                              pred.type="re", ci.lvl=NA)
distPlot

## combine plots
library("cowplot")
library(ggpubr)
# install.packages("ggpubr")

tempSl <- plot_grid(tempPlot, divPlot, distPlot,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
tempSl
file.name1 <- paste0(out.dir, "interactions_random_slope_temp.jpg")
ggsave(tempSl, filename=file.name1, dpi=300, height=8, width=10)

### random slope for diversity


mem_mixed2 <- lmerMultiMember::lmer(Sync ~  (diversity2 +distance+annual_avg)*Connectivity
                                    + (1 + diversity2| Region ) + ## add predictors here to get random effect per region
                                      (1 + diversity2| RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)
# mem_mixed2
summary(mem_mixed2, ddf = "Satterthwaite") 
# ranova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) ## 0.565

### plots 
class(mem_mixed2) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
estsDiv <- sjPlot::plot_model(mem_mixed2, 
                           show.values=TRUE, show.p=TRUE)

estsDiv
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_random_slope_diversity.jpg")
ggsave(estsDiv, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed2, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.3,  #To change axis title size
          axis.textsize.x = 1.2,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

# ?set_theme

tempPlot2 <- sjPlot::plot_model(mem_mixed2, type="pred", terms=c("annual_avg","Region", "Connectivity"),
                               axis.title = c("Temperature Synchrony", "Thermal Trait Synchrony"), pred.type="re", 
                               ci.lvl=NA)

divPlot2 <- sjPlot::plot_model(mem_mixed2, type="pred", terms=c("diversity2","Region", "Connectivity"),
                              axis.title = c("Trait Diversity", "Thermal Trait Synchrony"), 
                              pred.type="re", ci.lvl=NA)
divPlot2

distPlot2 <- sjPlot::plot_model(mem_mixed2, type="pred", terms=c("distance","Region", "Connectivity"),
                               axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                               pred.type="re", ci.lvl=NA)
distPlot2

## combine plots
library("cowplot")
library(ggpubr)
install.packages("ggpubr")

divsl <- plot_grid(tempPlot2, divPlot2, distPlot2,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

file.name1 <- paste0(out.dir, "interactions_random_slope_diversity.jpg")
ggsave(divsl, filename=file.name1, dpi=300, height=8, width=10)

### random slope for distance


mem_mixed3 <- lmerMultiMember::lmer(Sync ~  (diversity2 +distance+annual_avg)*Connectivity
                                    + (1 + distance| Region ) + ## add predictors here to get random effect per region
                                      (1 + distance| RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed3, ddf = "Satterthwaite") 
# ranova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed3) ## 0.46

### plots 
class(mem_mixed3) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
estsDist <- sjPlot::plot_model(mem_mixed3, 
                           show.values=TRUE, show.p=TRUE)

estsDist
file.name1 <- paste0(out.dir, "effect_sizes_diversity2_random_slope_distance.jpg")
ggsave(estsDist, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed3, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.2,  #To change axis title size
          axis.textsize.x = 1,    #To change x axis text size
          # axis.angle.x = 60,      #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 0.9)  #To change y axis text size

# ?set_theme

tempPlot3 <- sjPlot::plot_model(mem_mixed3, type="pred", terms=c("annual_avg","Region", "Connectivity"),
                                axis.title = c("Temperature Synchrony", "Thermal Trait Synchrony"), pred.type="re", 
                                ci.lvl=NA)

divPlot3 <- sjPlot::plot_model(mem_mixed3, type="pred", terms=c("diversity2","Region", "Connectivity"),
                               axis.title = c("Trait Diversity", "Thermal Trait Synchrony"), 
                               pred.type="re", ci.lvl=NA)
divPlot3

distPlot3 <- sjPlot::plot_model(mem_mixed3, type="pred", terms=c("distance","Region", "Connectivity"),
                                axis.title = c("Trait Distance", "Thermal Trait Synchrony"), 
                                pred.type="re", ci.lvl=NA)
distPlot3

## combine plots
library("cowplot")
library(ggpubr)
install.packages("ggpubr")

distsl <- plot_grid(tempPlot3, divPlot3, distPlot3,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

file.name1 <- paste0(out.dir, "interactions_random_slope_distance.jpg")
ggsave(distsl, filename=file.name1, dpi=300, height=8, width=10)


allplot <- plot_grid(tempPlot, divPlot2, distPlot3,
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
allplot

file.name1 <- paste0(out.dir, "interactions_random_slope_all_combined.jpg")
ggsave(allplot, filename=file.name1, dpi=300, height=8, width=10)

### ggplot

## get fitted values

allsyncx$fittedVals <- predict(mem_mixed1)

## make variables long
allsyncx1 <- allsyncx %>%
  pivot_longer(c(diversity2, annual_avg, distance), names_to = "Variable", values_to = "Value")

## get model coefs
model_coefs <- coef(mem_mixed1)$Region %>% 
  rename(Intercept = `(Intercept)`,Connectivity = "ConnectivityWithin Basin") %>% 
  rownames_to_column("Region") %>%
  pivot_longer(diversity2:Connectivity, names_to = "Variable", values_to ="Slopes" )

model_coefs

unique(model_coefs$Variable)
unique(allsyncx1$Variable)
  
allsyncx1Coefs <- left_join(allsyncx1, model_coefs, by = c("Region", "Variable"))

names(allsyncx1Coefs)
  
  model_coef_plot <- ggplot(data = allsyncx1Coefs, 
                            mapping = aes(x = Value, 
                                          y = Sync, 
                                          colour = Region, lty = Variable)) +
    # geom_point(na.rm = T, alpha = 0.5) +
    geom_abline(aes(intercept = Intercept, 
                    slope = Slopes,
                    colour = Region, lty = Variable),size = 0.5) +
    facet_wrap(~Connectivity)+
    scale_y_continuous() +
    scale_x_continuous() +
    theme(legend.position = "top")
  
  # see the plot
  model_coef_plot

## plot with fitted values
ggplot(allsyncx1,aes(Value, Sync, group=Region, col=Region, shape=Variable)) + 
  facet_grid(~Connectivity) +
  geom_line(aes(y=fittedVals), size=0.8) +
  # geom_point(alpha = 0.3) + 
  # geom_hline(yintercept=0, linetype="dashed") +
  theme_bw()



effects_diversity <- effects::effect(term= "diversity2:Connectivity", mod= mem_mixed1)
summary(effects_diversity) #output of what the values are


x_div <- as.data.frame(effects_diversity) 
names(x_div)[3] <- "Sync"
names(x_div)[1] <- "Value"
x_div$Variable <- "diversity2"

effects_annualAvg <- effects::effect(term= "annual_avg:Connectivity", mod= mem_mixed1)
summary(effects_annualAvg) #output of what the values are

x_temp <- as.data.frame(effects_annualAvg)
names(x_temp)[3] <- "Sync"
names(x_temp)[1] <- "Value"
x_temp$Variable <- "annual_avg"
x_temp

effects_distance <- effects::effect(term= "distance:Connectivity", mod= mem_mixed1)
summary(effects_distance) #output of what the values are

x_dist <- as.data.frame(effects_distance)
names(x_dist)[3] <- "Sync"
names(x_dist)[1] <- "Value"
x_dist$Variable <- "distance"
x_dist

x_all <- bind_rows(x_div, x_temp, x_dist)
x_all
## add column for variable name
## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity,  "1" = "Within Basin", "0" = "Between Basin")

x_all$Connectivity <- as.factor(x_all$Connectivity)
x_all$Connectivity <- recode_factor(x_all$Connectivity,  "1" = "Within Basin", "0" = "Between Basin")

head(allsyncx)

allsyncx1 <- allsyncx %>%
  pivot_longer(c(diversity2, annual_avg, distance), names_to = "Variable", values_to = "Value")

allsyncx1

div1<-ggplot() +
  # geom_point(data = allsyncx1, aes(x=Value, y=Sync, colour = Variable), size = 0.5) +
  geom_line(data=x_all, aes(x=Value, y=Sync, colour = Variable)) +
  geom_ribbon(data= x_all, aes(x=Value, ymin=lower, ymax=upper, colour = Variable), alpha= 0.3) +
  facet_wrap(~Connectivity) +
  # # scale_color_hue(labels = c("Within Basin", "Between Basin")) +
  # scale_color_discrete(labels = c("Within Basin", "Between Basin")) +
  labs(x="Variable", y="Trait Synchrony")

div1

file.name1 <- paste0(out.dir, "model_figure.jpg")
ggsave(div1, filename=file.name1, dpi=300, height=5, width=6)






mem_mixed2 <- lmerMultiMember::lmer(Sync ~  (diversity2 +distance+annual_avg)*Connectivity
                                    + (1 + annual_avg + diversity2 +distance| Region ) + ## add predictors here to get random effect per region
                                      (1 + annual_avg+ diversity2 +distance| RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed2, ddf = "Satterthwaite")
ranova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) ## singular! some variance is zero, conditional r2 is NA
check_singularity(mem_mixed2) ## true

mem_mixed3 <- lmerMultiMember::lmer(Sync ~  (diversity2 +distance+annual_avg)*Connectivity
                                    + (1 + diversity2 | Region ) + ## add predictors here to get random effect per region
                                      (1 + diversity2 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed3, ddf = "Satterthwaite")
ranova(mem_mixed3, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed3) ## 0.541
check_singularity(mem_mixed3) ## False

mem_mixed4 <- lmerMultiMember::lmer(Sync ~  (diversity2 +distance+annual_avg)*Connectivity
                                    + (1 + diversity2 + annual_avg| Region ) + ## add predictors here to get random effect per region
                                      (1 + diversity2 + annual_avg| RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed4, ddf = "Satterthwaite")
ranova(mem_mixed4, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed4) ## NA
check_singularity(mem_mixed4) ## True

lattice::qqmath(mem_mixed1)
plot(mem_mixed1, type=c("p","smooth"), col.line=1)
# check_model(mem_mixed1)

### plots 
class(mem_mixed1) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_diversity2.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed1, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab



mem_mixed_log_sync <- lmerMultiMember::lmer(log(Sync) ~ (distance + annual_avg  + log(diversity)) * Connectivity
                                            + (1 | Region) + ## add predictors here to get random effect per region
                                              (1 | RegionXSiteName), 
                                            memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                            REML = T,
                                            data = allsyncx)


# save(mem_mixed, file = "output_data/models/06_multimembership_model.RData")
class(mem_mixed_log_sync)
summary(mem_mixed_log_sync, ddf = "Satterthwaite")
anova(mem_mixed_log_sync, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed_log_sync) 
check_singularity(mem_mixed_log_sync) ## False

## diagnistics
lattice::qqmath(mem_mixed_log_sync)
plot(mem_mixed_log_sync, type=c("p","smooth"), col.line=1)

### plots 
class(mem_mixed_log_sync) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed_no_con)
ests <- sjPlot::plot_model(mem_mixed_log_sync, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_MS_logged_log_sync.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=5, width=6)

estsTab <- sjPlot::tab_model(mem_mixed_log_sync, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab


## model for correlation
mod_cor = lmer(Sync ~ (distance + annual_avg  + log(diversity)) + Connectivity
               + (1 | Region) + ## add predictors here to get random effect per region
                 (1 | RegionXSiteName), 
               memberships = list(Region = Wa, RegionXSiteName = Waj), 
               REML = T,
               data = allsyncx)

check_collinearity(mod_cor) 

syncCor <- cor(data.frame(allsyncx$distance, allsyncx$diversity, allsyncx$annual_avg))
syncCor


## lm to check direction of annual temp

lm1 <- lm(Sync~ (distance + annual_avg  + diversity)*Connectivity, data = allsyncx)
summary(lm1)

lm2 <- lm(Sync~  annual_avg , data = allsyncx)
summary(lm2)

