### final models

## analysis: functional synchrony ~ distance * connectivity + functional diversity + environmental synchrony

# packages

library(tidyverse)
library(tidylog)
library("easystats")
# install.packages("PopGenReport")
library(PopGenReport)

## directory for figures
out.dir <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V2/Figures/"

## function to normalise data
min_max_norm <- function(x) {
  (x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x)))
}

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

# upload data -------------------------------------------------------------

## functional synchrony

load(file = "output_data/sync/03_sync_traits_CWM_CWV_distances.RData")
funcsync <- syncDF %>%
  rename(Sync = synchrony) %>%
  select(-X) %>%
  filter(!Site_ID1 == Site_ID2)

head(funcsync)
unique(funcsync$Trait)
## environment synchrony
load(file="output_data/sync/03_sync_temp_distances.RData")
tempsync <- syncDF %>%
  # mutate(SyncType = "TSync")  %>%
  pivot_wider(names_from = Metric, values_from = synchrony) %>%
  filter( !Site_ID1 == Site_ID2)## remove pairs comrised of same sites
head(tempsync)

## remove for now  - fix later if needed
tempsync <- na.omit(tempsync)

load(file="output_data/sync/03_sync_flow_distances.RData")

flowsync <- syncDF %>%
  # mutate(SyncType = "QSync") %>%
  pivot_wider(names_from = Metric, values_from = synchrony) %>%
  filter( !Site_ID1 == Site_ID2) ## remove pairs comrised of same sites

flowsync <- na.omit(flowsync)

head(flowsync)

## join temp and flow together 

names(flowsync)
names(tempsync)

envsync <- full_join(tempsync, flowsync, by = c("Pair","Site_ID1", "Site_ID2", "Region", "Connectivity", "Euclid_Dist_Meters",
                                                "Similarity", "MeanLat", "MeanLon", "MaxDist"))
## missing rows from missing flow sites (i think)

head(envsync)

## NAs where sites are missing - different for flow and temp
# which(is.na(envsync$Metric.x))
# test <- envsync[2465413:2465421, ]

### join functinal synchrony

allsync <- full_join(funcsync, envsync, by = c("Pair","Site_ID1", "Site_ID2", "Region", "Connectivity", "Euclid_Dist_Meters",
                                               "Similarity", "MeanLat", "MeanLon", "MaxDist")) %>%
  mutate(DistKM = Euclid_Dist_Meters/1000) %>%
  filter(Trait == c( "AVG_MXL","Tp_pref" )) %>%
  select(-qmean_raw, - qmax_raw, -qmin_raw) ## add back if needed

head(allsync)

object.size(allsync)

unique(allsync$Trait)

save(allsync, file = "output_data/sync/04a_all_sync_for_models.RData")

## split into regions
EurSync <- allsync %>% filter(Region == "Europe")
AusSync <- allsync %>% filter(Region == "Oceania")
USASync <- allsync %>% filter(Region == "USA")

## split into traits - max length
sizeSyncEU <- EurSync %>%
  filter(Trait == "AVG_MXL")

# save(sizeSyncEU, file = "output_data/02_europe_body_size_data_for_model_explo.RData")
# load( file = "output_data/02_europe_body_size_data_for_model_explo.RData") ## sizeSyncEU

sizeSyncAU <- AusSync %>%
  filter(Trait == "AVG_MXL")

sizeSyncUS <- USASync %>%
  filter(Trait == "AVG_MXL")

## split into traits - temp preference
TPSyncEU <- EurSync %>%
  filter(Trait == "Tp_pref")


TPSyncAU <- AusSync %>%
  filter(Trait == "Tp_pref")

TPSyncUS <- USASync %>%
  filter(Trait == "Tp_pref")


# Multicolinearality ------------------------------------------------------

### body sixe
modEur <- lm(Sync~(Connectivity+diversity+DistKM+distance)+annual_avg, data = sizeSyncEU)
summary(modEur)
anova(modEur)
check_collinearity(modEur)
# performance::check_collinearity(modEur)

modUS <- lm(Sync~(Connectivity+diversity+DistKM+distance)+annual_avg, data = sizeSyncUS)
summary(modUS)
anova(modUS)
check_collinearity(modUS) ## high VIF in diversity and annual temp, moderate VIF in DistKM

modAU <- lm(Sync~(Connectivity+diversity+DistKM+distance)+annual_avg, data = sizeSyncAU)
summary(modAU)
anova(modAU)
check_collinearity(modAU)

## temp preference
modEur <- lm(Sync~(Connectivity+diversity+DistKM+distance)+annual_avg, data = TPSyncEU)
summary(modEur)
anova(modEur)
check_collinearity(modEur)

modUS <- lm(Sync~(Connectivity+diversity+DistKM+distance)+annual_avg, data = TPSyncUS)
summary(modUS)
anova(modUS)
check_collinearity(modUS) ## high VIF in diversity and annual temp, moderate VIF in DistKM

modAU <- lm(Sync~(Connectivity+diversity+DistKM+distance)+annual_avg, data = TPSyncAU)
summary(modAU)
anova(modAU)
check_collinearity(modAU)

# full model --------------------------------------------------------------

sizesyncall <- bind_rows(sizeSyncEU, sizeSyncAU, sizeSyncUS)

modAll <- lm(Sync~(Connectivity+diversity+DistKM+distance)+annual_avg, data = sizesyncall)
# summary(modAU)
# anova(modAU)
check_collinearity(modAll)

## Europe
sizeSyncallLog <- sizesyncall %>%
  mutate(fs = log(Sync+1), 
         d = log(DistKM+1),
         dv = diversity,
         di = log(distance+1),
         tav = annual_avg,
         c = as.factor(Connectivity))  %>%
  mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance


## interactions with environmental synchrony
mod2a <- lm(fs~(d+c+dv_Norm+di_Norm)*tav, data = sizeSyncallLog)

summary(mod2a) ## all significant, except functional diversity

## r2 = 0.16

anova(mod2a) ### all significant, single effects and interactions

## interactions with connectivity
mod2a <- lm(fs~(d+tav+dv_Norm+di_Norm)*c, data = sizeSyncallLog)

summary(mod2a) ## all significant, except functional diversity

## full model
mod2a <- lm(fs~(d+dv_Norm+di_Norm)*c*tav, data = sizeSyncallLog)

summary(mod2a) ## all significant, except functional diversity

#produce added variable plots
 avPlots(mod2a, ~ tav+dv_Norm)
 
 avPlots.invis <- function(MODEL, ...) {
   
   ff <- tempfile()
   png(filename = ff)
   OUT <- car::avPlots(MODEL, ...)
   dev.off()
   unlink(ff)
   OUT }
 
 ggAVPLOTS  <- function(MODEL, YLAB = NULL) {
   
   #Extract the information for AV plots
   AVPLOTS <- avPlots.invis(MODEL)
   K       <- length(AVPLOTS)
   
   #Create the added variable plots using ggplot
   GGPLOTS <- vector('list', K)
   for (i in 1:K) {
     DATA         <- data.frame(AVPLOTS[[i]])
     GGPLOTS[[i]] <- ggplot2::ggplot(aes_string(x = colnames(DATA)[1], 
                                                y = colnames(DATA)[2]), 
                                     data = DATA) +
       geom_point(colour = 'blue') + 
       geom_smooth(method = 'lm', se = FALSE, 
                   color = 'red', formula = y ~ x, linetype = 'dashed') +
       xlab(paste0('Predictor Residual \n (', 
                   names(DATA)[1], ' | others)')) +
       ylab(paste0('Response Residual \n (',
                   ifelse(is.null(YLAB), 
                          paste0(names(DATA)[2], ' | others'), YLAB), ')')) }
   
   #Return output object
   GGPLOTS }
 
 ggAVPLOTS(mod2a)
 avPlots.invis(mod2a)
 getwd()
 library(gridExtra)
 PLOTS <- ggAVPLOTS(mod2a)
 K     <- length(PLOTS)
 NCOL  <- ceiling(sqrt(K))
 AVPLOTS <- do.call("arrangeGrob", c(PLOTS, ncol = NCOL, top = 'Added Variable Plots'))
 ggsave('AV Plots - Trucking.jpg', width = 10, height = 10)

file.name1 <- paste0(out.dir, "avplot_body_size_all_regions.jpg")
ggsave(sizeallplot, filename=file.name1, dpi=300, height=5, width=6)

## extract fitted values

sizeSyncallFit <- cbind(mod2a$model, mod2a$fitted.values)
dim(sizeSyncallFit)

sizeSyncallFit <- sizeSyncallFit %>%
  rename(SyncFit = "mod2a$fitted.values") %>%
  mutate(Trait = sizeSyncallLog$Trait[1],
         Region = "All") %>%
  pivot_longer(fs:di_Norm, names_to = "Variable", values_to = "Values")

head(sizeSyncallFit)

### temp pref

TPall <- bind_rows(TPSyncEU, TPSyncAU, TPSyncUS)

modAll <- lm(Sync~(Connectivity+diversity+DistKM+distance)+annual_avg, data = TPall)
# summary(modAU)
# anova(modAU)
check_collinearity(modAll)

## Europe
TPallLog <- TPall %>%
  mutate(fs = log(Sync+1), 
         d = log(DistKM+1),
         dv = diversity,
         di = log(distance+1),
         tav = annual_avg,
         c = Connectivity)  %>%
  mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance


## interactions with environmental synchrony
mod2a <- lm(fs~(d+c+dv_Norm+di_Norm)*tav, data = TPallLog)

summary(mod2a) ## all significant, except functional diversity

## r2 = 0.16

anova(mod2a) ### all significant, single effects and interactions

## interactions with connectivity
mod2a <- lm(fs~(d+tav+dv_Norm+di_Norm)*c, data = TPallLog)

summary(mod2a) ## all significant, except functional diversity

## full model
mod2a <- lm(fs~(d+dv_Norm+di_Norm)*c*tav, data = TPallLog)

summary(mod2a) #

ggplotRegression(mod2a)

library(car)

#produce added variable plots
TPallplot <- avPlots(mod2a, layout= c(5,3))

print(TPallplot)

file.name1 <- paste0(out.dir, "avplot_temp_pref_all_regions.jpg")
ggsave(TPallplot, filename=file.name1, dpi=300, height=5, width=6)


TPSyncallFit <- cbind(mod2a$model, mod2a$fitted.values)
dim(TPSyncallFit)

TPSyncallFit <- TPSyncallFit %>%
  rename(SyncFit = "mod2a$fitted.values") %>%
  mutate(Trait = TPallLog$Trait[1],
         Region = "All") %>%
  pivot_longer(fs:di_Norm, names_to = "Variable", values_to = "Values")

head(TPSyncallFit)

### bind both traits

SyncFitAllRegs <- bind_rows(TPSyncallFit, sizeSyncallFit)

# head(SyncFitAllRegs)
# 
# allsyncRegs <- SyncFitAllRegs %>%
#   pivot_wider(names_from = Variable, values_from = Values) %>%
#   pivot_longer(c(tav, dv_Norm, di_Norm), names_to = "Variable", values_to = "Values")
# 
# head(allsync)
# 
# object.size(SyncFitAll)
# ### save out 
# 
# save(allsync, file = "output_data/sync/04a_fitted_sync_values_for_figures.RData")



# LM: body size----------------------------------------------------------------------


## per region and connectivity?

# sizeSyncEULogW <- sizeSyncEULog %>%
#   filter(Connectivity == 1)
# 
# sizeSyncEULogB <- sizeSyncEULog %>%
#   filter(Connectivity == 0)
# 
# sizeSyncAULogW <- sizeSyncAULog %>%
#   filter(Connectivity == 1)
# 
# sizeSyncAULogB <- sizeSyncAULog %>%
#   filter(Connectivity == 0)
# 
# 
# ## interactions with environmental synchrony
# mod2W <- lm(fs~(d+dv_Norm+di_Norm)*tav, data = sizeSyncEULogW)
# summary(mod2W)
# anova(mod2W)
# 
# mod2B <- lm(fs~(d+dv_Norm+di_Norm)*tav, data = sizeSyncEULogB)
# summary(mod2B)
# anova(mod2B)
# 
# mod2W <- lm(fs~(d+dv_Norm+di_Norm)*tav, data = sizeSyncAULogW)
# summary(mod2W)
# anova(mod2W)
# 
# mod2B <- lm(fs~(d+dv_Norm+di_Norm)*tav, data = sizeSyncAULogB)
# summary(mod2B)
# anova(mod2B)


## Europe
sizeSyncEULog <- sizeSyncEU %>%
  mutate(fs = log(Sync+1), 
         d = log(DistKM+1),
         dv = diversity,
         di = log(distance+1),
         tav = annual_avg,
         c = Connectivity)  %>%
  mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance

head(sizeSyncEULog)
length(unique(sizeSyncEU$Pair))
sum(is.na(unique(sizeSyncEULog$DistKM)))
sum(sizeSyncEULog$DistKM == 0)

## interactions with environmental synchrony
mod2a <- lm(fs~(d+c+dv_Norm+di_Norm)*tav, data = sizeSyncEULog)

summary(mod2a) ## all significant, except functional diversity

## r2 = 0.16

anova(mod2a) ### all significant, single effects and interactions

## interactions with connectivity
mod2a <- lm(fs~(d+tav+dv_Norm+di_Norm)*c, data = sizeSyncEULog)

summary(mod2a) ## all significant, except functional diversity

## full model
mod2a <- lm(fs~(d+dv_Norm+di_Norm)*c*tav, data = sizeSyncEULog)

summary(mod2a) ## all significant, except functional diversity

## r2 = 0.16

anova(mod2a) ### all significant, single effects and interactions

## extract fitted values

sizeSyncEUFit <- cbind(mod2a$model, mod2a$fitted.values)
dim(sizeSyncEUFit)

sizeSyncEUFit <- sizeSyncEUFit %>%
  rename(SyncFit = "mod2a$fitted.values") %>%
  mutate(Trait = sizeSyncEULog$Trait[1],
         Region = sizeSyncEULog$Region[1]) %>%
  pivot_longer(fs:di_Norm, names_to = "Variable", values_to = "Values")

head(sizeSyncEUFit)


## Australia

sizeSyncAULog <- sizeSyncAU %>%
  mutate(fs = log(Sync+1), 
         d = log(DistKM+1),
         dv = diversity,
         di = log(distance+1),
         tav = annual_avg,
         c = Connectivity)  %>%
  mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance

### LM
## interactions with environmental synchrony
mod2b <- lm(fs~(d+c+dv_Norm+di_Norm)*tav, data = sizeSyncAULog)

summary(mod2b) ##
# plot(mod2) ## qq plot heavily tailed

## r2 = 0.03

anova(mod2b) ###

## interactions with connectivity
mod2b <- lm(fs~(d+tav+dv_Norm+di_Norm)*c, data = sizeSyncAULog)

summary(mod2b) ## 
# plot(mod2) ## qq plot heavily tailed

## r2 = 0.05

anova(mod2b) ### 

mod2b <- lm(fs~(d+dv_Norm+di_Norm)*c*tav, data = sizeSyncAULog)

summary(mod2b) ## 
## extract fitted values

sizeSyncAUFit <- cbind(mod2b$model, mod2b$fitted.values)
dim(sizeSyncAUFit)

sizeSyncAUFit <- sizeSyncAUFit %>%
  rename(SyncFit = "mod2b$fitted.values") %>%
  mutate(Trait = sizeSyncAULog$Trait[1],
         Region = sizeSyncAULog$Region[1]) %>%
  pivot_longer(fs:di_Norm, names_to = "Variable", values_to = "Values")

head(sizeSyncAUFit)

## USA

sizeSyncUSLog <- sizeSyncUS %>%
  mutate(fs = log(Sync+1), 
         d = log(DistKM+1),
         dv = diversity,
         di = log(distance+1),
         tav = annual_avg,
         c = Connectivity)  %>%
  mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance

### LM
## interactions with environmental synchrony
mod2c <- lm(fs~(d+c+dv_Norm+di_Norm), data = sizeSyncUSLog) ## remove temp due to colinearity

summary(mod2c) ##  functional diversity, distance, connectivity 
# plot(mod2c) ## qq plot heavily tailed

## r2 = 0.29

anova(mod2c) ### functional diversity, distance, connectivity 

## interactions with connectivity
mod2c <- lm(fs~(d+dv_Norm+di_Norm)*c, data = sizeSyncUSLog) ## remove temp due to colinearity

summary(mod2c) 
# plot(mod2c) ## qq plot heavily tailed

## r2 = 0.29

anova(mod2c) ### functional distance almost siginificant (0.07)

## extract fitted values

sizeSyncUSFit <- cbind(mod2c$model, mod2c$fitted.values)
dim(sizeSyncUSFit)

sizeSyncUSFit <- sizeSyncUSFit %>%
  rename(SyncFit = "mod2c$fitted.values") %>%
  mutate(Trait = sizeSyncUSLog$Trait[1],
         Region = sizeSyncUSLog$Region[1]) %>%
  pivot_longer(fs:di_Norm, names_to = "Variable", values_to = "Values")

head(sizeSyncUSFit)

## bind all fitted values

sizeSyncallFit <- bind_rows(sizeSyncUSFit, sizeSyncEUFit, sizeSyncAUFit)
head(sizeSyncallFit)



# LM: temperature preference ----------------------------------------------

## check fs for skewness

library(e1071)
skewness(na.omit(TPSyncEU$Sync)) ## -0.291 - mod skewed, log transformed not ok
skewness(na.omit(TPSyncAU$Sync)) ## -0.246 - mod skewed, log transformed not ok
skewness(na.omit(TPSyncUS$Sync)) ## -0.198 - mod skewed, log transformed not ok

## Europe
TPSyncEULog <- TPSyncEU %>%
  mutate(fs = log(Sync+1), 
         d = log(DistKM+1),
         dv = diversity,
         di = log(distance+1),
         tav = annual_avg,
         c = Connectivity)  %>%
  mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance

head(TPSyncEULog)
sum(is.na(unique(TPSyncEULog$DistKM)))
sum(TPSyncEULog$DistKM == 0)

## interactions with environmental synchrony
mod3a <- lm(fs~(d+c+dv_Norm+di_Norm)*tav, data = TPSyncEULog)

summary(mod3a) ## all significant, except functional distance, and func dist * temp

## r2 = 0.03

anova(mod3a) ### all significant, except functional distance, and func dist * temp

## interactions with connectivity
mod3a <- lm(fs~(d+tav+dv_Norm+di_Norm)*c, data = TPSyncEULog)

summary(mod3a) ## all significant, except functional distance, and func dist * c

## r2 = 0.03

anova(mod3a) ### all significant, except functional distance, and func dist * c & d*c

## full model

mod3a <- lm(fs~(d+dv_Norm+di_Norm)*c*tav, data = TPSyncEULog)

summary(mod3a) ## 

## extract fitted values

TPSyncEUFit <- cbind(mod3a$model, mod3a$fitted.values)

TPSyncEUFit <- TPSyncEUFit %>%
  rename(SyncFit = "mod3a$fitted.values") %>%
  mutate(Trait = TPSyncEULog$Trait[1],
         Region = TPSyncEULog$Region[1]) %>%
  pivot_longer(fs:di_Norm, names_to = "Variable", values_to = "Values")

head(TPSyncEUFit)

## Australia

TPSyncAULog <- TPSyncAU %>%
  mutate(fs = log(Sync+1), 
         d = log(DistKM+1),
         dv = diversity,
         di = log(distance+1),
         tav = annual_avg,
         c = Connectivity)  %>%
  mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance

### LM
## interactions with environmental synchrony
mod3b <- lm(fs~(d+c+dv_Norm+di_Norm)*tav, data = TPSyncAULog)

summary(mod3b) ## distance km, temp and d*t interaction
# plot(mod2) ## qq plot heavily tailed

## r2 = 0.006

anova(mod3b) ###  temp single, temp*d interaction

## interactions with connectivity
mod3b <- lm(fs~(d+tav+dv_Norm+di_Norm)*c, data = TPSyncAULog)

summary(mod3b) ## distance km, temp and d*t interaction
# plot(mod2) ## qq plot heavily tailed

## r2 = 0.008

anova(mod3b) ###  temp, d*c, diversity*c

## full model

mod3b <- lm(fs~(d+dv_Norm+di_Norm)*c*tav, data = TPSyncAULog)

summary(mod3b)

## extract fitted values

TPSyncAUFit <- cbind(mod3b$model, mod3b$fitted.values)

TPSyncAUFit <- TPSyncAUFit %>%
  rename(SyncFit = "mod3b$fitted.values") %>%
  mutate(Trait = TPSyncAULog$Trait[1],
         Region = TPSyncAULog$Region[1]) %>%
  pivot_longer(fs:di_Norm, names_to = "Variable", values_to = "Values")

## USA

TPSyncUSLog <- TPSyncUS %>%
  mutate(fs = log(Sync+1), 
         d = log(DistKM+1),
         dv = diversity,
         di = log(distance+1),
         tav = annual_avg,
         c = Connectivity)  %>%
  mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance

### LM
## interactions with environmental synchrony
mod3c <- lm(fs~(d+c+dv_Norm+di_Norm), data = TPSyncUSLog) ## remove temp due to colinearity

summary(mod3c) ##  all significant except functional distance
# plot(mod2c) ## qq plot heavily tailed

## r2 = 0.05

anova(mod3c) ### KM, connectivity, func diversity

## interactions with connectivity
mod3c <- lm(fs~(d+dv_Norm+di_Norm)*c, data = TPSyncUSLog) ## remove temp due to colinearity

summary(mod3c) ##  all significant except functional distance
# plot(mod2c) ## qq plot heavily tailed

## r2 = 0.06

anova(mod3c) ### func dist, diversity*c, func dist*c not significant

## extract fitted values

TPSyncUSFit <- cbind(mod3c$model, mod3c$fitted.values)

TPSyncUSFit <- TPSyncUSFit %>%
  rename(SyncFit = "mod3c$fitted.values") %>%
  mutate(Trait = TPSyncUSLog$Trait[1],
         Region = TPSyncUSLog$Region[1]) %>%
  pivot_longer(fs:di_Norm, names_to = "Variable", values_to = "Values")

head(TPSyncUSFit)

## bind all fitted values

TPSyncallFit <- bind_rows(TPSyncUSFit, TPSyncEUFit, TPSyncAUFit)
head(TPSyncallFit)

### bind both traits

SyncFitAll <- bind_rows(TPSyncallFit, sizeSyncallFit, SyncFitAllRegs)

head(SyncFitAll)

allsync <- SyncFitAll %>%
  pivot_wider(names_from = Variable, values_from = Values) %>%
  pivot_longer(c(tav, dv_Norm, di_Norm), names_to = "Variable", values_to = "Values")
  
  head(allsync)

object.size(SyncFitAll)
### save out 

save(allsync, file = "output_data/sync/04a_fitted_sync_values_for_figures.RData")


# Full model --------------------------------------------------------------



# GLM: body size ----------------------------------------------------------

# Europe
moda <- glm(fs~(d+c+di_Norm+dv_Norm)*tav ,data=sizeSyncEULog,family=Gamma(link = "log"))
summary(moda) # -457326 - lowest so far, 
# anova(mod)

with(summary(moda), 1 - deviance/null.deviance) ##  0.1372422

## Australia
modb <- glm(fs~(c+dv_Norm+d+di_Norm)*tav ,data=sizeSyncAULog,family=Gamma(link = "log"))
summary(modb) ##  -6050.9, nothing significant
anova(modb)

with(summary(modb), 1 - deviance/null.deviance) ## 0.02763475

## USA
modc <- glm(fs~(c+dv_Norm+d+di_Norm) ,data=sizeSyncUSLog,family=Gamma(link = "log"))
summary(modc) ##  -6032.9 - connetivity, distance, functional diversity
anova(modc)

with(summary(modc), 1 - deviance/null.deviance) ## 0.2331024



# GLM: Temp pref ----------------------------------------------------------

# Europe
moda <- glm(fs~(d+c+di_Norm+dv_Norm)*tav ,data=TPSyncEULog,family=Gamma(link = "log"))
summary(moda) # 

with(summary(moda), 1 - deviance/null.deviance) ##  0.027

## Australia
modb <- glm(fs~(c+dv_Norm+d+di_Norm)*tav ,data=TPSyncAULog,family=Gamma(link = "log"))
summary(modb) 
anova(modb)

with(summary(modb), 1 - deviance/null.deviance) ## 0.01

## USA
modc <- glm(fs~(c+dv_Norm+d+di_Norm) ,data=TPSyncUSLog,family=Gamma(link = "log"))
summary(modc) ##
anova(modc)

with(summary(modc), 1 - deviance/null.deviance) ## 0.04
