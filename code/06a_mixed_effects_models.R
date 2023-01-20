### mixed effects models


library(tidyverse)
library(tidylog)
library(ecodist)

library(performance)
library(lme4)
library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions
library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
library("easystats") ## multicolinearality
library(lmerMultiMember)
library("scales")

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
## scale variables

allsyncx <- allsyncx %>%
  mutate(distance= rescale(distance, to=c(0,1))) %>%
  mutate(diversity2= rescale(diversity2, to=c(0,1)))

## some NAs in sync. find and remove - from sites with 1 spoecies, cannot compute synchrony
ind <- which(is.na(allsyncx$Sync))
allsyncx[ind,]
allsyncx <- allsyncx[-ind,]
sum(is.na(allsyncx))

allsyncx <- na.omit(allsyncx)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V3/Figures/"

# Membership mixed model --------------------------------------------------

names(allsyncx)
head(allsyncx)
dim(allsyncx)

# hist(allsyncx$Sync)
# hist(allsyncx$distance)
hist(log(allsyncx$diversity))
hist(allsyncx$diversity2)
# hist(allsyncx$annual_avg)

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "0" = "Between Basin", "1" = "Within Basin") 

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

mem_mixed <- lmerMultiMember::lmer(Sync ~ (annual_avg*log(diversity) + distance*log(diversity)) * Connectivity
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

mem_mixed1 <- lmerMultiMember::lmer(Sync ~  (diversity2 +distance+annual_avg)*Connectivity
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

# (annual_avg +distance+diversity2))*Connectivity
summary(mem_mixed1, ddf = "Satterthwaite")
anova(mem_mixed1, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed1) 
check_singularity(mem_mixed1) ## False

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

##plot
summary(mem_mixed1)


lmmod <-  lm(Sync ~ 1 + as.factor(Region), data=allsyncx)
summary(lmmod)


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
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "0" = "Between Basin", "1" = "Within Basin")

x_all$Connectivity <- as.factor(x_all$Connectivity)
x_all$Connectivity <- recode_factor(x_all$Connectivity, "0" = "Between Basin", "1" = "Within Basin")

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



## stepwise model runs

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

summary(mem_mixed2, ddf = "Satterthwaite") ## negative diversity
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

summary(mem_mixed2, ddf = "Satterthwaite") ## neg temp, pos connectivity
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) ## lower r2

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  distance + diversity2 
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## pos dist, neg div
anova(mem_mixed2, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed2) 

mem_mixed2 <- lmerMultiMember::lmer(Sync ~  distance + diversity2 +annual_avg
                                    + (1 | Region ) + ## add predictors here to get random effect per region
                                      (1 | RegionXSiteName), 
                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                    REML = T,
                                    data = allsyncx)

summary(mem_mixed2, ddf = "Satterthwaite") ## neg temp & div, pos dist
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

summary(mem_mixed2, ddf = "Satterthwaite") ## neg dist, temp and div, pos connectivity
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

# mixed model with random slopes ------------------------------------------



mod_mixed = lmer(Sync ~ (distance + annual_avg  + diversity2) *  Connectivity + 
                   (annual_avg | Region) + # random slopes
                   (distance| Region) +
                   (diversity2| Region) +
                   (Connectivity| Region), 
                 REML = T, ## estimates p vals and F test
                 data = allsyncx)
mod_mixed
summary(mod_mixed, ddf = "Satterthwaite")
anova(mod_mixed, ddf = "Satterthwaite")
r2_nakagawa(mod_mixed) 
check_singularity(mod_mixed) ## False

lattice::qqmath(mod_mixed)
plot(mod_mixed, type=c("p","smooth"), col.line=1)
check_model(mod_mixed)

### plots 
# class(mod_mixed) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mod_mixed, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_random_slopes.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mod_mixed, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab



mem_mixed1a <- lmer(Sync ~ (distance + annual_avg  + diversity2) *  Connectivity 
                    + (1 + Region| Region ) + ## add predictors here to get random effect per region
                      # (1 | RegionXSiteName), 
                      # memberships = list(Region = Wa, RegionXSiteName = Waj), 
                      REML = T,
                    data = allsyncx)



mem_mixed <- lmerMultiMember::lmer(Sync ~  ( annual_avg*diversity2  +distance*diversity2)  *Connectivity  #   ++ log(diversity)
                                   + (1 | Region) + ## add predictors here to get random effect per region
                                     (1 | RegionXSiteName), 
                                   memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                   REML = T,
                                   data = allsyncx)

mem_mixed <- lmerMultiMember::lmer(Sync ~  (distance*Connectivity)  #   ++ log(diversity)
                                   + (1 | Region) + ## add predictors here to get random effect per region
                                     (1 | RegionXSiteName), 
                                   memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                   REML = T,
                                   data = allsyncx)


# mem_mixed_glm <- lmerMultiMember::glmer(Sync ~  ( annual_avg*log(diversity)  +distance*log(diversity))  *Connectivity  #   ++ log(diversity)
#                                    + (1 | Region) + ## add predictors here to get random effect per region
#                                      (1 | RegionXSiteName), 
#                                    memberships = list(Region = Wa, RegionXSiteName = Waj), 
#                                    family = gaussian(link = "identity"),
#                                    data = allsyncx)


## model without connectivity

mem_mixed_no_con <- lmerMultiMember::lmer(Sync ~ (distance + annual_avg  + log(diversity)) 
                                          + (1 | Region) + ## add predictors here to get random effect per region
                                            (1 | RegionXSiteName), 
                                          memberships = list(Region = Wa, RegionXSiteName = Waj), 
                                          REML = T,
                                          data = allsyncx)

range(allsyncx$Sync)

# save(mem_mixed, file = "output_data/models/06_multimembership_model.RData")
class(mem_mixed_no_con)
summary(mem_mixed_no_con, ddf = "Satterthwaite")
anova(mem_mixed_no_con, ddf = "Satterthwaite")
r2_nakagawa(mem_mixed_no_con) 
check_singularity(mem_mixed_no_con) ## False

### plots 
class(mem_mixed_no_con) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed_no_con)
ests <- sjPlot::plot_model(mem_mixed_no_con, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_MS_logged_no_con.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=5, width=6)

estsTab <- sjPlot::tab_model(mem_mixed_no_con, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

## model with logged synchrony


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

