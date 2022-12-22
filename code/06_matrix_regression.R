### matrix regression

library(tidyverse)
library(tidylog)
library(ecodist)

# install.packages("remotes")
# remotes::install_github("reumandc/mms")
getwd()
# install.packages("PopGenReport")
# library(PopGenReport)

load(file = "output_data/sync/04_temp_pref_env_dist_no_dupl_pairs.RData")
head(allsyncx)

## some NAs in sync. find and remove - check in sync code later!!!!!
ind <- which(is.na(allsyncx$Sync))
allsyncx[ind,]
allsyncx <- allsyncx[-ind,]
sum(is.na(allsyncx))

allsyncx <- na.omit(allsyncx)

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V3/Figures/"


# Create matrices ---------------------------------------------------------

head(allsyncx)
# ind <- which(is.na(allsyncx$distance))

## temp pref
bio_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "USA") %>%
  dplyr::select(Site_ID1, Site_ID2,  Sync) %>%
  pivot_wider(names_from = Site_ID2, values_from = Sync) %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows

# ## get mirror of lower/upper triangle
bio_wide[upper.tri(bio_wide)] <- t(bio_wide)[upper.tri(bio_wide)]
bio_wide <- as.matrix(bio_wide)
dim(bio_wide)
## make into distance object
bio_wide <- as.dist(bio_wide)
class(bio_wide)
# write.csv(bio_wide, "output_data/temp_pref_sync_matrix_eu.csv")

## temp sync
temp_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "USA") %>%
  dplyr::select(Site_ID1, Site_ID2,  annual_avg) %>%
  pivot_wider(names_from = Site_ID2, values_from = annual_avg)  %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows

## get mirror of lower/upper triangle
temp_wide[upper.tri(temp_wide)] <- t(temp_wide)[upper.tri(temp_wide)]
# # ## convert NAs to 1
temp_wide[is.na(temp_wide)] <- 1
# 
temp_wide <- as.matrix(temp_wide)
dim(temp_wide)
temp_wide <- as.dist(temp_wide)

# as.data.frame(temp_wide[1:5, 1:5]) ##  check visually

MRM(as.dist(bio_wide) ~ as.dist(temp_wide),  nperm=10)

# View(temp_wide)


# write.csv(temp_wide, "output_data/temp_sync_matrix_eu.csv")

## diversity
div_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "USA") %>%
  dplyr::select(Site_ID1, Site_ID2,  diversity) %>%
  pivot_wider(names_from = Site_ID2, values_from = diversity) %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows

div_wide <- as.dist(div_wide)

## get mirror of lower/upper triangle
# div_wide[upper.tri(div_wide)] <- t(div_wide)[upper.tri(div_wide)]
# as.data.frame(div_wide) ##  check visually site-site pairs aren't 1?
div_wide <- as.matrix(div_wide)
# View(div_wide)

# write.csv(div_wide, "output_data/temp_diversity_matrix_eu.csv")

## distance
dist_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region =="USA") %>%
  dplyr::select(Site_ID1, Site_ID2,  distance) %>%
  pivot_wider(names_from = Site_ID2, values_from = distance)   %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows


dist_wide <- as.dist(dist_wide)

## get mirror of lower/upper triangle
# dist_wide[upper.tri(dist_wide)] <- t(dist_wide)[upper.tri(dist_wide)]
# as.data.frame(dist_wide) ##  check visually site-site pairs aren't 1?
dist_wide <- as.matrix(dist_wide)
# View(dist_wide)

# write.csv(dist_wide, "output_data/temp_diversity_matrix_eu.csv")

## connectivity
con_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "USA") %>%
  dplyr::select(Site_ID1, Site_ID2,  Connectivity) %>%
  pivot_wider(names_from = Site_ID2, values_from = Connectivity)   %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows

con_wide <- as.dist(con_wide)
class(con_wide)
## get mirror of lower/upper triangle
# dist_wide[upper.tri(dist_wide)] <- t(dist_wide)[upper.tri(dist_wide)]
# as.data.frame(dist_wide) ##  check visually site-site pairs aren't 1?
con_wide <- as.matrix(con_wide)
# View(dist_wide)

# write.csv(con_wide, "output_data/connectivity_matrix_eu.csv")

## geo distance
km_wide <- allsyncx %>%
  ungroup() %>%
  filter(Region == "USA") %>%
  dplyr::select(Site_ID1, Site_ID2,  DistKM) %>%
  pivot_wider(names_from = Site_ID2, values_from = DistKM)   %>% ## pivot sites to get matrix
  remove_rownames %>% column_to_rownames(var="Site_ID1") ## site 1 names to rows

km_wide <- as.dist(km_wide)
names(allsyncx)


# Leave one out matrix regression - Walter et al 2017 ---------------------
##
# list matrices

mats <- list(bio_wide, temp_wide, div_wide, dist_wide, con_wide)
names(mats) <- c("bio_wide", "temp_wide", "div_wide", "dist_wide", "con_wide")
mats

clean.dat(mats, resp, preds, 2)
lno.score(mats, resp = 1, pred = 2:length(mats), 2, 10)
pred = 2:length(mats)
resp = 1
n=10
maxruns = 10
# Matrix regression -------------------------------------------------------

## make into dist objects
# matrix regression. so far does not include interaction
names(allsyncx)
dim(temp_wide)

## environmental model
MRM(as.dist(bio_wide) ~ as.dist(temp_wide),  nperm=999)

## spatial model
MRM(as.dist(bio_wide) ~ as.dist(con_wide) + as.dist(km_wide),  nperm=999)

## trait diversity model
MRM(as.dist(bio_wide) ~  as.dist(div_wide),  nperm=999)

## combined model

MRM(as.dist(bio_wide) ~ as.dist(temp_wide) + as.dist(con_wide),  nperm=999)
combMat <- MRM(as.dist(bio_wide) ~ as.dist(temp_wide) + as.dist(con_wide) + as.dist(div_wide)+ as.dist(km_wide),  nperm=999)
combMat
?MRM
estsTab <- sjPlot::tab_model(combMat, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "Temp", "Connectivity", "Diversity", "Distance (km)"),
                             dv.labels= "Drivers of Thermal Synchrony")
?tab_model
library(MDMR)
data(mdmrdata)

dim(Y.mdmr)
D <- dist(Y.mdmr, method = "manhattan")
class(Y.mdmr)
View(as.matrix(D))
class(X.mdmr)
names(allsyncx)
newX <- allsyncx %>%  filter(Region == "USA") %>% dplyr::select(distance, diversity, DistKM, annual_avg) %>% as.matrix()
class(newX)
dim(newX)
summary(mdmr.res)

bio_wide[upper.tri(bio_wide)] <- t(bio_wide)[upper.tri(bio_wide)]
X.mdmr
D <- as.dist(bio_wide)
newX <- temp_wide
mdmr.res <- mdmr(X = newX, D = bio_wide)

# Mixed effect models -----------------------------------------------------

# install.packages("effects")
# install.packages("MuMIn")
# library(MuMIn)
library(lme4)
library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions
library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
library("easystats") ## multicolinearality

head(allsyncx)

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "0" = "Between Basin", "1" = "Within Basin") 

mod_mixed = lmer(Sync ~ (distance + annual_avg  + diversity2) *  Connectivity + 
                   (1 + Region | Site_ID1) + ## nested siteID1 within region
                   (1 + Region | Site_ID2), ## nested SiteID2 within region
                    REML = T, ## estimates p vals and F test
                    data = allsyncx)

summary(mod_mixed)
anova(mod_mixed)
r2_nakagawa(mod_mixed) ## throwing NA
# r.squaredGLMM(mod_mixed) ## package not working for R version
check_singularity(mod_mixed) ## False

mod_cor = lmer(Sync ~ (distance + annual_avg  + diversity)  + 
                   (1 | Region/Site_ID1) + ## nested siteID1 within region
                   (1 | Region/Site_ID2), ## nested SiteID 2 within region
                 REML = T, ## estimates p vals and F test
                 data = allsyncx)

check_collinearity(mod_cor) ## all ok, connectivity and interactions not included as not valid

### plots 

ests <- sjPlot::plot_model(mod_mixed, 
                   show.values=TRUE, show.p=TRUE,
                   title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_v2.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=5, width=6)

estsTab <- sjPlot::tab_model(mod_mixed, 
                  show.re.var= TRUE, 
                  pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                  dv.labels= "Drivers of Thermal Synchrony")

estsTab



effects_diversity <- effects::effect(term= "diversity:Connectivity", mod= mod_mixed)
summary(effects_diversity) #output of what the values are


lmmod <-  lm(Sync ~ 1 + as.factor(Region), data=allsyncx)
summary(lmmod)

x_div <- as.data.frame(effects_diversity) 
names(x_div)[3] <- "Sync"

## connectiovity as a factor and change name 
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "0" = "Between Basin", "1" = "Within Basin")

x_div$Connectivity <- as.factor(x_div$Connectivity)
x_div$Connectivity <- recode_factor(x_div$Connectivity, "0" = "Between Basin", "1" = "Within Basin")


div1<-ggplot() +
      geom_point(data = allsyncx, aes(x=diversity, y=Sync, colour =Connectivity), size = 0.5) +
        geom_line(data=x_div, aes(x=diversity, y=Sync, colour =Connectivity)) +
          geom_ribbon(data= x_div, aes(x=diversity, ymin=lower, ymax=upper, fill = Connectivity), alpha= 0.3) +
  # # scale_color_hue(labels = c("Within Basin", "Between Basin")) +
  # scale_color_discrete(labels = c("Within Basin", "Between Basin")) +
  labs(x="Thermal Diversity", y="Trait Synchrony")

div1
                      
file.name1 <- paste0(out.dir, "diversity_sync_con.jpg")
ggsave(div1, filename=file.name1, dpi=300, height=5, width=6)

effects_annualAvg <- effects::effect(term= "annual_avg:Connectivity", mod= mod_mixed)
summary(effects_annualAvg) #output of what the values are

x_div <- as.data.frame(effects_annualAvg)
names(x_div)[3] <- "Sync"
x_div

## connectivity as a factor and change name for plot legend
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "0" = "Between Basin", "1" = "Within Basin")

x_div$Connectivity <- as.factor(x_div$Connectivity)
x_div$Connectivity <- recode_factor(x_div$Connectivity, "0" = "Between Basin", "1" = "Within Basin")

T1 <- ggplot() +
  geom_point(data = allsyncx, aes(x=annual_avg, y=Sync, colour =Connectivity), size = 0.5) +
  geom_line(data=x_div, aes(x=annual_avg, y=Sync, colour =Connectivity)) +
  geom_ribbon(data= x_div, aes(x=annual_avg, ymin=lower, ymax=upper, fill = Connectivity), alpha= 0.3) +
  # # scale_color_hue(labels = c("Within Basin", "Between Basin")) +
  # scale_color_discrete(labels = c("Within Basin", "Between Basin")) +
  labs(x="Temperature Synchrony", y="Trait Synchrony")
  

T1

file.name1 <- paste0(out.dir, "temp_sync_con.jpg")
ggsave(T1, filename=file.name1, dpi=300, height=5, width=6)

effects_distance <- effects::effect(term= "distance:Connectivity", mod= mod_mixed)
summary(effects_diversity) #output of what the values are

x_div <- as.data.frame(effects_distance)
names(x_div)[3] <- "Sync"
x_div

## connectivity as a factor and change name for plot legend
allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)
allsyncx$Connectivity <- recode_factor(allsyncx$Connectivity, "0" = "Between Basin", "1" = "Within Basin")

x_div$Connectivity <- as.factor(x_div$Connectivity)
x_div$Connectivity <- recode_factor(x_div$Connectivity, "0" = "Between Basin", "1" = "Within Basin")

D1 <-ggplot() +
  geom_point(data = allsyncx, aes(x=distance, y=Sync, colour =Connectivity), size = 0.5) +
  geom_line(data=x_div, aes(x=distance, y=Sync, colour =Connectivity)) +
  geom_ribbon(data= x_div, aes(x=distance, ymin=lower, ymax=upper, fill = Connectivity), alpha= 0.3) +
  # # scale_color_hue(labels = c("Within Basin", "Between Basin")) +
  # scale_color_discrete(labels = c("Within Basin", "Between Basin")) +
  labs(x="Thermal Distance", y="Trait Synchrony")

D1

file.name1 <- paste0(out.dir, "distance_sync_con.jpg")
ggsave(D1, filename=file.name1, dpi=300, height=5, width=6)

names(allsyncx)


library(merTools)
# install.packages("merTools")

predictInterval(mod_mixed)   # for various model predictions, possibly with new data

REsim(mod_mixed)             # mean, median and sd of the random effect estimates

plotREsim(REsim(mod_mixed))  # plot the interval estimates

library(car)
install.packages("car")
install.packages("emmeans")
anova(mod_mixed)


# Membership mixed model --------------------------------------------------
# library(devtools)
# install_github("jvparidon/lmerMultiMember")
# install.packages("performance")
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
?rescale
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
check_model(mem_mixed1)

### plots 
class(mem_mixed1) <- "lmerModLmerTest"
# sjPlot::plot_model(mem_mixed) 
ests <- sjPlot::plot_model(mem_mixed1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony")

ests
file.name1 <- paste0(out.dir, "effect_sizes_diverisity2.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

estsTab <- sjPlot::tab_model(mem_mixed1, 
                             show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony")

estsTab

### mixed model with random slopes

mod_mixed = lmer(Sync ~ (distance + annual_avg  + diversity2) *  Connectivity + 
                   (annual_avg | Region) + # random slopes
                   (distance| Region) +
                   (diversity2| Region) +
                   (Connectivity| Region), 
                 REML = T, ## estimates p vals and F test
                 data = allsyncx)

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




# Mixed effects model with mean and spatial cor ----------------------------
library(nlme)
library(car)
library(emmeans)
head(allsyncx)
names(allsyncx)

## calculate means per site
syncMeans <- allsyncx %>%
  pivot_longer(c(Sync, distance, diversity, annual_avg), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable, Region, Trait, Site_ID1, Connectivity) %>%
  summarise(MeanVals = mean(Value)) %>%
  pivot_wider(names_from = Variable, values_from = MeanVals) %>%
  rename(SiteID = Site_ID1)

## upload fish abundance and site data
originaldata <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")
head(originaldata)

## take only sites and basins
sites <- originaldata %>%
  dplyr::select(SiteID, HydroBasin, Latitude, Longitude) %>%
  distinct()

## join coords to main df by site

meanSyncx <- left_join(syncMeans, sites, by = "SiteID")

head(sites)

head(meanSyncx)

d <- cbind(meanSyncx, dummy=rep(1, nrow(meanSyncx)))

# d <- d %>% filter(Connectivity == "Between Basin")

## between/within basins have same coord for each site, make dummy coords by adding 2 to between basin. Then
## spatial correlation will work
d <- d %>% group_by(Connectivity) %>%
  mutate(LongitudeD = ifelse(Connectivity == "Between Basin", Longitude +2, Longitude)) %>%
  mutate(LatitudeD = ifelse(Connectivity == "Between Basin", Latitude +2, Latitude))

## mixed model wth sptail correlation
modS.mn <- lme(fixed = Sync ~ (distance + annual_avg  + diversity) *  Connectivity, 
               data = d, 
               control = lmeControl(opt = "optim"),
               random = ~ 1 | Region,
               correlation = corGaus(1, form = ~ LatitudeD + LongitudeD)) 

summary(modS.mn)
Anova(modS.mn)


### plots 

ests <- sjPlot::plot_model(modS.mn, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of Thermal Synchrony (Spatial Correlation)")

ests
file.name1 <- paste0(out.dir, "effect_sizes_SC.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=5, width=6)

estsTab <- sjPlot::tab_model(modS.mn, 
                             # show.re.var= TRUE, 
                             pred.labels =c("(Intercept)", "DistKM", "annual avg", "diversity", "Connectivity"),
                             dv.labels= "Drivers of Thermal Synchrony (Spatial Correlation)")
?plot_model
estsTab

lmmod <-  lm(Sync ~ 1 + as.factor(Region), data=syncMeans)
summary(lmmod)
