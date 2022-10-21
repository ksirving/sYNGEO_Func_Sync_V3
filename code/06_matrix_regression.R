### matrix regression

library(tidyverse)
library(tidylog)
library(ecodist)

# install.packages("remotes")
# remotes::install_github("reumandc/mms")
getwd()
# install.packages("PopGenReport")
library(PopGenReport)

load(file = "output_data/sync/04_temp_pref_env_dist_no_dupl_pairs.RData")

out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V3/Figures/"


# Create matrices ---------------------------------------------------------

head(allsyncx)

## some NAs in sync. find and remove - check in sync code later!!!!!
ind <- which(is.na(allsyncx$Sync))
allsyncx[ind,]
allsyncx <- allsyncx[-ind,]
sum(is.na(allsyncx))

allsyncx <- na.omit(allsyncx)
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
class(D)
class(bio_wide)
write.csv(bio_wide, "output_data/temp_pref_sync_matrix_eu.csv")

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


write.csv(temp_wide, "output_data/temp_sync_matrix_eu.csv")

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


write.csv(div_wide, "output_data/temp_diversity_matrix_eu.csv")

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

write.csv(dist_wide, "output_data/temp_diversity_matrix_eu.csv")

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

write.csv(con_wide, "output_data/connectivity_matrix_eu.csv")

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
library(cowplot) #for manuscript ready figures
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions

allsyncx$Connectivity <- as.factor(allsyncx$Connectivity) 

mod_mixed = lmer(Sync ~ (annual_avg  + diversity) *  Connectivity + (1 | Site_ID1), data = allsyncx)

summary(mod_mixed)

confint(mod_mixed)

library(merTools)
# install.packages("merTools")

predictInterval(mod_mixed)   # for various model predictions, possibly with new data

REsim(mod_mixed)             # mean, median and sd of the random effect estimates

plotREsim(REsim(mod_mixed))  # plot the interval estimates

library(car)
install.packages("car")

anova(mod_mixed)

## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
install.packages("lmerTest")
library(lmerTest)
L

Model.REML = lmer(Sync ~ (annual_avg  + diversity) *  Connectivity + (1 | Site_ID1), REML = T, data = allsyncx)

anova(Model.REML) #Performs F test on fixed effects using Satterthwaite approximation

summary(Model.REML) #gives model output with estimated df and p values using Satterthwaite

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



effects_diversity <- effects::effect(term= "diversity", mod= mod_mixed)
summary(effects_diversity) #output of what the values are

x_div <- as.data.frame(effects_diversity)
x_div

allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)

ggplot() + 
  #2
  geom_point(data=allsyncx, aes(diversity, Sync, colour = Connectivity)) + 
  #3
  geom_point(data=x_div, aes(x=diversity, y=fit), color="blue") +
  #4
  geom_line(data=x_div, aes(x=diversity, y=fit), color="blue") +
  #5
  geom_ribbon(data= x_div, aes(x=diversity, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
 
  #6
  labs(x="Thermal Diversity", y="Trait Synchrony")


effects_annualAvg <- effects::effect(term= "annual_avg*Connectivity", mod= mod_mixed)
summary(effects_annualAvg) #output of what the values are

x_div <- as.data.frame(effects_annualAvg)
x_div

allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)

T1 <- ggplot() + 
  #2
  geom_point(data=allsyncx, aes(annual_avg, Sync)) + 
  #3
  geom_point(data=x_div, aes(x=annual_avg, y=fit)) +
  #4
  geom_line(data=x_div, aes(x=annual_avg, y=fit)) +
  #5
  geom_ribbon(data= x_div, aes(x=annual_avg, ymin=lower, ymax=upper), alpha= 0.3, color = "blue") +
  
  facet_wrap(~Connectivity) +
  
  #6
  labs(x="Temperature Synchrony", y="Trait Synchrony")

file.name1 <- paste0(out.dir, "temp_sync_con.jpg")
ggsave(T1, filename=file.name1, dpi=300, height=5, width=6)

effects_diversity <- effects::effect(term= "diversity*Connectivity", mod= mod_mixed)
summary(effects_diversity) #output of what the values are

x_div <- as.data.frame(effects_diversity)
x_div

allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)

D1 <- ggplot() + 
  #2
  geom_point(data=allsyncx, aes(diversity, Sync)) + 
  #3
  geom_point(data=x_div, aes(x=diversity, y=fit)) +
  #4
  geom_line(data=x_div, aes(x=diversity, y=fit)) +
  #5
  geom_ribbon(data= x_div, aes(x=diversity, ymin=lower, ymax=upper), alpha= 0.3, color = "blue") +
  
  facet_wrap(~Connectivity) +
  
  #6
  labs(x="Thermal Diversity", y="Trait Synchrony")

file.name1 <- paste0(out.dir, "diversity_sync_con.jpg")
ggsave(D1, filename=file.name1, dpi=300, height=5, width=6)

effects_diversity <- effects::effect(term= "DistKM*Connectivity", mod= mod_mixed)
summary(effects_diversity) #output of what the values are

x_div <- as.data.frame(effects_diversity)
x_div

allsyncx$Connectivity <- as.factor(allsyncx$Connectivity)

KM1 <- ggplot() + 
  #2
  geom_point(data=allsyncx, aes(DistKM, Sync)) + 
  #3
  geom_point(data=x_div, aes(x=DistKM, y=fit)) +
  #4
  geom_line(data=x_div, aes(x=DistKM, y=fit)) +
  #5
  geom_ribbon(data= x_div, aes(x=DistKM, ymin=lower, ymax=upper), alpha= 0.3, color = "blue") +
  
  facet_wrap(~Connectivity) +
  
  #6
  labs(x="Distance (KM)", y="Trait Synchrony")

KM1

file.name1 <- paste0(out.dir, "Euc_dust_con.jpg")
ggsave(KM1, filename=file.name1, dpi=300, height=5, width=6)
