## Figures

# packages
library(tidylog)
library(tidyverse)


getwd()

## set up folder to save figures
out.dir <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V2/figures/"

# All synchrony -----------------------------------------------------------
## data 

load(file = "output_data/sync/03_sync_data_funcgroup_traitgroup_similarity_euclidean_dist_interpolated.RData") ## syncDF
head(SyncBasin)

SyncBasin <- syncDF %>%
  select(-X) %>%
  mutate(Euclid_Dist_KM = Euclid_Dist_Meters/1000) 

# New facet label names for inbasin
# supp.labs <- c("Within Same Basin", "Outside Basin")
# names(supp.labs) <- c(1, 0)


## within basin synchrony per trait group & feeding group

sm1 <- ggplot(subset(SyncBasin, Connectivity == 1), aes(x=Euclid_Dist_KM, y=Correlation, color = Trait)) +
  geom_smooth(method = "gam") +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony", limits=c(-1, 1)) +
  scale_x_continuous(name="Eucliean Distance (km)") 
sm1


file.name1 <- paste0(out.dir, "Within_basin_sync_single_traits.jpg")
ggsave(sm1, filename=file.name1, dpi=300, height=5, width=6)

## between basin synchrony

sm2 <- ggplot(subset(SyncBasin, Connectivity == 0), aes(x=Euclid_Dist_KM, y=Correlation, color = Trait)) +
  geom_smooth(method = "gam") +
  # facet_grid(cols = vars(TraitGroup), rows = vars(FeedingGroup)) +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony", limits=c(-1, 1)) +
  scale_x_continuous(name="Eucliean Distance (km)") 
sm2


file.name2 <- paste0(out.dir, "Between_basin_sync_single_traits.jpg")
ggsave(sm2, filename=file.name2, dpi=300, height=5, width=6)

# reduce distances

SyncBasin600 <- SyncBasin %>%
  filter(Connectivity == 1, Euclid_Dist_KM <= 600)

### log distance, use lm model 

## within basin synchrony per trait group & feeding group

sm1a <- ggplot(subset(SyncBasin600, Connectivity == 1), aes(x=Euclid_Dist_KM, y=Correlation, color = Trait)) +
  geom_smooth(method = "lm") +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony") +
  scale_x_log10(name="Log Eucliean Distance (km)") 
sm1a


file.name1 <- paste0(out.dir, "Within_basin_sync_single_traits_distance_lm_600_logx.jpg")
ggsave(sm1a, filename=file.name1, dpi=300, height=5, width=6)


## between basin synchrony

sm2a <- ggplot(subset(SyncBasin, Connectivity == 0), aes(x=Euclid_Dist_KM, y=Correlation, color = Trait)) +
  geom_smooth(method = "lm") +
  # facet_grid(cols = vars(TraitGroup), rows = vars(FeedingGroup)) +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony") +
  scale_x_log10(name="Log Eucliean Distance (km)")
sm2a


file.name2 <- paste0(out.dir, "Between_basin_sync_single_traits_distance_lm_logx.jpg")
ggsave(sm2a, filename=file.name2, dpi=300, height=5, width=6)

### points to get distribution

## calculate means, too many points for one figure
names(SyncBasin600)

## within basin below 600km
MeanSync600 <- SyncBasin600 %>%
  group_by(Trait, Site_ID2, Connectivity, TraitGroup) %>%
  summarise(Mean_Cor = mean(Correlation), Mean_Dist = mean(Euclid_Dist_KM), Mean_Lat =mean(MeanLat))

s1<- ggplot(MeanSync600, aes(x=Mean_Dist, y=Mean_Cor, color = Trait)) + 
  geom_point() +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony") +
  scale_x_log10(name="Log Eucliean Distance (km)")

s1

file.name1 <- paste0(out.dir, "Mean_within_basin_sync_single_traits_points.jpg")
ggsave(s1, filename=file.name1, dpi=300, height=5, width=6)


## all pairs
MeanSync <- SyncBasin %>%
  group_by(Trait, Site_ID2, Connectivity, TraitGroup) %>%
  summarise(Mean_Cor = mean(Correlation), Mean_Dist = mean(Euclid_Dist_Meters), Mean_Lat =mean(MeanLat))

## between basin
s2 <- ggplot(subset(MeanSync, Connectivity == 0), aes(x=Mean_Dist, y=Mean_Cor, color = Trait)) + 
  geom_point() +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony") +
  scale_x_log10(name="Log Eucliean Distance (km)")

s2

file.name1 <- paste0(out.dir, "Mean_between_basin_sync_single_traits_points.jpg")
ggsave(s2, filename=file.name1, dpi=300, height=5, width=6)


# All synchrony with temp and flow ----------------------------------------

## uses all distances, may need to remove >600km within basin!!!!

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

head(SyncBasinTF)

## within basin synchrony v temp

t1 <- ggplot(subset(SyncBasinTF, Connectivity == 1), aes(x=TempCor, y=Correlation, color = Trait)) +
  geom_smooth(method = "lm") +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony") +
  scale_x_continuous(name="Synchrony of Max Temp") 
t1


file.name1 <- paste0(out.dir, "Within_basin_sync_v_temp_sync_single_traits.jpg")
ggsave(t1, filename=file.name1, dpi=300, height=5, width=6)

## between basin synchrony

t2 <- ggplot(subset(SyncBasinTF, Connectivity == 0), aes(x=TempCor, y=Correlation, color = Trait)) +
  geom_smooth(method = "lm") +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony") +
  scale_x_continuous(name="Synchrony of Max Temp") 
t2


file.name2 <- paste0(out.dir, "Between_basin_sync_v_temp_sync_single_traits.jpg")
ggsave(t2, filename=file.name2, dpi=300, height=5, width=6)

## within basin synchrony v  flow

f1 <- ggplot(subset(SyncBasinTF, Connectivity == 1), aes(x=FlowCor, y=Correlation, color = Trait)) +
  geom_smooth(method = "lm") +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony") +
  scale_x_continuous(name="Synchrony of Max Q") 
f1


file.name1 <- paste0(out.dir, "Within_basin_sync_v_flow_sync_single_traits.jpg")
ggsave(f1, filename=file.name1, dpi=300, height=5, width=6)

## between basin synchrony

f2 <- ggplot(subset(SyncBasinTF, Connectivity == 0), aes(x=FlowCor, y=Correlation, color = Trait)) +
  geom_smooth(method = "lm") +
  facet_wrap(~TraitGroup) +
  scale_y_continuous(name="Synchrony") +
  scale_x_continuous(name="Synchrony of Max Q") 
f2


file.name2 <- paste0(out.dir, "Between_basin_sync_v_flow_sync_single_traits.jpg")
ggsave(f2, filename=file.name2, dpi=300, height=5, width=6)

head(SyncBasinTF)


# All together now... -----------------------------------------------------

## make wide
SyncBasinTF_wide <- SyncBasinTF %>%
  pivot_wider(names_from = "Trait", values_from = "Correlation") %>%
  pivot_longer(TempCor:Q_pref, names_to="Variable", values_to = "Synchrony")
names(SyncBasinTF_wide)


a2 <- ggplot(subset(SyncBasinTF_wide, Connectivity == 1), aes(x=Euclid_Dist_KM, y=Synchrony, color = Variable)) +
  geom_smooth(method = "glm") +
  facet_wrap(~TraitGroup) +
  scale_x_log10(name="Log Eucliean Distance (km)") +
  scale_y_continuous(name="Synchrony")
a2

file.name1 <- paste0(out.dir, "Within_basin_single_traits_sync_temp_flow_lm_log.jpg")
ggsave(a2, filename=file.name1, dpi=300, height=5, width=6)

a3 <- ggplot(subset(SyncBasinTF_wide, Connectivity == 0), aes(x=Euclid_Dist_KM, y=Synchrony, color = Variable)) +
  geom_smooth(method = "glm") +
  facet_wrap(~TraitGroup) +
  scale_x_log10(name="Log Eucliean Distance (km)") +
  scale_y_continuous(name="Synchrony")
a3

file.name1 <- paste0(out.dir, "Between_basin_single_traits_sync_temp_flow_lm_log.jpg")
ggsave(a3, filename=file.name1, dpi=300, height=5, width=6)

unique(SyncBasinTF_wide$Variable)


# Flow and temp figures ---------------------------------------------------

# keep all summary stats

TempSync <- read.csv("output_data/sync/02_temperature_between_all_sites_biogeographic_regions.csv")
head(TempSync)

unique(TempSync$env_var)

## format for join
TempSync <- rename(TempSync, Pair = X)
TempSync <- TempSync %>%
  # filter(env_var == "clim_max_raw") %>%
  select(Correlation, Site_ID1, Site_ID2, env_var) %>%
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  select(-Site_ID2, - Site_ID1) %>%
  mutate(Env = "Temperature" )

FlowSync <- read.csv("output_data/sync/02_flow_between_all_sites_biogeographic_regions.csv")
head(FlowSync)

## format for join
FlowSync <- rename(FlowSync, Pair = X)
FlowSync <- FlowSync %>%
  # filter(env_var == "qmax_raw") %>%
  select(Correlation, Site_ID1, Site_ID2, env_var) %>%
  mutate(Pair = paste0(Site_ID1, ".", Site_ID2)) %>%
  select(-Site_ID2, - Site_ID1) %>% 
  mutate(Env = "Flow")

SyncBasinTF <- rbind(TempSync, FlowSync)
head(SyncBasinTF)

SyncBasinTF <- left_join(SyncBasinTF, SyncBasin, by = "Pair")

### remove traits, keep connectivity and distance etc
SyncBasinTF_wide <- SyncBasinTF %>%
  dplyr::select(-c(Trait:Correlation.y)) %>%
  rename(Synchrony = Correlation.x, Variable = env_var)
  # pivot_wider(names_from = "env_var", values_from = "Correlation.x") %>%
  # pivot_longer(TempCor:Q_pref, names_to="Variable", values_to = "Synchrony")

head(SyncBasinTF_wide)

a2 <- ggplot(subset(SyncBasinTF_wide, Connectivity == 1), aes(x=Euclid_Dist_KM, y=Synchrony, color = Variable)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Env) +
  scale_x_log10(name="Log Eucliean Distance (km)") +
  # scale_x_continuous(name = "Euclidean Distance (km)") +
  scale_y_continuous(name="Synchrony")
a2

file.name1 <- paste0(out.dir, "Within_basin_temp_flow_logx_lm.jpg")
ggsave(a2, filename=file.name1, dpi=300, height=5, width=6)

a3 <- ggplot(subset(SyncBasinTF_wide, Connectivity == 0), aes(x=Euclid_Dist_KM, y=Synchrony, color = Variable)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Env) +
  scale_x_log10(name="Log Eucliean Distance (km)") +
  # scale_x_continuous(name = "Euclidean Distance (km)") +
  scale_y_continuous(name="Synchrony")
a3

file.name1 <- paste0(out.dir, "Between_basin_temp_flow_logx_lm.jpg")
ggsave(a3, filename=file.name1, dpi=300, height=5, width=6)



# Leave one out synchrony -------------------------------------------------

## data 

load(file = "output_data/sync/04_sync_single_traits_LOO.RData") # ssLOO_join

head(ssLOO_join)
names(ssLOO_join)
dim(ssLOO_join)

##  convert distance to KM and year to factor
SyncBasin <- ssLOO_join %>%
  # select(-X) %>%
  mutate(Euclid_Dist_KM = Euclid_Dist_Meters/1000, YearRemoved = as.factor(YearRemoved)) 

# SyncBasin <- SimSync

head(SyncBasin)

## within basin boxplot

smbx <- ggplot(subset(SyncBasin, Connectivity == 1), aes(x=YearRemoved, y=SyncDiff, color = Trait)) +
  geom_boxplot() +
  facet_wrap(~TraitGroup) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(name="Change in Synchrony") 
# scale_x_log10(name="Year Removed")

smbx

file.name1 <- paste0(out.dir, "Within_basin_single_traits_traitgrp_LOO_boxplots.jpg")
ggsave(smbx, filename=file.name1, dpi=300, height=5, width=6)



## between basin

smbx2 <- ggplot(subset(SyncBasin, Connectivity == 0), aes(x=YearRemoved, y=SyncDiff, color = Trait)) +
  geom_boxplot() +
  facet_wrap(~TraitGroup) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_y_continuous(name="Change in Synchrony") 
# scale_x_log10(name="Year Removed")

smbx2

file.name1 <- paste0(out.dir, "Between_basin_single_traits_traitgrp_LOO_boxplots.jpg")
ggsave(smbx2, filename=file.name1, dpi=300, height=5, width=6)

