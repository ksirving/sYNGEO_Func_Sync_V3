# Figures
# 1 - synchrony over distance - temp, traits
# 2 - func sync v env sync
# 3 - func sync v func distance
# 4 - func sync v func diversity

# packages

library(tidyverse)
library(tidylog)
library("easystats")
library(scales)

## directory for figures
out.dir <- "/Users/katieirving/Documents/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V2/Figures/"

## raw data

# load(file = "output_data/sync/04a_all_sync_for_figures.RData") #allsync

## fitted data

load(file = "output_data/sync/04a_fitted_sync_values_for_figures.RData") ## allsync

head(allsync)

# all fitted values on one plot -------------------------------------------

## define labels
reg.labs <- c("Europe", "Oceania", "USA")
names(reg.labs) <- c("Europe", "Oceania", "USA")

supp.labs <- c("Within Basin", "Between Basin")
names(supp.labs) <- c("1","0")

trait.labs <- c("Body Size", "Temp. Preference")
names(trait.labs) <- c("AVG_MXL","Tp_pref")

allsync <- allsync %>%
  mutate(Connectivity = factor(c, levels = c(1, 0)))

unique(allsync$Variable)

ggplot(data = subset(allsync, Trait == "AVG_MXL"), aes(y=SyncFit, x=Values, colour = Variable)) +
  geom_smooth(method = "lm") +
  facet_grid(Region ~ Connectivity , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
             scales = "free_x") +
  scale_colour_discrete(name  ="Variable",
                        breaks=c("dv_Norm", "di_Norm", "tav"),
                        labels=c( "Functional Diversity", "Functional Distance", 
                                  "Annual Temp")) +
  scale_y_continuous(name="Synchrony") +
  scale_x_continuous(name="Variable") 

ggplot(data = subset(allsync, Trait == "Tp_pref"), aes(y=SyncFit, x=Values, colour = Variable)) +
  geom_smooth(method = "lm") +
  facet_grid(Region ~ Connectivity , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
             scales = "free_x") +
  scale_colour_discrete(name  ="Variable",
                        breaks=c("dv_Norm", "di_Norm", "tav"),
                        labels=c( "Functional Diversity", "Functional Distance", 
                                  "Annual Temp")) +
  scale_y_continuous(name="Synchrony") +
  scale_x_continuous(name="Variables") 

## all together

ggplot(data = allsync, aes(y=SyncFit, x=Values, colour = Variable)) +
  geom_smooth(method = "lm") +
  facet_grid(Trait ~ Connectivity , labeller = labeller(Connectivity = supp.labs, Region = trait.labs),
             scales = "free_x") +
  scale_colour_discrete(name  ="Variable",
                        breaks=c("dv_Norm", "di_Norm", "tav"),
                        labels=c( "Functional Diversity", "Functional Distance", 
                                  "Annual Temp")) +
  scale_y_continuous(name="Synchrony") +
  scale_x_continuous(name="Variable") 

# ggplot(data = subset(allsync, Trait == "Tp_pref"), aes(y=SyncFit, x=Values, colour = Variable)) +
#   geom_smooth(method = "lm") +
#   facet_grid(Region ~ Connectivity , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
#              scales = "free_x") +
#   scale_colour_discrete(name  ="Variable",
#                         breaks=c("dv_Norm", "di_Norm", "tav"),
#                         labels=c( "Functional Diversity", "Functional Distance", 
#                                   "Annual Temp")) +
#   scale_y_continuous(name="Synchrony") +
#   scale_x_continuous(name="Variables") 


# plotting over distance----------------------------------------------------------------

head(allsync)
unique(allsync$Region)

allsync <- allsync %>% rename(Sync = SyncFit) %>%
  mutate(d = exp(d)) ## reverse log transformation

## summary stats
# syncVals <- allsync %>%
#   group_by(Region, Trait, Connectivity) %>%
#   summarise(minSync = min(na.omit(Sync)),
#             maxSync = max(na.omit(Sync)),
#             meanSync = mean(na.omit(Sync)),
#             medianSync = median(na.omit(Sync)))

## define labels
reg.labs <- c("Europe", "Oceania", "USA", "All")
names(reg.labs) <- c("Europe", "Oceania", "USA", "All")

supp.labs <- c("Within Basin", "Between Basin")
names(supp.labs) <- c("1","0")

allsync <- allsync %>%
  mutate(Connectivity = factor(c, levels = c(1, 0)))

## trait sync over distance
S1 <- ggplot(allsync, aes(x=d, y=Sync, color = Trait)) +
  geom_smooth(method = "lm") +
  geom_point(size = 0.05) +
  facet_grid(Connectivity ~ Region , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
             scales = "free_x") +
  theme(legend.position = "none") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c( "Max. Length", "Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Synchrony") +
  scale_x_log10(name="Eucliean Distance (km)", labels = comma) 

S1

file.name1 <- paste0(out.dir, "Trait_sync_distance_regions_fitted.jpg")
ggsave(S1, filename=file.name1, dpi=300, height=5, width=6)

## trait sync over distance
# S1i <- ggplot(allsync, aes(x=d, y=Sync, color = Trait)) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~Connectivity, labeller = as_labeller(supp.labs)) +
#   # facet_grid(Region ~Connectivity , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
#   #            scales = "free_x") +
#   scale_color_discrete(name = "Trait", 
#                        labels = c( "Max. Length", "Temp Preference")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_y_continuous(name="Synchrony") +
#   scale_x_log10(name="Eucliean Distance (km)", labels = comma) 
# 
# S1i
# 
# file.name1 <- paste0(out.dir, "Trait_sync_distance_all_fitted.jpg")
# ggsave(S1i, filename=file.name1, dpi=300, height=5, width=6)

## make all vars long for all on one plot

head(allsync)

allsync_one <- allsync %>%
  pivot_wider(names_from = Variable, values_from = Values) %>%
  pivot_wider(names_from = Trait, values_from = Sync) %>%
  pivot_longer(c(tav, AVG_MXL, Tp_pref), names_to="Variable", values_to = "Sync")  %>%
  mutate(Variable = factor(Variable, levels = c( "AVG_MXL", "Tp_pref", "tav"))) #%>%
  # mutate(Type = case_when(Variable %in% c("annual_avg") ~ "Environment",
  #                         # Variable %in% c("qmean_raw", "qmax_raw", "qmin_raw") ~ "Flow",
  #                         Variable %in% c( "AVG_MXL", "Tp_pref") ~ "Trait"))



# ## all sync over distance
# S1a <- ggplot(allsync_one, aes(x=d, y=Sync, color = Variable)) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~Connectivity, labeller = as_labeller(supp.labs)) +
#   scale_colour_discrete(name  ="Variable",
#                         breaks=c("AVG_MXL", "Tp_pref", "annual_avg"),
#                         labels=c( "Max. Length", "Temp Preference", 
#                                  "Annual Temp")) +
#   # scale_linetype_discrete(name  ="Type",
#   #                         breaks=c("Trait", "Temperature"),
#   #                         labels=c("Trait", "Temperature")) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_y_continuous(name="Synchrony") +
#   scale_x_log10(name="Eucliean Distance (km)", labels = comma) 
# 
# S1a
# 
# 
# file.name1 <- paste0(out.dir, "All_sync_distance_all_fitted.jpg")
# ggsave(S1a, filename=file.name1, dpi=300, height=5, width=6)

## all sync over distance - per region
S1ai <- ggplot(allsync_one, aes(x=d, y=Sync, color = Variable)) +
  geom_smooth(method = "lm") +
  facet_grid(Region ~ Connectivity , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
             scales = "free_x") +
  scale_colour_discrete(name  ="Variable",
                        breaks=c("AVG_MXL", "Tp_pref", "annual_avg"),
                        labels=c( "Max. Length", "Temp Preference", 
                                  "Annual Temp")) +
  # scale_linetype_discrete(name  ="Type",
  #                         breaks=c("Trait", "Temperature"),
  #                         labels=c("Trait", "Temperature")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Synchrony") +
  scale_x_log10(name="Eucliean Distance (km)", labels = comma) 

S1ai


file.name1 <- paste0(out.dir, "All_sync_distance_regions_fitted.jpg")
ggsave(S1ai, filename=file.name1, dpi=300, height=5, width=6)

# plotting func sync v env sync -------------------------------------------

load(file = "output_data/sync/04a_fitted_sync_values_for_figures.RData") ## allsync

head(allsync)
## format data - make env sync long

allsync_wide <- allsync %>%
  pivot_wider(names_from = Variable, values_from = Values) %>%
  rename(Sync = SyncFit)

unique(allsync_wide$Trait)
head(allsync_wide)

## define labels
reg.labs <- c("Europe", "Oceania", "USA", "All")
names(reg.labs) <- c("Europe", "Oceania", "USA", "All")

supp.labs <- c("Within Basin", "Between Basin")
names(supp.labs) <- c("1","0")

allsync_wide <- allsync_wide %>%
  mutate(Connectivity = factor(c, levels = c(1, 0)))

temp.labs <- c("Annual Temp", "Summer Temp", "Annual Mean Q")
names(temp.labs) <- c("annual_avg", "summer_avg", "qmean_raw")

# flow.labs <- c( "Annual Mean Q", "Annual Max Q", "Annual Min Q")
# names(flow.labs) <- c( "qmean_raw", "qmax_raw", "qmin_raw")

## functional sync vs env sync

S2a <- ggplot(allsync_wide,
              aes(x=tav, y=Sync, color = Trait)) +
  geom_point(size = 0.05) +
  geom_smooth(method = "lm") +
  facet_grid(Connectivity ~ Region , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
             scales = "free_x") +
  theme(legend.position = "none") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c("Max. Length","Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Functional Synchrony") +
  scale_x_continuous(name="Environmental Synchrony") 

S2a

file.name1 <- paste0(out.dir, "Trait_env_sync_region.jpg")
ggsave(S2a, filename=file.name1, dpi=300, height=5, width=6)

S2a <- ggplot(allsync_wide,
              aes(x=tav, y=Sync, color = Connectivity)) +
  geom_point(size = 0.05) +
  geom_smooth(method = "lm") +
  facet_grid(Trait ~ Region , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
             scales = "free_x") +
  # theme(legend.position = "none") +
  scale_color_discrete(name = "Connectivity",
                       labels = c("Within Basin","Between Basine")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Functional Synchrony") +
  scale_x_continuous(name="Environmental Synchrony") 

S2a

# file.name1 <- paste0(out.dir, "Trait_env_sync_region.jpg")
# ggsave(S2a, filename=file.name1, dpi=300, height=5, width=6)

# S2ai <- ggplot(allsync_wide,
#               aes(x=tav, y=Sync, color = Trait)) +
#   geom_point(size = 0.05) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~Connectivity, labeller = as_labeller(supp.labs)) +
#   scale_color_discrete(name = "Trait", 
#                        labels = c("Max. Length","Temp Preference")) +
#   scale_y_continuous(name="Functional Synchrony") +
#   scale_x_continuous(name="Environmental Synchrony") 
# 
# S2ai
# 
# file.name1 <- paste0(out.dir, "Trait_env_sync_all.jpg")
# ggsave(S2ai, filename=file.name1, dpi=300, height=5, width=6)


# Synchrony vs func dist and diversity ------------------------------------

head(allsync_wide)

## for normalising
min_max_norm <- function(x) {
  (x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x)))
}


# allsyncLOG <- allsync %>%
#   group_by(Region, Trait) %>%
#   mutate(dv = diversity,
#          di = log(distance+1))  %>%
#   mutate(di_Norm = min_max_norm(di), dv_Norm = min_max_norm(dv)) ## normalise diversity & distance


## define labels
trait.labs <- c( "Max. Length", "Temp Preference")
names(trait.labs) <- c( "AVG_MXL",   "Tp_pref")

# S3a <- ggplot(allsync_wide,aes(x=di_Norm, y=Sync, color = Trait)) +
#   geom_smooth(method = "lm") +
#   geom_point(size = 0.05) +
#   facet_wrap(~Connectivity, labeller = as_labeller(supp.labs)) +
#   # theme(legend.position = "none") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_color_discrete(name = "Trait", 
#                        labels = c("Max. Length","Temp Preference")) +
#   scale_y_continuous(name="Functional Synchrony") +
#   scale_x_continuous(name="Functional distance", labels = comma) 
# 
# S3a
# 
# file.name1 <- paste0(out.dir, "func_sync_dist_all.jpg")
# ggsave(S3a, filename=file.name1, dpi=300, height=5, width=6)

## define labels
reg.labs <- c("Europe", "Oceania", "USA", "All")
names(reg.labs) <- c("Europe", "Oceania", "USA", "All")



S3ai <- ggplot(allsync_wide,aes(x=di_Norm, y=Sync, color = Trait)) +
  geom_smooth(method = "lm") +
  # geom_point(size = 0.05) +
  facet_grid(Connectivity ~ Region , labeller = labeller(Connectivity = supp.labs, Region = reg.labs)) +
  theme(legend.position = "none") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c("Max. Length","Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10(name="Functional Synchrony", labels = comma) +
  scale_x_log10(name="Functional distance", labels = comma) 

S3ai

file.name1 <- paste0(out.dir, "func_sync_dist_region.jpg")
ggsave(S3ai, filename=file.name1, dpi=300, height=5, width=6)

S3b <-  ggplot(allsync_wide,aes(x=dv_Norm, y=Sync, color = Trait)) +
  geom_smooth(method = "lm") +
  # geom_point(size = 0.05) +
  facet_grid(Connectivity ~ Region, labeller = labeller(Connectivity = supp.labs, Region = reg.labs)) +
  theme(legend.position = "none") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c("Max. Length","Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10(name="Functional Synchrony") +
  scale_x_continuous(name="Functional diversity", labels = comma) 

S3b

file.name1 <- paste0(out.dir, "func_sync_div_log_region.jpg")
ggsave(S3b, filename=file.name1, dpi=300, height=5, width=6)

# S3bi <-  ggplot(allsync_wide,aes(x=dv_Norm, y=Sync, color = Trait)) +
#   geom_smooth(method = "lm") +
#   geom_point(size = 0.05) +
#   facet_wrap(~Connectivity, labeller = as_labeller(supp.labs)) +
#   # theme(legend.position = "none") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_color_discrete(name = "Trait", 
#                        labels = c("Max. Length","Temp Preference")) +
#   scale_y_continuous(name="Functional Synchrony") +
#   scale_x_continuous(name="Functional diversity", labels = comma) 
# 
# S3bi
# 
# file.name1 <- paste0(out.dir, "func_sync_div_log_all.jpg")
# ggsave(S3bi, filename=file.name1, dpi=300, height=5, width=6)


head(allsync)




# Multivariate plot per model --------------------------------------------

head(allsync_wide)

load(file = "output_data/sync/04a_fitted_sync_values_for_figures.RData") ## allsync

head(allsync)
unique(allsync$Variable)
## format data - make env sync long

allsync<- allsync %>%
  pivot_wider(names_from = Variable, values_from = Values) %>%
  pivot_longer(d:di_Norm, names_to = "Variable", values_to = "Values") %>%
  rename(Sync = SyncFit) %>%
  mutate(Connectivity = factor(c, levels = c(1, 0)))


## define labels
reg.labs <- c("Europe", "Oceania", "USA")
names(reg.labs) <- c("Europe", "Oceania", "USA")

supp.labs <- c("Within Basin", "Between Basin")
names(supp.labs) <- c("1","0")

trait.labs <- c("Body Size", "Temp. Preference")
names(trait.labs) <- c("AVG_MXL","Tp_pref")

var.labs <- c("Env. Sync", "Func. Diversity", "Func. Distance", "Euc. Distance")
names(var.labs) <- c("tav","dv_Norm" ,"di_Norm", "d" )

var.labs
head(allsync)

M1 <- ggplot(filter(allsync, Trait == "AVG_MXL"),
              aes(x=Values, y=Sync, color = Region)) +
  # geom_point(alpha = .6) +
  geom_smooth(method = "lm") +
  facet_grid(Connectivity ~ Variable , labeller = labeller(Connectivity = supp.labs, Variable = var.labs),
             scales = "free_x") +
  # theme(legend.position = "none") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c("Max. Length","Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Functional Synchrony") +
  scale_x_continuous(name="Environmental Synchrony") 

M1

file.name1 <- paste0(out.dir, "body_size_all_preds_regions.jpg")
ggsave(M1, filename=file.name1, dpi=300, height=5, width=6)

