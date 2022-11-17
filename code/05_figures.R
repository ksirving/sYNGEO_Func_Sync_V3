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
getwd()
## directory for figures
out.dir <- "/Users/katieirving/OneDrive - SCCWRP/Documents - Katie’s MacBook Pro/git/sYNGEO_Func_Sync_V3/Figures/"

## raw data

load(file = "output_data/sync/04_temp_pref_env_dist_no_dupl_pairs.RData") #allsyncx
str(allsyncx)

# all fitted values on one plot -------------------------------------------

## define labels
reg.labs <- c("Europe", "Oceania", "USA")
names(reg.labs) <- c("Europe", "Oceania", "USA")

supp.labs <- c("Within Basin", "Between Basin")
names(supp.labs) <- c("1","0")

trait.labs <- c("Body Size", "Temp. Preference")
names(trait.labs) <- c("AVG_MXL","Tp_pref")

allsync <- allsyncx %>%
  mutate(Connectivity = as.factor(Connectivity)) %>%
  mutate(Connectivity = factor(Connectivity, levels = c(1, 0))) %>%
  # pivot_longer(c(distance:diversity, annual_avg), values_to = "Values", names_to = "Variable") %>%
  rename(SyncFit = Sync)

unique(allsync$Variable)


# ggplot(data = subset(allsync, Trait == "Tp_pref"), aes(y=SyncFit, x=Values, colour = Variable)) +
#   geom_smooth(method = "lm") +
#   facet_grid(Region ~ Connectivity , labeller = labeller(Connectivity = supp.labs, Region = reg.labs),
#              scales = "free_x") +
#   scale_colour_discrete(name  ="Variable",
#                         breaks=c("diversity", "distance", "annual_avg"),
#                         labels=c( "Functional Diversity", "Functional Distance", 
#                                   "Annual Temp")) +
#   scale_y_continuous(name="Synchrony") +
#   scale_x_continuous(name="Variables") 



# plotting over distance----------------------------------------------------------------

head(allsync)
unique(allsync$Region)

# allsync <- allsync %>% rename(Sync = SyncFit) %>%
#   mutate(d = exp(d)) ## reverse log transformation

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
  mutate(Connectivity = factor(Connectivity, levels = c(1, 0)))

## trait sync over distance
S1 <- ggplot(allsync, aes(x=DistKM, y=SyncFit, color = Trait)) +
  geom_smooth(method = "lm") +
  # geom_point(size = 0.05) #+
  facet_wrap(vars(Connectivity), labeller = as_labeller(supp.labs),
             scales = "free_x") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Synchrony") +
  scale_x_log10(name="Eucliean Distance (km)", labels = comma) 

S1

file.name1 <- paste0(out.dir, "Sync_Over_Distance.jpg")
ggsave(S1, filename=file.name1, dpi=300, height=5, width=6)


# plotting func sync v env sync -------------------------------------------


## functional sync vs env sync

S2a <- ggplot(allsync,
              aes(x=annual_avg, y=SyncFit, color = Trait)) +
  # geom_point(size = 0.05) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(Connectivity), labeller = as_labeller(supp.labs),
             scales = "free_x") +
  theme(legend.position = "none") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c("Max. Length","Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Trait Synchrony") +
  scale_x_continuous(name="Environmental Synchrony") 

S2a

file.name1 <- paste0(out.dir, "Trait_env_sync.jpg")
ggsave(S2a, filename=file.name1, dpi=300, height=5, width=6)


# Synchrony vs func dist and diversity ------------------------------------

S3ai <- ggplot(allsync,aes(x=distance, y=SyncFit, color = Trait)) +
  geom_smooth(method = "lm") +
  # geom_point(size = 0.05) +
  facet_wrap(vars(Connectivity), labeller = as_labeller(supp.labs),
             scales = "free_x") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c("Max. Length","Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10(name="Functional Synchrony", labels = comma) +
  scale_x_log10(name="Functional distance", labels = comma) 

S3ai

file.name1 <- paste0(out.dir, "trait_sync_dist.jpg")
ggsave(S3ai, filename=file.name1, dpi=300, height=5, width=6)

S3b <-  ggplot(allsync,aes(x=diversity, y=SyncFit, color = Trait)) +
  geom_smooth(method = "lm") +
  # geom_point(size = 0.05) +
  facet_wrap(vars(Connectivity), labeller = as_labeller(supp.labs),
             scales = "free_x") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c("Max. Length","Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_log10(name="Functional Synchrony") +
  scale_x_continuous(name="Functional diversity", labels = comma) 

S3b

file.name1 <- paste0(out.dir, "trait_sync_div.jpg")
ggsave(S3b, filename=file.name1, dpi=300, height=5, width=6)



# Multivariate plot per model --------------------------------------------

head(allsync)
names(allsync)
## format data - make env sync long

allsync<- allsync %>%
  pivot_longer(c(distance, diversity, annual_avg), names_to = "Variable", values_to = "Values") %>%
  rename(Sync = SyncFit) %>%
  mutate(Connectivity = factor(Connectivity, levels = c(1, 0)))


## define labels

supp.labs <- c("Within Basin", "Between Basin")
names(supp.labs) <- c("1","0")


var.labs <- c("Env. Sync", "Func. Diversity", "Func. Distance")
names(var.labs) <- c("annual_avg","diversity" ,"distance" )

var.labs
head(allsync)

M1 <- ggplot(allsync,
              aes(x=Values, y=Sync, color = Variable)) +
  # geom_point(alpha = .6) +
  geom_smooth(method = "lm") +
  facet_grid(Connectivity ~ Variable , labeller = labeller(Connectivity = supp.labs, Variable = var.labs),
             scales = "free_x") +
  # theme(legend.position = "none") +
  # scale_color_discrete(name = "Trait", 
  #                      labels = c("Max. Length","Temp Preference")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(name="Functional Synchrony") +
  scale_x_continuous(name="Predictors") 

M1

file.name1 <- paste0(out.dir, "trait_all_vars.jpg")
ggsave(M1, filename=file.name1, dpi=300, height=5, width=6)



# Boxplots ----------------------------------------------------------------

head(allsync)

## scale predictors

allsyncScaled <- allsync %>%
  group_by(Region, Connectivity, Variable) %>%
  filter( !Site_ID1 == Site_ID2)  %>% ## remove pairs comprised of same sites
  pivot_wider(names_from = Variable, values_from = Values) %>%
  mutate(DivScaled = (diversity-min(na.omit(diversity)))/
           (max(na.omit(diversity))-min(na.omit(diversity)))) %>%
  mutate(DistScaled = (distance-min(na.omit(distance)))/
           (max(na.omit(distance))-min(na.omit(distance)))) %>%
  pivot_longer(c(Sync, annual_avg, distance, diversity, DivScaled, DistScaled), names_to = "Variable", values_to = "Values") 
  # mutate(ValueScaled = (Values-min(na.omit(Values)))/
  #          (max(na.omit(Values))-min(na.omit(Values))))  

# test <- allsync %>%
#   filter(Region == "Europe", Connectivity ==1, Variable == "annual_avg")

test

allsyncScaled

b1 <- ggplot(filter(allsyncScaled, !Variable %in%  c("DivScaled", "DistScaled")), aes(x = Variable, y = Values, colour = Variable)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(Connectivity~Region, labeller = labeller(Connectivity = supp.labs, Variable = var.labs))

b1

b2 <- ggplot(allsyncScaled, aes(x = Variable, y = Values, colour = Variable)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(Connectivity~Region, labeller = labeller(Connectivity = supp.labs, Variable = var.labs))

b2


b3 <- ggplot(filter(allsyncScaled, !Variable %in% c("diversity", "distance")), aes(x = Variable, y = Values, colour = Variable)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(Connectivity~Region, labeller = labeller(Connectivity = supp.labs, Variable = var.labs))

b3

file.name1 <- paste0(out.dir, "boxplots.jpg")
ggsave(b3, filename=file.name1, dpi=300, height=5, width=6)


# Distance between sites --------------------------------------------------


