### testing temp preferences vs site temp

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
library(here)
library(tidylog)

# upload temp preferences
trt1 <- read.csv("input_data/Bio/matSpecies_imputedCloseTaxa.csv")

head(trt1)
dim(trt1) ## 272, 21 traits

## only using temp pref
## remove species with no trait data

trt <- trt1 %>%
  dplyr::select(Species, Tp_pref) %>%
  drop_na()

dim(trt) ## 271
trt

## upload fish abundances

fish_ab <- read.csv("input_data/Bio/fishdata_selection_basins_same_time_window_10262020.csv")

head(fish_ab)
str(fish_ab)


## upload raw temp data

melt_clim_raw <- read.csv(file="input_data/Env/air_annual_and_summer_avg.csv") %>%
  rename(SiteID = siteid) %>% mutate(site_year = paste0(SiteID, "_", year))
head(melt_clim_raw)


range(trt$Tp_pref) ## 12.25 32.06,
mean(trt$Tp_pref) ## 25.19255
median(trt$Tp_pref) ## 25.42


range(na.omit(melt_clim_raw$annual_avg)) ## -2.85 28.25
mean(na.omit(melt_clim_raw$annual_avg)) ## 9.192402
median(na.omit(melt_clim_raw$annual_avg)) ## 9.275

## plot cwm and temp by site_year

temp_cwm <- read.csv(here("output_data/01_trt_single_traits_interpolated_cwm_cmv.csv"))
head(temp_cwm)

length(unique(temp_cwm$site_year))
length(unique(melt_clim_raw$site_year))
## join by site year

temps <- inner_join(temp_cwm, melt_clim_raw, by = "site_year") %>%
  mutate(BiogeoRegion = case_when(Country %in% c("FIN", "SWE", "GBR", "ESP", "FRA") ~ "Europe",
                                  Country== "AUS" ~ "Oceania", Country == "USA" ~ "USA"))
head(temps)

s1 <- ggplot(temps, aes(x=annual_avg, y=CWM, color = BiogeoRegion)) +
  geom_point() #+
  # facet_wrap(~Country)
s1
file.name1 <- paste0(out.dir, "pref_site_temp_scatter.jpg")
ggsave(s1, filename=file.name1, dpi=300, height=5, width=6)

temps_long <- temps %>%
  pivot_longer(c(CWM, annual_avg), names_to = "Variable", values_to = "Values")

## boxplot

b3 <- ggplot(temps_long, aes(x = Variable, y = Values, colour = Variable)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  facet_wrap(~BiogeoRegion)

b3

file.name1 <- paste0(out.dir, "pref_site_temp_boxplots.jpg")
ggsave(b3, filename=file.name1, dpi=300, height=5, width=6)
