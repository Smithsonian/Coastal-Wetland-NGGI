## Coastal Wetland NGGI
## wolfejax@si.edu

# workflow prepares synthesis data for use in calculating CAR from sediment accretion rates in reported values table

## Prepare Workspace ####

library(tidyverse)

source("scripts/nggi_utils.R")

# pull cores and depthseries tables

# V1.0.0
cores <- read_csv("https://ndownloader.figshare.com/files/42289539")
ds <- read_csv("https://ndownloader.figshare.com/files/42289545", guess_max = 51000)

# V1.1.0
cores <- read_csv("https://ndownloader.figshare.com/files/43251309", guess_max = 11000)
ds <- read_csv("https://ndownloader.figshare.com/files/43251963", guess_max = 72000)

## Calculate carbon interval stock 

# Gapfill depthseries data and calculate carbon density

ds_stocks <- ds %>%
  left_join(cores %>% select(study_id, site_id, core_id, habitat, country)) %>%
  filter(country == "United States") %>% 
  # transfer representative depths to actual depths (need to adjust this later)
  mutate(depth_min = ifelse(is.na(depth_min) & !is.na(representative_depth_min), 
                            representative_depth_min, depth_min),
         depth_max = ifelse(is.na(depth_max) & !is.na(representative_depth_max), 
                            representative_depth_max, depth_max)) %>% 
  # need to straighten out the habitats
  # this should really happen in the synthesis
  # mutate(habitat = case_when(habitat == "scrub shrub" ~ "scrub/shrub",
  #                            # one of carlin's is actually mudflat
  #                            study_id %in% c("Boyd_et_al_2017", "Watson_and_Byrne_2013", "Thom_1992",
  #                                            "Drexler_et_al_2019", "Carlin_et_al_2021") ~ "marsh",
  #                            core_id == "W4" ~ "marsh",
  #                            core_id %in% c("W1", "W2", "W2") ~ "swamp",
  #                            # "Krauss_et_al_2018" and "Ensign_et_al_2020" is a mix of marsh and swamp
  #                            T ~ habitat)) %>% 
  carbonStock(.)

# nrow(ds_stocks %>% drop_na(carbon_density))/nrow(ds_stocks)
# 66% coverage  - gets better with more habitats assigned
# Nahlic and Fennessey marsh outlier DBD sample_id NWCA11-3583-1-11-MD-019-016-1
# density
ggplot(ds_stocks %>% filter(carbon_density < 0.8), aes(carbon_density, habitat)) + geom_boxplot()
# ggplot(ds_stocks, aes(fraction_organic_matter, fraction_carbon, col = habitat)) + geom_point()
# stocks

# Standardize calculations for specified interval(s)
standard_stocks <- standardizeDepths(ds_stocks)

# standard_stocks %>% 
#   drop_na(stock_gCm2) %>% 
#   mutate(stock_MgHa = stock_gCm2 * (10^4/10^6)) %>% 
#   # filter(stock_MgHa < 1500) %>% 
#   left_join(ds_stocks %>% distinct(study_id, core_id, habitat)) %>% 
#   ggplot(aes(stock_MgHa, habitat)) + geom_boxplot()

# write to intermediate folder for use in the accretionToCAR fxn
# write_csv(standard_stocks, "coastal_wetland_nggi/data/derived/soil_carbon_calculations.csv")

# # calculate whole 1m core stocks and scale to Mg/Ha
# core_stocks <- ds_stocks %>% 
#   drop_na(stock_gCm2) %>% 
#   group_by(study_id, site_id, core_id, habitat) %>%
#   summarize(stock_gCm2 = mean(stock_gCm2)) %>% 
#   #   # convert gC m-2 to MgC ha-1
#   mutate(stock_MgHa = stock_gCm2 * (10^4/10^6))
# 
# ggplot(core_stocks, aes(stock_gCm2, habitat)) + geom_boxplot()


## Reported Values ####

# # read in reported values
# reported <- read_csv("coastal_wetland_nggi/data/original/NGGI_reported_values_20230818 - literature_values.csv",
#                      na = c("---", "NA", "")) %>% 
#   filter(study_id != "Drexler_et_al_2013") %>% # temporary
#   mutate(pb210_rate = as.numeric(pb210_rate)) %>% 
#   # temp recode core_ids to test 
#   mutate(core_id = case_when(study_id == "Allen_et_al_2022" ~ paste0("NRM_", core_id),
#                           T ~ core_id)) %>% 
#   # cols to ignore
#   select(-c(contains("radiocarbon"), contains("SET"), contains("marker"), contains("stock"),
#             notes, reported_rates, publication_table))
#   
# # apply accretion to CAR function
# # the more core_id's we get to match the better this works
# reported_derived <- accretionToCAR(reported)

# Reported values tables needs to straighten out:
# Core IDs (make sure they align with the Data Library)
# Depth intervals/Max depths for given rates
# Summarizing that one core ("2L-Pb")
# Drexler studies (get rid of ranges)
# averaging Pb210 CRS and CIC accretion rates (Luk 2020)
# what to do with site-level rates (multiple cores)
# propagating error

# nohab <- ds_stocks %>% filter(is.na(habitat)) %>% distinct(study_id) %>% pull(study_id)
# nohab[which(nohab %in% unique(reported$study_id))]

