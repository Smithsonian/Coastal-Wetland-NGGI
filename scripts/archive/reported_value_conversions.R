## CCN-Data-Analytics
# Workflow for standardizing C accumulation values for the 2022 NGGI
# contact: Henry Betts, BettsH@si.edu



## Prepare workspace ####
library(tidyverse)
library(readxl)

# source synthesis
source("resources/pull_synthesis.R")
source("coastal_wetland_nggi/scripts/nggi_utils.R")

# read in reported values
reported_raw <- read_csv("coastal_wetland_nggi/data/original/NGGI_reported_values_20230818 - literature_values.csv")

# read in data from previous NGGI
mengdat <- read_csv("coastal_wetland_nggi/data/original/US-BC-Analysis-1-105.csv")

# read in synthesis data
cores <- getSynthesisData("cores")
depth <- getSynthesisData("depthseries")



## Meng data ####
# Prepare Meng's 2017 synthesis data for merge with reported values
meng <- mengdat %>% 
  rename(latitude = Latitude,
         longitude = Longitude,
         habitat = Ecosystem,
         ecosystem = CCAP_Class,
         climate_zone = Climate_Zone,
         carbon_stock = SOC1, 
         carbon_stock_unit = SOC1units,
         depth_max = SOC1depth, 
         salinity_class = Salinity_Regime, 
         pb210_rate = delSOC1Pb, 
         pb210_rate_unit = delSOC1units,
         cs137_rate = delSOC2Cs,
         cs137_rate_unit = delSOC2units) %>% 
  mutate(study_id = paste(First_Author, Year, sep = "_"),
         pb210_rate = as.numeric(pb210_rate),
         salinity_class = case_when(salinity_class == "0-0.3" ~ "fresh",
                                    salinity_class == "0-0.4" ~ "fresh",
                                    salinity_class == "<0.5" ~ "fresh",
                                    salinity_class == "0-15" ~ "mesohaline",
                                    salinity_class == ">18" ~ "polyhaline",
                                    salinity_class == "<18" ~ "mesohaline",
                                    salinity_class == "20-35" ~ "mixoeuhaline",
                                    salinity_class == "oligo" ~ "oligiohaline",
                                    salinity_class == "poly" ~ "polyhaline",
                                    salinity_class == "meso" ~ "mesohaline",
                                    T ~ salinity_class),
         depth_min = ifelse(!is.na(depth_max), 0, NA_character_)) 
 


## Prepare climate zone list ####
climate_zones <- read_xlsx("coastal_wetland_nggi/data/original/Climate zones.xlsx") %>% 
  mutate(State = recode(State, "Delamare" = "Delaware")) %>% 
  rename(admin_division = State, 
         climate_zone = `Climate zone`)

# isolate US cores from the library with dating information and join climate zones
dated_us_cores <- cores %>% 
  filter(country == "United States") %>% 
  drop_na(dates_qual_code) %>% # isolate dated cores
  left_join(climate_zones) %>% 
  mutate(habitat = recode(habitat, "scrub shrub" = "scrub/shrub"),
         vegetation_class = case_when(core_id == "MR3" ~ "marsh",
                                      T ~ vegetation_class))



## Standardize reported values ####

# look for disaggregated core data. average across depth intervals. (McTigue_et_al_2020: 2L-Pb)
# reported_raw %>%
#   drop_na(core_id) %>% add_count(site_id, core_id) %>% filter(n > 1) 
reported_avg <- reported_raw %>% 
  filter(core_id == "2L-Pb") %>% 
  mutate(core_id = "TB_2L") %>% 
  group_by(study_id, site_id, core_id, pb210_rate_unit) %>%
  summarise(pb210_rate = mean(as.numeric(pb210_rate)),
            pb210_rate_se = sd(as.numeric(pb210_rate))) 

# Drexler_et_al_2013 has pb210_rate as a range - remove and replace with summarized data
reported_drex <- reported_raw %>% 
  filter(grepl("Drexler_et_al_2013", study_id)) %>% 
  separate(pb210_rate, into = c("pb210_rate", "pb210_rate_max"), sep = "-", convert = T) %>% 
  separate(depth_max_cm, into = c("depth_min_cm", "depth_max_cm"), sep = "-", convert = T) %>% 
  mutate(pb210_rate = (as.numeric(pb210_rate) + as.numeric(pb210_rate_max))/2) %>% 
  mutate_at(c("pb210_rate_se", "cs137_rate", "cs137_rate_se", "carbon_stock"), as.numeric)


reported <- reported_raw %>%
  separate(carbon_stock, into = c("carbon_stock", "carbon_stock_se"), sep = " [+][/][-] | [+][/][-]", convert = T) %>% 

  # remove the disaggregated core and replace with summary data
  filter(!grepl("2L-Pb", core_id)) %>%  
  # Drexler_et_al_2013 has pb210_rate as a range - remove and replace with summarized data
  filter(!grepl("Drexler_et_al_2013", study_id)) %>% 
  
  mutate_at(c("pb210_rate", "pb210_rate_se", "cs137_rate", "cs137_rate_se", "depth_min_cm", "depth_max_cm"), as.numeric) %>% 
  full_join(reported_avg) %>% 
  full_join(reported_drex) %>% 
  
  # only work with observations of Pb210 or Cs137 rate data
  filter(!is.na(pb210_rate) | !is.na(cs137_rate)) %>% 
  
  # when CRS and CIC Pb210 rates are available, replace with their average
  mutate(pb210_rate = ifelse(!is.na(pb210_rate_cic & pb210_rate), (pb210_rate_cic + pb210_rate)/2, 
                             ifelse(!is.na(pb210_rate_cic), pb210_rate_cic, pb210_rate)),
         pb210_rate_se = ifelse(!is.na(pb210_rate_cic_se & pb210_rate_se), (pb210_rate_cic_se + pb210_rate_se)/2, 
                                ifelse(!is.na(pb210_rate_cic_se), pb210_rate_cic_se, pb210_rate_se)),
         
         # convert OM to OC depending on habitat type
         carbon_stock = case_when(habitat == "marsh" & accumulation_type == "carbon accumulation" ~ 0.486 * carbon_stock - 0.0160 * carbon_stock^2,
                                  habitat == "mangrove | scrub/shrub | swamp" & accumulation_type == "carbon accumulation" ~ 0.427 * carbon_stock + 0.0635 * carbon_stock^2,
                                  T ~ carbon_stock),
         carbon_stock_se = ifelse(habitat == "marsh" & accumulation_type == "carbon accumulation", 0.486 * carbon_stock_se - 0.0160 * carbon_stock_se^2,
                                  ifelse(habitat == "mangrove | scrub/shrub | swamp" & accumulation_type == "carbon accumulation", 0.427 * carbon_stock_se + 0.0635 * carbon_stock_se^2,
                                         carbon_stock_se)),
         
         # convert pb210 and cs137 rates to cm/yr (accretion rate) or g/m2/yr (accumulation rate), and carbon stock to Mg/ha
         carbon_stock = case_when(carbon_stock_unit == "gramsPerSquareMeter" ~ carbon_stock * .01, 
                                  carbon_stock_unit == "kilogramsPerSquareMeter" ~ carbon_stock * 10,
                                  T ~ carbon_stock),
         carbon_stock_se = case_when(carbon_stock_unit == "gramsPerSquareMeter" ~ carbon_stock_se * .01, 
                                  carbon_stock_unit == "kilogramsPerSquareMeter" ~ carbon_stock_se * 10,
                                  T ~ carbon_stock_se),
         pb210_rate = case_when(pb210_rate_unit == "millimeterPerYear" ~ pb210_rate * .1, 
                                pb210_rate_unit == "gramsPerSquareCentimeterPerYear" ~ pb210_rate * 0.0001,
                                pb210_rate_unit == "kilogramsPerSquareMeterPerYear" ~ pb210_rate * 1000,
                                T ~ pb210_rate),
         pb210_rate_se = case_when(pb210_rate_unit == "millimeterPerYear" ~ pb210_rate_se * .1, 
                                pb210_rate_unit == "gramsPerSquareCentimeterPerYear" ~ pb210_rate_se * 0.0001,
                                pb210_rate_unit == "kilogramsPerSquareMeterPerYear" ~ pb210_rate_se * 1000,
                                T ~ pb210_rate_se),
         cs137_rate = case_when(cs137_rate_unit == "millimeterPerYear" ~ cs137_rate * .1,
                                cs137_rate_unit == "gramsPerSquareCentimeterPerYear" ~ cs137_rate * 0.0001,
                                cs137_rate_unit == "kilogramsPerSquareMeterPerYear" ~ cs137_rate * 1000,
                                T ~ cs137_rate),
         cs137_rate_se = case_when(cs137_rate_unit == "millimeterPerYear" ~ cs137_rate_se * .1,
                                cs137_rate_unit == "gramsPerSquareCentimeterPerYear" ~ cs137_rate_se * 0.0001,
                                cs137_rate_unit == "kilogramsPerSquareMeterPerYear" ~ cs137_rate_se * 1000,
                                T ~ cs137_rate_se),
         
         # make accumulation/accretion flag for cs137 and pb210 to convert accretion rate to accumulation rate after depthseries table merge 
         cs_accretion = ifelse(cs137_rate_unit == "millimeterPerYear", "accretion",
                                         ifelse(cs137_rate_unit == "centimeterPerYear", "accretion", NA_character_)),
         cs_accumulation = ifelse(cs137_rate_unit == "gramsPerSquareMeterPerYear", "accumulation",
                                  ifelse(cs137_rate_unit == "gramsPerSquareCentimeterPerYear", "accumulation",
                                         ifelse(cs137_rate_unit == "kilogramsPerSquareMeterPerYear", "accumulation", NA_character_))),
         pb_accretion = ifelse(pb210_rate_unit == "millimeterPerYear", "accretion",
                               ifelse(pb210_rate_unit == "centimeterPerYear", "accretion", NA_character_)),
         pb_accumulation = ifelse(pb210_rate_unit == "gramsPerSquareCentimeterPerYear", "accumulation",
                                  ifelse(pb210_rate_unit == "gramsPerSquareMeterPerYear", "accumulation",
                                         ifelse(pb210_rate_unit == "kilogramsPerSquareMeterPerYear", "accumulation", NA_character_))),
         
         # create column for accumulation rates
         pb_CAR = ifelse(!is.na(pb_accumulation), pb210_rate,
                         ifelse(!is.na(pb_accretion) & !is.na(carbon_stock), pb210_rate * carbon_stock, 
                                NA_character_)),
         pb_CAR_se = ifelse(!is.na(pb_accumulation), pb210_rate_se,
                         ifelse(!is.na(pb_accretion) & !is.na(carbon_stock), pb210_rate_se * carbon_stock, 
                                NA_character_)),
         cs_CAR = ifelse(!is.na(cs_accumulation), cs137_rate,
                         ifelse(!is.na(cs_accretion) & !is.na(carbon_stock), cs137_rate * carbon_stock, 
                                NA_character_)),
         cs_CAR_se = ifelse(!is.na(pb_accumulation), cs137_rate_se,
                            ifelse(!is.na(cs_accretion) & !is.na(carbon_stock), cs137_rate_se * carbon_stock, 
                                   NA_character_))) %>% 
  mutate_at(c("pb_CAR", "pb_CAR_se", "cs_CAR", "cs_CAR_se"), as.numeric) 


## Align dataframe to match depthseries data from library
reported_merge <- reported %>% 
  
  # remove study_ids duplicated in the meng data:
  # reported_test <- reported %>% 
  #   mutate(study_id = gsub("_et_al_|_and_.*_", "_", study_id)) # match the meng study_id structure
  # intersect(meng$study_id, reported_test$study_id)
  filter(study_id != "Craft_2007") %>% # 3 obs
  filter(study_id != "Noe_et_al_2016") %>% # 4 obs
  filter(study_id != "Boyd_and_Sommerfield_2016") %>%  # 12 obs; these core_ids were wrongly assigned, corrected below
  filter(study_id != "Callaway_et_al_2019") %>% # 34 obs; Callaway_2019 in reported values = Callaway_2012 in meng 

  # rename core_ids to match library data
  mutate(depth_min_cm = ifelse(!is.na(depth_max_cm) & is.na(depth_min_cm), 0, depth_min_cm),
         core_id = case_when(grepl("Nanticoke", site_id) ~ paste("NRM_", core_id, sep = ""), # Allen_et_al_2022
                             grepl("Typha", core_id) ~ gsub("Typha", "Typha ", core_id), # Arias-Ortiz_et_al_2021
                             grepl("Mouth", site_id) ~ "snipe_creek_M", # Arriola_and_Cable_2017
                             grepl("MC", site_id) ~ "snipe_creek_MC",
                             grepl("HE", site_id) ~ "snipe_creek_HE",
                             grepl("HI", site_id) ~ "snipe_creek_HI",
                             grepl("Upstream", site_id) ~ "snipe_creek_S",
                             grepl("Forest", site_id) ~ "snipe_creek_F",
                             grepl("CRMS|CRMS0", core_id) ~ gsub("CRMS|CRMS0", "", core_id), # Baustian_et_al_2021
                             study_id == "Boyd_2012" ~ gsub("-", "", core_id), # Boyd_2012
                             study_id == "Boyd_et_al_2017" ~ gsub("_", "", core_id), # Boyd_et_al_2017
                             grepl("Mudflat", core_id) ~ "ELM1812-MFA1", # Carlin_et_al_2021
                             grepl("Spartina", core_id) ~ "ELM1812-SPA2",
                             grepl("Pickleweed", core_id) ~ "ELM1812-PWA1",
                             grepl("Choptank_non_tidal", core_id) ~ "ChoptankNontidal", # Ensign_et_al_2020
                             grepl("Choptank_head_of_tide", core_id) ~ "ChoptankUpperTidal",
                             grepl("Choptank_TFFW", core_id) ~ "ChoptankMiddleTidal",
                             grepl("Choptank_marsh", core_id) ~ "ChoptankLowerTidal",
                             grepl("Pocomoke_non_tidal", core_id) ~ "PocomokeNontidal",
                             grepl("Pocomoke_head_of_tide", core_id) ~ "PocomokeUpperTidal",
                             grepl("Pocomoke_TFFW", core_id) ~ "PocomokeMiddleTidal",
                             grepl("Pocomoke_marsh", core_id) ~ "PocomokeLowerTidal",
                             grepl("SM02-VC1", core_id) ~ "SM02_VC1", # Johnson_et_al_2007
                             grepl("Traps_Bay", site_id) ~ paste("TB_", core_id, sep = ""), # McTigue_et_al_2020
                             grepl("2J", core_id) ~ "1_2J", # Poppe_and_Rybczyk_2018
                             grepl("3F", core_id) ~ "2_3F",
                             grepl("4F", core_id) ~ "3_4F",
                             study_id == "Poppe_and_Rybczyk_2018" & grepl("5B", core_id) ~ "4_5B",
                             grepl("5F", core_id) ~ "5_5F",
                             grepl("5J", core_id) ~ "6_5J",
                             grepl("St_Augustine_Salt_marsh", core_id) ~ "St_Augustine_Salt_Marsh", # Vaughn_et_al_2020
                             grepl("Waccasassa_Bay_Salt_marsh", core_id) ~ "Waccasassa_Bay_Salt_Marsh",
                             grepl("Luk_et_al_2020", study_id) ~ site_id, # Luk_et_al_2020
                             T ~ core_id),
         
         # rename site_ids to match library data where only site-level values are reported
         site_id = case_when(grepl("Wardlaw 1", site_id) ~ "Wardlaw_Shallow_Managed", # Drexler_et_al_2013
                             grepl("Hasty Point", site_id) ~ "Hasty_Point_Managed_Shallow",
                             grepl("Coastal EDU", site_id) ~ "Coastal_EDU_Field",
                             grepl("Bird Field Managed", site_id) ~ "Bird_Field",
                             grepl("Sandy Island Tidal", site_id) ~ "Sandy_Island_Natural",
                             grepl("Sandy Island Managed", site_id) ~ "Sandy_Island_Deeply_Flooded",
                             grepl("Poppe_and_Rybczyk_2019", study_id) ~ paste("Stillaguamish_", site_id, sep = ""), # Poppe_and_Rybczyk_2019
                             site_id == "QM" ~ "Snohomish_Quilceda_Marsh", # Poppe_et_al_2019
                             site_id == "HP" ~ "Snohomish_Heron_Point",
                             site_id == "OI" ~ "Snohomish_Otter_Island",
                             site_id == "NE" ~ "Snohomish_North_Ebey",
                             site_id == "SP" ~ "Snohomish_Spencer_Island",
                             site_id == "MA" ~ "Snohomish_Marysville_Mitigation",
                             site_id == "US" ~ "Snohomish_Union_Slough",
                             site_id == "SS" ~ "Snohomish_Smith_Island_City",
                             site_id == "QW" ~ "Snohomish_Qwuloolt",
                             site_id == "SN" ~ "Snohomish_Smith_Island_County",
                             site_id == "WW" ~ "Snohomish_WDFW_Wetland",
                             site_id == "WF" ~ "Snohomish_WDFW_Forest",
                             T ~ site_id),
         # note that Krauss_et_al_2018, Watson_and_Byrne_2013, and Abbott_et_al_2019 did not have Pb210 or Cs137 rates and were left out
         
         # correct this core_id assignment
         study_id = case_when(core_id %in% c("PM1", "PM2", "PM3", "PM4", "PM5", "PM6", 
                                             "PM7", "PM8", "PM9", "PM10", "PM11", "PM12") ~ "Boyd_and_Sommerfield_2016",
                              T ~ study_id)) 



## Pull depthseries summaries ####

## Core-level data
# Subset depthseries data into intervals represented by dating info
depth_weights_core <- reported_merge %>% 
  full_join(depth, by = c("study_id", "core_id")) %>% # site_id is inconsistently assigned, but core_id is corrected for in reported_merge

  # remove observations not matched in core-level reported values 
  filter(core_id %in% unique(reported_merge$core_id)) %>% 
  filter(study_id != "Luk_et_al_2020") %>% 

  # extract represented core length 
  group_by(study_id, core_id) %>% 
  summarise(rate_depth_max = max(depth_max_cm), # represented depth from reported values
            rate_depth_min = min(depth_min_cm),
            horizon_max = max(depth_max), # represented depth from depthseries table
            horizon_min = min(depth_min)) %>% 

  # determine the weight for dating rates (reported value depth interval / full depthseries depth interval,
  #                                        or: 1m depth interval / full depthseries depth interval,
  #                                        or: 1 [full dating rate])
  mutate(weight = case_when(!is.na(rate_depth_max) ~ (rate_depth_max - rate_depth_min) / (horizon_max - horizon_min), # if rate depth increment is present, use the interval to weight how much of the depthseries stock data to use
                            horizon_max > 100 ~ 100 / (horizon_max - horizon_min), # otherwise, if no rate depth increment is present and the core is deeper than 1m, weight the first meter of the core
                            T ~ 1)) # otherwise, use the entire rate given

# Extract DBD/FOM/FOC core summaries
depth_summary_core <- reported_merge %>% 
  full_join(depth, by = c("study_id", "core_id")) %>% 
  filter(core_id %in% unique(reported_merge$core_id)) %>% 
  filter(study_id != "Luk_et_al_2020") %>% 
  filter(!is.na(dry_bulk_density) | !is.na(fraction_carbon) | !is.na(fraction_organic_matter)) %>% 
  group_by(study_id, core_id) %>% 
  summarise(dry_bulk_density = mean(dry_bulk_density, na.rm = T), 
            fraction_organic_matter = mean(fraction_organic_matter, na.rm = T),
            fraction_carbon = mean(fraction_carbon, na.rm = T)) %>% 
  full_join(depth_weights_core)


## Site-level data
# Repeat the same as above for Luk_et_al_2020 site-level summaries
depth_weights_site <- reported_merge %>% 
  full_join(depth, by = c("study_id", "site_id")) %>% 
  filter(study_id == "Luk_et_al_2020") %>% 
  group_by(study_id, site_id) %>% 
  summarise(rate_depth_max = max(depth_max_cm), 
            rate_depth_min = min(depth_min_cm),
            horizon_max = max(depth_max), 
            horizon_min = min(depth_min)) %>% 
  mutate(weight = case_when(!is.na(rate_depth_max) ~ (rate_depth_max - rate_depth_min) / (horizon_max - horizon_min), 
                            horizon_max > 100 ~ 100 / (horizon_max - horizon_min), 
                            T ~ 1)) 

depth_summary_site <- reported_merge %>% 
  full_join(depth, by = c("study_id", "site_id")) %>% 
  filter(study_id == "Luk_et_al_2020") %>% 
  group_by(study_id, site_id) %>% 
  summarise(dry_bulk_density = mean(dry_bulk_density, na.rm = T), 
            fraction_organic_matter = mean(fraction_organic_matter, na.rm = T),
            fraction_carbon = mean(fraction_carbon, na.rm = T)) %>% 
  full_join(depth_weights_site) %>% 
  mutate(core_id = site_id) %>% 
  full_join(depth_summary_core) 
  



## Create the finalized table ####
synthdat <- reported_merge %>% 
  full_join(depth_summary_site, by = c("study_id", "core_id")) %>% 
  left_join(dated_us_cores) %>% 
  assignEcosystem() %>% 
  mutate(weighted_dry_bulk_density = dry_bulk_density * weight, 
         weighted_fraction_organic_matter = fraction_organic_matter * weight, 
         weighted_fraction_carbon = fraction_carbon * weight, 
         
         # convert FOM to FOC depending on habitat type (same conversion as above in reported)
         OM_density = weighted_dry_bulk_density * weighted_fraction_organic_matter, 
         foc_calc = ifelse(grepl("marsh", habitat) & weighted_fraction_organic_matter > 0, 0.486 * OM_density - 0.0160 * OM_density^2,
                           ifelse(grepl("swamp", habitat) & weighted_fraction_organic_matter > 0, 0.427 * OM_density + 0.0635 * OM_density^2,
                                  NA)),
         
         # use the FOC calculated above to convert accretion rate to C accumulation rate
         pb_CAR = ifelse(grepl("accretion", pb_accretion) & !is.na(weighted_fraction_carbon), pb210_rate * weighted_dry_bulk_density * weighted_fraction_carbon, # if FOC is present in the library # all present in pb_CAR
                         ifelse(grepl("accretion", pb_accretion) & !is.na(foc_calc), pb210_rate * foc_calc, # if FOC isn't present in the library, use calculated c_stock from FOM 
                                pb_CAR)), # the remainder are the already calculated/existing accumulation rates
         pb_CAR_se = ifelse(grepl("accretion", pb_accretion) & !is.na(weighted_fraction_carbon), pb210_rate_se * weighted_dry_bulk_density * weighted_fraction_carbon, 
                         ifelse(grepl("accretion", pb_accretion) & !is.na(foc_calc), pb210_rate_se * foc_calc, pb_CAR_se)), 
         cs_CAR = ifelse(grepl("accretion", cs_accretion) & !is.na(weighted_fraction_carbon), cs137_rate * weighted_dry_bulk_density * weighted_fraction_carbon,
                         ifelse(grepl("accretion", cs_accretion) & !is.na(foc_calc), cs137_rate * foc_calc, cs_CAR)),
         cs_CAR_se = ifelse(grepl("accretion", cs_accretion) & !is.na(weighted_fraction_carbon), cs137_rate_se * weighted_dry_bulk_density * weighted_fraction_carbon,
                            ifelse(grepl("accretion", cs_accretion) & !is.na(foc_calc), cs137_rate_se * foc_calc, cs_CAR_se))) %>% 
  filter(!is.na(pb_CAR|cs_CAR)) %>% # remove instances without accumulation rates
  distinct() %>% 
  select(study_id, site_id, core_id, core_count, pb_CAR, pb_CAR_se, cs_CAR, cs_CAR_se, climate_zone, ecosystem)

  


