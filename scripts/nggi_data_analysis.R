## CCN-NGGI-Analysis
## Authors: Jaxine Wolfe, Henry Betts, Rose Cheney, James Holmquist

## Workflow for standardizing C accumulation values for the 2022 NGGI
## And merging the table with the previous NGGI table from 2017

## Prepare Workspace ####

# load libraries
library(tidyverse)
library(readxl)

# source helper functions
source("scripts/nggi_utils.R")

# read in reported values table (new accumulation/accretion rate data since 2017)
reported_raw <- read_csv("data/NGGI_literature_reported_values.csv", na = c("---", "NA", ""))

# read in data from previous NGGI
mengdat <- read_csv("data/US-BC-Analysis-1-105.csv")

# read in climate zones (provided by Lisa Beers)
climate_zones <- read_xlsx("data/Climate zones.xlsx") %>%
  mutate(State = recode(State, "Delamare" = "Delaware")) %>%
  rename(admin_division = State, climate_zone = `Climate zone`) %>%
  select(-Color)

# pull CCN synthesis data
# Database V1.1.0
cores <- read_csv("https://ndownloader.figshare.com/files/43251309", guess_max = 11000)
ds <- read_csv("https://ndownloader.figshare.com/files/43251963", guess_max = 72000)

## Prepare the Synthesis Data ####

# Gapfill CCN depthseries data and calculate carbon stock
ds_stocks <- ds %>%
  left_join(cores %>% select(study_id, site_id, core_id, habitat, country)) %>%
  filter(country == "United States") %>% 
  # transfer representative depths to depth min and max column if they are NA
  mutate(depth_min = ifelse(is.na(depth_min) & !is.na(representative_depth_min), 
                            representative_depth_min, depth_min),
         depth_max = ifelse(is.na(depth_max) & !is.na(representative_depth_max), 
                            representative_depth_max, depth_max)) %>% 
  carbonStock(.)

# ggplot(ds_stocks %>% filter(carbon_density < 0.8), aes(carbon_density, habitat)) + geom_boxplot()

# Standardize calculations for specified interval(s)
standard_stocks <- standardizeDepths(ds_stocks)

## Standardize Reported Values ####

reported_aligned <- reported_raw %>% 
  
  # remove study_ids duplicated in the meng data:
  filter(study_id != "Craft_2007") %>% # 3 obs
  filter(study_id != "Noe_et_al_2016") %>% # 4 obs
  filter(study_id != "Johnson_et_al_2007") %>% # 1 obs
  filter(study_id != "Boyd_and_Sommerfield_2016") %>%  # 12 obs; this was duplicated by Boyd_2012
  # filter(study_id != "Weston_et_al_2023") %>% # 10 obs; data not in depthseries table
  filter(study_id != "Callaway_et_al_2019") %>% # 34 obs; Callaway_2019 in reported values = Callaway_2012 in meng
  # filter(!(study_id == "Boyd_2012" & site_id == "PM")) %>% 
  
  # rename core_ids to match library data
  mutate(source = "reported_value", # create a flag to remove unmatched cores after the depthseries data merge
         depth_min_cm = ifelse(!is.na(depth_max_cm) & is.na(depth_min_cm), 0, depth_min_cm),
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
                             # grepl("Traps_Bay", site_id) ~ paste("TB_", core_id, sep = ""), # McTigue_et_al_2020
                             grepl("2J", core_id) ~ "1_2J", # Poppe_and_Rybczyk_2018
                             grepl("3F", core_id) ~ "2_3F",
                             grepl("4F", core_id) ~ "3_4F",
                             study_id == "Poppe_and_Rybczyk_2018" & grepl("5B", core_id) ~ "4_5B",
                             grepl("5F", core_id) ~ "5_5F",
                             grepl("5J", core_id) ~ "6_5J",
                             grepl("St_Augustine_Salt_marsh", core_id) ~ "St_Augustine_Salt_Marsh", # Vaughn_et_al_2020
                             grepl("Waccasassa_Bay_Salt_marsh", core_id) ~ "Waccasassa_Bay_Salt_Marsh",
                             study_id == "Weston_et_al_2023" ~ paste0(site_id, "_", core_id),
                             core_id == "2L-Pb" ~ "TB_2L",
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
         
         # correct this core_id assignment
         study_id = case_when(
           # core_id %in% c("PM1", "PM2", "PM3", "PM4", "PM5", "PM6", 
           #                                   "PM7", "PM8", "PM9", "PM10", "PM11", "PM12") ~ "Boyd_and_Sommerfield_2016",
                              # study_id == "Weston_et_al_2023" ~ "Weston_et_al_2020",
                              T ~ study_id),
         
         pb210_rate_unit = ifelse(study_id %in% c("Arias-Ortiz_et_al_2021", "Carlin_et_al_2021"), 
                                  "gramsPerSquareMeterPerYear", pb210_rate_unit)) 

# check mismatches in core IDs
mismatch <- anti_join(reported_aligned %>% select(study_id, core_id), cores) %>% drop_na(core_id)

## Handle outliers in the reported values table

# 2L-Pb core needs to be averaged across its depth intervals
core_2L_Pb <- reported_aligned %>% 
  # drop_na(core_id) %>% add_count(site_id, core_id) %>% filter(n > 1) %>% 
  filter(core_id == "TB_2L") %>% 
  group_by(study_id, site_id, core_id, habitat, management, accumulation_type, pb210_rate_unit) %>%
  summarise(pb210_rate_se = se(as.numeric(pb210_rate)),
            pb210_rate = mean(as.numeric(pb210_rate)))

# Drexler_et_al_2013 has pb210_rate as a range - remove and replace with summarized data
reported_drex <- reported_aligned %>% 
  filter(grepl("Drexler_et_al_2013", study_id)) %>% 
  separate(pb210_rate, into = c("pb210_rate_min", "pb210_rate_max"), sep = "-", convert = T) %>% 
  separate(depth_max_cm, into = c("depth_min_cm", "depth_max_cm"), sep = "-", convert = T) %>% 
  mutate(pb210_rate = (as.numeric(pb210_rate_min) + as.numeric(pb210_rate_max))/2) %>% 
  select(-pb210_rate_min, -pb210_rate_max) %>% 
  select_if(~!all(is.na(.)))

# Join climate zone data with core table
core_categories <- cores %>%
  filter(country == "United States") %>% 
  left_join(climate_zones) %>% 
  rename(habitat_ccn = habitat) %>% # so it won't be confused with the reported values table one
  select(study_id, core_id, habitat_ccn, admin_division, climate_zone)

# core_categories %>% filter(study_id %in% unique(mismatch$study_id)) %>% 
  # distinct(study_id, climate_zone) %>% arrange(climate_zone)

# Create the finalized reported 
reported_clean <- reported_aligned %>% 
  # Drexler 2013 is going to replaced by reported_drex, and Gerlach is one marker rate from a disturbed site so we'll leave it out
  filter(!(study_id %in% c("Drexler_et_al_2013", "Gerlach_et_al_2017"))) %>%
  filter(is.na(core_id) | core_id != "TB_2L") %>% 
  # carbon stock throws NA coercion (we're not using this column in this analysis however)
  mutate(across(c(depth_min_cm, depth_max_cm, pb210_rate, pb210_rate_se, marker_rate, carbon_stock), as.numeric)) %>% 
  bind_rows(reported_drex) %>% 
  bind_rows(core_2L_Pb) %>% 
  left_join(core_categories) %>% 
  mutate(inventory_year = "nggi_2022",
         accumulation_type = ifelse(study_id == "Kemp_et_al_2020", "sediment accretion", accumulation_type)) %>% 
  
  # patch in some climate zones for studies that only had site-level rates
  mutate(climate_zone = case_when(study_id %in% c("Poppe_et_al_2019", "Poppe_and_Rybczyk_2019", "Thom_1992",
                                                  "Peck_et_al_2020", "McTigue_et_al_2020", "Krauss_et_al_2018",
                                                  "Drexler_et_al_2019", "Drexler_et_al_2013", "Boyd_2012", "Boyd_and_Sommerfield_2016") ~ "Warm Temperate",
                                  study_id %in% c("Breithaupt_et_al_2014", "Abbott_et_al_2019", "Piazza_et_al_2021") ~ "Subtropical",
                                  study_id %in% c("Drexler_et_al_2009", "Watson_and_Byrne_2013") ~ "Mediterranean",
                                  study_id == "Luk_et_al_2020" ~ "Cold Temperate",
                                  T ~ climate_zone)) %>% 
  
  # small habitat assignment fixes from the latest database update
  mutate(habitat = case_when(core_id == "ELM1812-MFA1" ~ "mudflat",
                             core_id %in% c("BACHI", "FW", "TT") ~ "swamp",
                             core_id %in% c("St_Augustine_Mangrove", "Waccasassa_Bay_Mangrove") ~ "scrub/shrub",
                             T ~ habitat)) %>% 
  
  # assign ecosystem based on habitat
  mutate(ecosystem = case_when(habitat == "marsh" ~ "Estuarine Emergent Wetland",
                               grepl("scrub|shrub", habitat) ~ "Estuarine Scrub/Shrub Wetland",
                               habitat == "swamp" ~ "Palustrine Forested Wetland",
                               habitat == "mangrove" ~ "Estuarine Forested Wetland",
                               habitat == "seagrass" ~ "Estuarine Aquatic Bed",
                               T ~ NA_character_)) %>% 
  
  mutate(management = case_when(study_id %in% c("Giblin_and_Forbrich_2018", "Weston_et_al_2023") ~ "natural",
                                T ~ management)) %>% 
  
  # when CRS and CIC Pb210 rates are available, replace with their average
  mutate(pb210_rate = ifelse(!is.na(pb210_rate_cic & pb210_rate), (pb210_rate_cic + pb210_rate)/2,
                             ifelse(!is.na(pb210_rate_cic), pb210_rate_cic, pb210_rate)),
         pb210_rate_se = ifelse(!is.na(pb210_rate_cic_se & pb210_rate_se), (pb210_rate_cic_se + pb210_rate_se)/2,
                                ifelse(!is.na(pb210_rate_cic_se), pb210_rate_cic_se, pb210_rate_se)),
         pb210_rate_unit = ifelse(is.na(pb210_rate), NA, pb210_rate_unit),
         cs137_rate_unit = ifelse(is.na(cs137_rate), NA, cs137_rate_unit)
         ) %>% 
  
  # finally, join the stocks data derived from the CCN depthseries
  left_join(standard_stocks %>% select(-site_id)) %>% 
  arrange(study_id, site_id, core_id) 

## Convert accretion rates to CAR
reported_standardized <- accretionToCAR(reported_clean)

# studies with sediment accretion
# 1 Allen_et_al_2022         
# 2 Boyd_2012                
# 3 Boyd_and_Sommerfield_2016 # removed
# 4 Boyd_et_al_2017          
# 5 Drexler_et_al_2009  # radiocarbon rate
# 6 Gerlach_et_al_2017  # not included 
# 7 Lagomasino_et_al_2020    
# 8 Luk_et_al_2020           
# 9 Smith_et_al_2015         
# 10 Thom_1992                
# 11 Weis_et_al_2001          
# 12 Vaughn_et_al_2020 

# Baustian is missing rate data for some cores?


## Merge all Reported Data ####

# Prepare Meng's 2017 synthesis data for merge with reported values
meng <- mengdat %>% 
  rename(latitude = Latitude,
         longitude = Longitude,
         habitat = Ecosystem,
         # ecosystem = CCAP_Class,
         management = Management,
         climate_zone = Climate_Zone,
         # carbon_stock = SOC1, 
         # carbon_stock_unit = SOC1units,
         depth_max_cm = SOC1depth, 
         salinity_class = Salinity_Regime, 
         pb210_CAR = delSOC1Pb, 
         pb210_CAR_unit = delSOC1units,
         cs137_CAR = delSOC2Cs,
         cs137_CAR_unit = delSOC2units,
         marker_rate = delSOC3Marker,
         marker_rate_unit = delSOC3units,
         SET_rate = delSOC4SET,
         SET_rate_unit = delSOC4units,
         radiocarbon_rate = delSOC5RadioC,
         radiocarbon_rate_unit = delSOC5units) %>% 
  mutate(inventory_year = "nggi_2017",
         study_id = paste(First_Author, Year, sep = "_"),
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
         
         # From Meng workflow
         # if Ecosystem='marsh' then Ecosystem='Estuarine_Emergent_Wetlands';
         # if Ecosystem='mangrove' and Stature='shrub' then Ecosystem='Estuarine_Emergent_Wetlands';
         # if Ecosystem='mangrove' and Stature ne 'shrub' then Ecosystem='Estuarine_Forested_Wetlands';
         # if Ecosystem='tidal_fresh_marsh' then Ecosystem='Palustrine_Emergent_Wetland';
         # if Ecosystem='tidal_fresh_forest' then Ecosystem='Palustrine_Forested_Wetland';
         ecosystem = case_when(habitat == "marsh" ~ "Estuarine Emergent Wetland",
                               habitat == "mangrove" & Stature == "shrub" ~ "Estuarine Emergent Wetland",
                               habitat == "mangrove" & Stature == "tree" ~ "Estuarine Forested Wetland",
                               habitat == "mangrove" & is.na(Stature) ~ "Estuarine Forested Wetland",
                               habitat == "tidal_fresh_marsh" ~ "Palustrine Emergent Wetland",
                               habitat == "tidal_fresh_forest" ~ "Palustrine Forested Wetland",
                               TRUE ~ NA), 
         climate_zone = recode(climate_zone,
                         "temperate_cold" = "Cold Temperate",
                         "temperate_warm" = "Warm Temperate"),
         pb210_CAR_unit = recode(pb210_CAR_unit, "gC_m2" = "gramsPerSquareMeterPerYear"),
         cs137_CAR_unit = recode(cs137_CAR_unit, "gC_m2" = "gramsPerSquareMeterPerYear"),
         carbon_stock  = case_when(SOC1units == "OCg_cc" ~ SOC1 * 10000, 
                                   SOC1units == "gC_m2" ~ SOC1 / 100, 
                                   SOC1units == "MgC_ha" ~ SOC1,
                                   T ~ NA),
         carbon_stock_unit = ifelse(!is.na(carbon_stock), "megagramsPerHectare", NA),
         depth_min_cm = ifelse(!is.na(depth_max_cm), 0, NA)) %>% 
  
  # recode management
  mutate(management = recode(management, 
                             "N" = "natural",
                             "R" = "restored",
                             "D" = "disturbed",
                             "M" = "impounded")) %>% 
  select_if(~!all(is.na(.)))


# bind all 
reported_all <- bind_rows(reported_standardized, meng) %>% 
  mutate(climate_zone = str_to_title(climate_zone)) %>% 
  select(names(reported_standardized), everything()) %>% 
  # remove columns which are unneccessary for calculations moving forwards
  select(-c(pb210_rate, pb210_rate_se, pb210_rate_unit, pb210_rate_cic, pb210_rate_cic_se,
            cs137_rate, cs137_rate_se, cs137_rate_unit, dry_bulk_density, fraction_organic_matter,
            fraction_carbon, reported_rates, notes, publication_table, contains("stock"), SOC1, SOC1units,
            Litter, BD, contains("Species"), contains("AGB"), contains("BGB"))) %>% 
  filter_at(vars(pb210_CAR, pb210_CAR_unit, cs137_CAR, cs137_CAR_unit), any_vars(!is.na(.))) %>% 
  select_if(~!all(is.na(.))) %>% 
  arrange(study_id)

# issues: 
# Boyd 2012, core NCGB1
# Lagomasino 2020, HG3, PP1, PP3

# test <- reported_all %>%
#   filter(inventory_year == "2022") %>%
#   filter(management == "natural") %>%
#   filter(habitat != "mudflat") %>%
#   mutate(core_count = ifelse(is.na(core_count), 1, core_count))
# sum(test$core_count)
# 
# nrow(meng %>% 
#   filter_at(vars(pb210_CAR, pb210_CAR_unit, cs137_CAR, cs137_CAR_unit), any_vars(!is.na(.))) %>% 
#   filter(management == "natural"))

## NGGI Summary Calculations ####

nggi_smry <- reported_all %>% 
  filter(management == "natural") %>% 
  filter(habitat != "mudflat") %>% 
  # filter(is.na(core_count) | core_count == 1) %>% # until we know how to handle a mean from multiple cores
  # convert from gC m-2 yr-1 to gC ha-1 yr-1 and calculate the geometric mean
  mutate(NewdelSOC1Pb = ifelse(pb210_CAR_unit != "centimeterPerYear", log10(pb210_CAR/100), NA),
         NewdelSOC2Cs = ifelse(cs137_CAR_unit != "centimeterPerYear", log10(cs137_CAR/100), NA)) %>% 
  group_by(ecosystem, climate_zone) %>% 
  summarize(MNewdelSOC1Pb = mean(NewdelSOC1Pb, na.rm = T),
            SeNewdelSOC1Pb = se(NewdelSOC1Pb, na.rm = T),
            SdNewdelSOC1Pb = sd(NewdelSOC1Pb, na.rm = T),
            MNewdelSOC2Cs = mean(NewdelSOC2Cs, na.rm = T),
            SeNewdelSOC2Cs = se(NewdelSOC2Cs, na.rm = T),
            SdNewdelSOC2Cs = sd(NewdelSOC2Cs, na.rm = T),
            pb210_n = sum(!is.na(NewdelSOC1Pb)),
            cs137_n = sum(!is.na(NewdelSOC2Cs))) %>% 
  ungroup() %>% 
  mutate(MGeodelSOC1Pb = 10^MNewdelSOC1Pb,
         MGeodelSOC2Cs = 10^MNewdelSOC2Cs,
         # confidence intervals
         LowerCIGeodelSOC1Pb = 10^(MNewdelSOC1Pb - (SeNewdelSOC1Pb*1.96)),
         LowerCIGeodelSOC2Cs = 10^(MNewdelSOC2Cs - (SeNewdelSOC2Cs*1.96)),
         UpperCIGeodelSOC1Pb = 10^(MNewdelSOC1Pb + (SeNewdelSOC1Pb*1.96)),
         UpperCIGeodelSOC2Cs = 10^(MNewdelSOC2Cs + (SeNewdelSOC2Cs*1.96)))


## Data Visualization ####

rmarkdown::render(input = "scripts/NGGI_2022_report.Rmd",
                  output_dir = "report")


## Write final NGGI table ####

nggi_smry_clean <- nggi_smry %>% 
  select(ecosystem, climate_zone, contains("Geo"), contains("_n")) %>% 
  select(ecosystem, climate_zone, contains("pb"), contains("cs")) %>% 
  mutate(across(everything(), ~replace_na(.x, NA)))

write_csv(nggi_smry_clean, "report/NGGI_2022_CAR.csv")

## Bibliography ####

bib <- read_csv("https://ndownloader.figshare.com/files/43251330")

library(RefManageR)

studies_for_bib <- reported_all %>% 
  filter(management == "natural") %>% 
  filter(habitat != "mudflat") %>% 
  distinct(study_id) %>% pull(study_id)

# Piazza_et_al_2020: https://pubs.usgs.gov/of/2011/1094/OF11-1094.pdf
# all disturbed, don't include (hurricane though, which is still natural..)
# piazza_article <- as.data.frame(GetBibEntryWithDOI("10.5066/P9D8WTQW")) %>%
#   mutate(study_id = "Piazza_et_al_2021",
#          bibliography_id = "Piazza_et_al_2020_data",
#          publication_type = "primary dataset") %>%
#   remove_rownames()

nggi_bibs <- bib %>% 
  filter(study_id %in% studies_for_bib) %>% 
  mutate(across(everything(), as.character)) %>% 
  arrange(study_id)
# filter(bibtype != "Misc")

write_csv(nggi_bibs, "report/NGGI_2022_CAR_bibliography.csv")

## Misc code ----

# studies with split core methods
# Okeefe-Suttles_et_al_2021_Cape
# Luk_et_al_2020
# Breithaupt_et_al_2014
# McTigue_et_al_2020

# No associated article: 
# Messerschmidt_and_Kirwan_2020 - CCN release, has interval sedimentation rate (mm yr-1). Need to derive %OC and calculate CAR
# O'keefe Suttles, other release - CAR needs to be extracted from data releases
# Weston_et_al_2020, CCN release - no accumulation rates, we have to calculate these or ask Nat
# Buffington_et_al_2020 - CCN release, 1 core annual accretion rate given as 4.8mm/yr, Need to derive %OC since %C is modeled total?
# Gonneea_et_al_2018 - other release, depth interval CAR in original 
# Breithaupt_et_al_2020 - CCN release, maybe a rate column? might have to age depth model it

# Study notes:
# Rodriguez_et_al_2022: https://doi.org/10.1038/s43247-022-00501-x - include next time


