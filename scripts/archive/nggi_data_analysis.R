## Workflow for processing the final table

# Overview

# Meng Lu Workflow
# 1.3 Calculate log10 value for each variable
# 1.4 Calculate the mean values for log10 values
# 1.5 Calculate the geomean values for each variable (with confidence interval)

library(tidyverse)

## Read in final table

nggi_table <- read_csv("coastal_wetland_nggi/data/intermediate/reported_values_clean.csv")

## NGGI Summary Calculations

nggi_smry <- nggi_table %>% 
  # calculate log10 of the activities in order to determine the geomean
  mutate(log10_pb210 = log10(pb210_standardized),
         log10_pb210_se = log10(pb210_standardized_se),
         log10_cs137 = log10(cs137_standardized),
         log10_cs137_se = log10(cs137_standardized_se))

# accumulation %>% 
#   drop_na(pb210_standardized) %>% 
#   ggplot(aes(pb210_standardized, col = accumulation_type)) +
#   geom_density() +
#   facet_wrap(~accumulation_type, scales = "free", dir = "v")
# 
# accumulation %>%
#   drop_na(pb210_standardized, cs137_standardized) %>%
#   ggplot(aes(cs137_standardized, pb210_standardized, col = accumulation_type)) +
#   geom_point()

# calculate the geometric mean for each accumulation rate: mean = (a1 x a2 x...x an)^(1/n)
geomean_pb210 <- nggi_smry %>% 
  drop_na(log10_pb210) %>% 
  group_by(ecosystem, climate_zone) %>%
  summarise(am = mean(log10_pb210, na.rm = T),
            stdev = sd(log10_pb210, na.rm = T),
            n = n(),
            ci_upper = am + (1.96*stdev)/sqrt(n),
            ci_lower = am - (1.96*stdev)/sqrt(n)
            # gm_test = prod(pb210_standardized)^(1/n)
  ) %>% 
  mutate(gm = 10^am,
         gm_ci_upper = 10^ci_upper,
         gm_ci_lower = 10^ci_lower
  )

ggplot(geomean_pb210, aes(ecosystem, gm, col = climate_zone)) + 
  geom_point() +
  coord_flip()

geomean_cs137 <- standardized %>% 
  drop_na(log10_cs137) %>% 
  group_by(ecosystem, climate_zone) %>%
  summarise(am = mean(log10_cs137, na.rm = T),
            stdev = sd(log10_cs137, na.rm = T),
            n = n(),
            ci_upper = am + (1.96*stdev)/sqrt(n),
            ci_lower = am - (1.96*stdev)/sqrt(n)) %>% 
  mutate(gm = 10^am,
         gm_ci_upper = 10^ci_upper,
         gm_ci_lower = 10^ci_lower
  )

# things to iron out
# total vs. organic carbon sequestration
# propagation of error from values that are already averaged to the site-level
# standardize to depth?
# 
# ## GapFilling Depthseries - IGNORE FOR NOW ####
# 
# gapfill_ds <- depthseries %>% 
#   # filter(study_id %in% unique(dated_us_cores$study_id))
#   filter(core_id %in% unique(dated_us_cores$core_id)) %>% 
#   left_join(dated_us_cores %>% distinct(core_id, habitat)) %>% # add habitat info from cores table
#   left_join(methods %>% distinct(study_id, method_id, fraction_carbon_type)) %>%   
#   select_if(function(x) {!all(is.na(x))}) %>% 
#   # flag rows by the presence/absence of things 
#   mutate(action_flag = case_when(!is.na(fraction_carbon) & fraction_carbon_type == "organic carbon" ~ "no action",
#                                  !is.na(fraction_carbon) & fraction_carbon_type == "total carbon" ~ "convert TC to OC",
#                                  is.na(fraction_carbon) & !is.na(fraction_organic_matter) ~ "convert LOI to C",
#                                  T ~ NA_character_)) %>% 
#   # gapfill fraction carbon 
#   # relationship to convert total carbon to organic?
#   # start with Craft 1991 relationship, increase complexity from there
#   mutate(soc = case_when(action_flag == "no action" ~ fraction_carbon * 100,
#                          action_flag == "convert LOI to C" & habitat == "marsh" ~ 0.4 * (100*fraction_organic_matter) + 0.0025*(100*fraction_organic_matter)^2,
#                          action_flag == "convert LOI to C" & habitat == "mangrove" ~ 0.415 * (100*fraction_organic_matter) + 2.89,
#                          action_flag == "convert TC to OC" ~ NA_real_,
#                          T ~ NA_real_)
#          # measured_or_modeled = case_when(study_id %in% c("Buffington_et_al_2020", "Drexler_et_al_2009",
#          #                                                     "Keshta_et_al_2020", "Radabaugh_et_al_2018",
#          #                                                     "Rodriguez_et_al_2022") ~ "modeled",
#          #                                 action_flag == "convert LOI to C" ~ "modeled",
#          #                                     T ~ "measured"),
#   ) %>%
#   drop_na(action_flag) %>% 
#   select(study_id, fraction_organic_matter, contains("carbon"), action_flag, soc, everything())
# 
# # prepare reported values table and merge with core table
# # do this after?
# # join_tables <- reported %>%
# #   # select(-contains("stock")) %>% 
# #   filter(management == "natural") %>% 
# #   left_join(cores) %>% 
# #   left_join(methods) %>% 
# #   left_join(impacts) %>% 
# #   select_if(function(x) {!all(is.na(x))}) 
# 
# # total vs organic carbon
# gapfill_ds %>%  
#   drop_na(fraction_carbon) %>% 
#   ggplot(aes(fraction_organic_matter*100, soc, col = fraction_carbon_type)) +
#   geom_point(alpha = 0.5, pch = 1) 
# 
# gapfill_ds %>%  
#   # filter(study_id == "Drexler_et_al_2009") %>% 
#   # filter(habitat == "swamp") %>% 
#   drop_na(soc) %>% 
#   ggplot(aes(fraction_organic_matter * 100, soc, col = action_flag)) +
#   geom_point(alpha = 0.5, pch = 1) +
#   facet_wrap(~action_flag)
# 
# gapfill_ds %>% drop_na(soc) %>% 
#   ggplot(aes(soc)) +
#   geom_density() + geom_rug()
# 
# 
