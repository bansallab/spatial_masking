library(tidyverse)
library(choroplethr)
library(choroplethrMaps)
library(lubridate)
library(ggpubr)
library(viridis)
library(RColorBrewer)

# before can impute these missing counties we need to know which are missing that matter, right?
nbrs <- read.csv("data/county_neighbors_fips.txt", header = FALSE) %>% rename(fips1 = V1, fips2 = V2)
bias <- read.csv("data/bias_correction_from_vaccination_08_23_rake_age_sex_educ_aprmay.csv", header = TRUE)
fb_fips <- read_csv("data/fb_resampled_data_from_raking_age_sex_educ_FINAL_corrected.csv", col_types = "iDiid") %>%
  select(fips) %>% distinct()

overlap_fips <- fb_fips %>% left_join(bias, by = "fips") %>% filter(is.na(bias_correction_logit))
# 45 counties to impute with change to april and may only (make sense, less data? yes some missing any obs in apr/may)

# time to automate
additions <- data.frame(fips = c(), bias_correction_logit = c())
counties <- overlap_fips %>% pull(fips) %>% sort()
for(i in 1:length(counties)){
  county_nbrs <- nbrs %>% filter(fips2 == counties[i])
  nbr_bias <- county_nbrs %>% rename(fips = fips1) %>% left_join(bias, by = "fips")
  county_bias <- sum(nbr_bias$bias_correction_logit, na.rm = TRUE)/nrow(nbr_bias)
  additions <- additions %>% bind_rows(data.frame(fips = counties[i], bias_correction_logit = county_bias))
}

bias <- bias %>% bind_rows(additions)

write.csv(bias, 
          "data/bias_correction_from_vaccination_08_23_complete_rake_age_sex_educ_aprmay_FINAL.csv", 
          row.names = FALSE)
