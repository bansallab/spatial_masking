### Imports and cleans ONM anf FB data for later analysis
### Need for FB and ONM binomial regression models without raking
### Juliana Taube
### updated 07/19/2022 to remove NYT and YouGov

library(tidyverse)
library(readxl)
library(lubridate)
library(tidycensus)
library(vroom)
library(Hmisc)
# INSERT CENSUS API KEY 

### READ IN DATA
# general conversion data sets
fips_to_zips <- read_excel("data/COUNTY_ZIP_062021.xlsx")
zips_to_fips <- read_excel("data/ZIP_COUNTY_062021.xlsx")
# these files have Wade Hampton, AK 2270 fips and Shannon, SD 46113 fips so we will need to make the switch
#   to Kusilvak (2158) and Oglala Lakota (46102) after we do the crosswalk

# from here: https://gist.github.com/dantonnoriega/bf1acd2290e15b91e6710b6fd3be0a53
state_fips <- read_csv("data/state_fips.csv", col_types = "cic")

# urb_rur_codes <- read_excel("NCHSURCodes2013.xlsx") %>% rename(fips = `FIPS code`)

# from here: https://covid19.census.gov/datasets/USCensus::average-household-size-and-population-density-county/explore?location=4.945434%2C0.315550%2C1.99&showTable=true
# pop_density <- read.csv("Average_Household_Size_and_Population_Density_-_County.csv") %>% 
  # select(GEOID, B01001_001E, B01001_001M, B01001_calc_PopDensity, B01001_calc_PopDensityM) %>% 
  # rename(region = GEOID, population = B01001_001E, population_moe = B01001_001M, density = B01001_calc_PopDensity, density_moe = B01001_calc_PopDensityM)

### ONM ------------------------------------------------------------------------
onm_df <- read_csv("data/onm_sm_data_6_28_21.csv", col_types = "cDddciiiiiiiiii")

# clean onm zip codes and aggregate by month
onm <- onm_df %>% mutate(zip = (ifelse(nchar(zip_code) == 4, paste0("0", zip_code), zip_code))) %>% # adds 0 in front
  mutate(zip = (ifelse(nchar(zip) == 3, paste0("00", zip), zip))) # adds two zeros in front

# ONM DATA IS TOO SPARSE FOR WEEKLY ANALYSIS, SO LET'S STICK TO THE MONTHLY LEVEL
# add county codes for zipcodes, add week and month columns for aggregation
onm <- onm %>% mutate(date = ymd(response_date)) %>% 
  mutate(week = floor_date(date, unit = "week")) %>% 
  mutate(month = floor_date(date, unit= "month"))

# convert likert scale data to proportions but also keep track of counts
# likely_wear_mask_xxx <= 2 gives a vector of booleans true or false for each entry,
#   then taking the mean of this will give me the proportion TRUE which is those who somewhat or more likely to mask
# putting in separate data frames because we actually use different files to convert
# exclude response == 5 because they didn't answer question, focusing on grocery shopping as that is most public
# grocery shopping question had fewest non responses ~11,000 rows
onm_props <- onm %>% 
  filter(likely_wear_mask_grocery_shopping != 5) %>%  # get rid of no responses
  group_by(zip, month) %>% 
  summarize(mask_grocery_very = mean(likely_wear_mask_grocery_shopping <= 1),
            mask_grocery_some = mean(likely_wear_mask_grocery_shopping <= 2),
            mask_grocery_some_only = mean(likely_wear_mask_grocery_shopping == 2),
            mask_grocery_notso_only = mean(likely_wear_mask_grocery_shopping == 3),
            mask_grocery_notatall_only = mean(likely_wear_mask_grocery_shopping ==4),
            count = n())

onm_counts <- onm %>% 
  filter(likely_wear_mask_grocery_shopping != 5) %>%  # get rid of no responses
  group_by(zip, month) %>% 
  summarize(mask_grocery_very_n = sum(likely_wear_mask_grocery_shopping <= 1),
            mask_grocery_some_n = sum(likely_wear_mask_grocery_shopping <= 2),
            count = n())

# ok now we've got a single value per zip code per month, time to convert to the county level
# still need to deal with zips that have no associated fips..removing them for now
onm_props_fips <- left_join(onm_props, fips_to_zips, by = "zip")
onm_counts_fips <- left_join(onm_counts, zips_to_fips, by = "zip")

# zip codes without corresponding fips codes
# na_zips_resp <- onm %>% filter(likely_wear_mask_grocery_shopping != 5) %>%
#   left_join(fips_to_zips, by = "zip") %>%
#   filter(is.na(county))
valid_resp2 <- onm %>% filter(likely_wear_mask_grocery_shopping != 5) %>%
  left_join(fips_to_zips, by = "zip") %>%
  filter(! is.na(county)) %>% 
  filter(date %within% interval(start = "2020-09-01", end = "2021-05-31")) %>% 
  select(response_id) %>% distinct() %>% 
  nrow()
# na_zips <- onm_props_fips %>% filter(is.na(county)) %>% select(zip) %>% distinct()
# na_zips$zip

# now converting to fips codes by taking weighted average across a fips code
# using tot_ratio from HUD crosswalk file
onm_props_wt <- onm_props_fips %>% filter(! is.na(county)) %>% 
  mutate(mask_grocery_very_wt = mask_grocery_very * tot_ratio,
         mask_grocery_some_wt = mask_grocery_some * tot_ratio,
         mask_grocery_some_only_wt = mask_grocery_some_only * tot_ratio,
         mask_grocery_notso_only_wt = mask_grocery_notso_only * tot_ratio,
         mask_grocery_notatall_only_wt = mask_grocery_notatall_only * tot_ratio)

onm_counts_allocation_unr <- onm_counts_fips %>% filter(! is.na(county)) %>% 
  mutate(mask_grocery_very_amt = mask_grocery_very_n * tot_ratio,
         mask_grocery_some_amt = mask_grocery_some_n * tot_ratio,
         count_amt = count * tot_ratio)
# going to round later instead

# group by fips and sum up all the zip codes in a fips
onm_mo_by_county_props <- onm_props_wt %>% ungroup() %>% group_by(county, month) %>% 
  summarize(mask_grocery_very_prop = sum(mask_grocery_very_wt),
            mask_grocery_some_prop = sum(mask_grocery_some_wt),
            mask_grocery_some_only_prop = sum(mask_grocery_some_only_wt),
            mask_grocery_notso_only_prop = sum(mask_grocery_notso_only_wt),
            mask_grocery_notatall_only_prop = sum(mask_grocery_notatall_only_wt))

onm_mo_by_county_counts_unr <- onm_counts_allocation_unr %>% ungroup() %>% group_by(county, month) %>% 
  summarize(mask_grocery_very_suc = sum(mask_grocery_very_amt),
            mask_grocery_some_suc = sum(mask_grocery_some_amt),
            total_obs = sum(count_amt))

# proportion and num successes / total observations are not identical, I think the later we round the better, so doing that here
onm_mo_by_county_counts <- onm_mo_by_county_counts_unr %>%
  mutate(across(c(mask_grocery_very_suc, mask_grocery_some_suc,
                  total_obs), round))

onm_mo_by_county <- left_join(onm_mo_by_county_props, onm_mo_by_county_counts, by = c("month", "county"))

 # prep for choroplethr: make fips numeric, label columns
# let mask_grocery be the main focus for now as it is a public space
# remove NAs and NaNs for now but may need to come back and address
# FIX FIPS FOR WADE HAMPTON, AK (2270) AND SHANNON, SD (46113)
# SHOULD INSTEAD BE KUSILVAK, AK (2158) AND OGLALA LAKOTA, SD (46102)
onm_mo_by_county <- onm_mo_by_county %>% ungroup() %>% mutate(fips = as.numeric(county)) %>% 
  mutate(fips = ifelse(fips == 2270, 2158, ifelse(fips == 46113, 46102, fips))) %>% 
  select(fips, month, mask_grocery_very_prop, mask_grocery_some_prop, mask_grocery_very_suc, mask_grocery_some_suc, total_obs,
         mask_grocery_some_only_prop, mask_grocery_notso_only_prop, mask_grocery_notatall_only_prop)

### CLEAN FB DATA FOR BINOMIAL REGRESSION MODELs WITHOUT RAKING 
### FB ------------------------------------------------------------------------
path <- dir('INSERT PATH TO FACEBOOK DATA')[6:16]
# 2020-09 through 2021-07
# only read in masking column, date, sample weight, fips code
cols <- cols_only(C14 = col_integer(),
                  C14a = col_integer(),
                  StartDatetime = col_datetime(),
                  weight = col_double(),
                  fips = col_integer(), 
                  wave = col_integer())
# load FB data
fb_df <- paste0('INSERT PATH TO FACEBOOK DATA', path) %>%
  map(vroom, # maps function to a vector, vroom is faster version of read_csv
      delim = ',',
      col_types = cols) %>% # map passes back a list of each of the separate paths
  bind_rows() 

fb <- fb_df %>% 
  filter(!(is.na(C14) & is.na(C14a)), !is.na(fips)) %>% 
  mutate(date = ymd_hms(StartDatetime),
         week = ymd(floor_date(date, unit = "week")), # need extra mutate to get in Date format instead of "POSIXt"
         month = ymd(floor_date(date, unit= "month")),
         response = ifelse(is.na(C14), C14a, C14)) %>% 
  select(-StartDatetime) %>% 
  filter(response != 6) %>% # these people haven't been in public
  filter(fips < 60000) # not including outside 50 + DC
# group by month since that is level at which we can really compare with ONM
# not doing anything with weight right now
# want both proportions and number of successes

fb <- fb %>% 
  group_by(fips, month) %>% 
  summarise(mask_prop_all = mean(response == 1),
            mask_n_all = sum(response == 1), # num successes
            mask_prop_most = mean(response <= 2),
            mask_n_most = sum(response <= 2), # num successes
            mask_prop_some = mean(response <= 3),
            mask_n_some = sum(response <= 3), # num successes
            mask_prop_most_only = mean(response == 2),
            mask_prop_some_only = mean(response == 3),
            mask_prop_little_only = mean(response == 4),
            mask_prop_none_only = mean(response == 5),
            count = n())

# keep fips as numeric for choroplethr
fb <- fb %>% mutate(fips = ifelse(fips == 2270, 2158, ifelse(fips == 46113, 46102, fips)))


# save data frames of interest
write.csv(onm_mo_by_county,"data/onm_processed.csv", row.names = FALSE)
write.csv(fb,"data/fb_processed.csv", row.names = FALSE)
