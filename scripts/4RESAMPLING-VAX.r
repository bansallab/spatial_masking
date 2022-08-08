library(tidyverse)
library(lubridate)

# instead of reading in raw FB vaccination data, read in data with assigned weights from raking
fb_weighted <- read_csv("data/fb_raking_weights_age_sex_educ_FINAL_corrected.csv",
                        col_types = "iicidiiiDDiifffiffd") %>% 
  filter(!is.na(updated_weight))

fb_data_to_use <- fb_weighted %>% filter(month %within% interval(start = "2021-04-01", end = "2021-05-01")) %>% 
  filter(!is.na(vax)) %>% # assume V1 response MCAR w/in responses
  mutate(vax_bool = vax == 1) %>% 
  filter(updated_weight != 0) # remove observations of weight 0 since they will never be sampled anyway
# instead of just summarising here we need to do resampling to get n_vax and n_response
counties <- fb_data_to_use %>% pull(fips) %>% unique() %>% sort()
weeks <- fb_data_to_use %>% pull(week) %>% unique()
resampled_fb_vax <- data.frame(fips = c(), week = c(), n_vax = c(), n_response = c())
for(i in 1:length(counties)){
  print(counties[i])
  for(j in 1:length(weeks)){
    doi <- fb_data_to_use %>% filter(week == weeks[j] & fips == counties[i])
    sample <- slice_sample(doi, n = nrow(doi), weight_by = doi$updated_weight, replace = TRUE)
    suc <- sample %>% filter(vax_bool == 1) %>% nrow()
    resampled_fb_vax <- resampled_fb_vax %>% bind_rows(data.frame(fips = counties[i], week = weeks[j], 
                                                                  n_vax = suc, n_response = nrow(sample)))
  }
}

resampled_fb_vax <- resampled_fb_vax %>% 
  mutate(p_fb = n_vax/n_response,
         week_rank = as.integer(as.numeric(week)/7), # had to change this to 7
         week_rank = week_rank - min(week_rank)+1, # orders the weeks from 1 to wth week so not a date but an integer
         fips = as.numeric(fips)) %>% 
  filter(n_response > 0) %>% 
  mutate(fips = ifelse(fips == 2270, 2158, 
                       ifelse(fips == 46113, 46102,
                              fips))) 

write_csv(resampled_fb_vax, "data/fb_resampled_vax_data_from_raking_age_sex_educ_FINAL_corrected_aprmay.csv")
