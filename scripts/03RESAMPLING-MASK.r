library(tidyverse)
library(lubridate)
# this needs to weight by updated_weight
fb_weighted2 <- read_csv("data/fb_raking_weights_age_sex_educ_FINAL_corrected.csv", 
                        col_types = "iicidiiiDDiifffiffd") %>% 
  filter(!is.na(self_response), self_response != 6) %>% 
  mutate(mask_most = ifelse(self_response <= 2, 1, 0)) %>% 
  filter(!is.na(updated_weight))
counties <- fb_weighted %>% pull(fips) %>% unique() %>% sort()
months <- fb_weighted %>% pull(month) %>% unique()
resampled_df <- data.frame(fips = c(), month = c(), successes = c(), trials = c())
for(i in 1:length(counties)){
  print(counties[i])
  for(j in 1:length(months)){
    doi <- fb_weighted %>% filter(month == months[j] & fips == counties[i])
    sample <- slice_sample(doi, n = nrow(doi), weight_by = doi$updated_weight, replace = TRUE)
    suc <- sample %>% filter(mask_most == 1) %>% nrow()
    resampled_df <- resampled_df %>% bind_rows(data.frame(fips = counties[i], month = months[j], 
                                                          successes = suc, trials = nrow(sample)))
  }
}
resampled_df <- resampled_df %>% mutate(mask_prop_most = successes/trials) %>% 
  filter(trials > 0)
write_csv(resampled_df, "data/fb_resampled_data_from_raking_age_sex_educ_FINAL_corrected.csv")
