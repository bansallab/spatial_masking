### This file generates estimates of the county-level bias in the FB responses
### by comparing true vaccination to FB reported vaccination with resampled
### vaccination observations using weights calculated via raking based on age, 
### sex, and education of respondents at county month level. 

# last run on 6/16/22 when decided on age, sex, and education, only using april and may
# last run on 5/19/22 when resampling with race variable, includes responses that don't answer mask, apr/may only
# last run on 5/16/22 when just averaging over april and may values
# last run on 5/4/22 with education as a raking variable
# last run on 4/20/22 fixing denominator for true vax data to pop 18+ instead of county pop of all ages
# last run on 4/15/22 adding resampling of vax observations
# last run on 3/17/22 fixing error that this was run on data through july 2021, should be june 2021

################

library(tidyverse)
library(lubridate)
library(vroom)
library(biglm)
library(lme4)

################
# acs over 18 population, produced in raking file
all_acs <- read_csv("data/acs_target_data_for_fb_raking.csv", col_types = "icidddddddddddddddddddddd") %>% 
  select(fips, NAME, all_18up)
library(tidycensus)
# INSERT CENSUS API KEY
total_pop <- get_acs(geography = "county", 
                     variables = c(all_ages_by_sex = "B01001_001"),
                     output = "wide",
                     year = 2018) %>% 
  mutate(fips = as.integer(GEOID))
pops <- all_acs %>% left_join(total_pop, by = "fips") %>% 
  rename(all_0up = all_ages_by_sexE) %>% 
  select(fips, all_0up, all_18up)

# true vax data
vax_path <- "data/archived_vaccination_data/data_master_county_20210823.csv"

vax <- read_csv(vax_path) %>% 
  mutate(week = floor_date(DATE, unit = 'week')) %>% 
  filter(CASE_TYPE == "Partial Coverage" | CASE_TYPE == 'Partial Coverage') %>% # includes completely vaxxed
  group_by(week, COUNTY, CASE_TYPE) %>% 
  summarize(p_vax_all = mean(CASES*.01)) %>% # cases is a percentage, we want proportion so divide by 100
  mutate(fips = as.numeric(COUNTY)) %>% 
  ungroup() %>% 
  select(week, fips, p_vax_all) %>% # p_vax_all is mean proportion vaccinated by county week with denominator of total population
  left_join(pops, by = "fips") %>% 
  mutate(p_vax_18up = p_vax_all * all_0up / all_18up)
# appears that new named counties have updated fips in this file, so maybe just had to adjust original fips county file

resampled_fb_vax <- read_csv("data/fb_resampled_vax_data_from_raking_age_sex_educ_FINAL_corrected_aprmay.csv",
                             col_types = "iDiidi")

clean_data_full <- left_join(filter(vax, week %within% interval(start = ymd("2021-04-01"), end = ymd("2021-05-31"))),
                             resampled_fb_vax) %>% 
  mutate(fips = as.factor(fips), # have to make factor for model
         diff = p_fb - p_vax_18up) %>% # get difference between observed p and true p, this is bad because it ignores the denominator
  filter(!is.na(n_vax), n_response > 0) %>% 
  filter(!is.na(fips), !is.na(diff)) %>% 
  droplevels() %>% # I think this is cause we might have gotten rid of a few fips so need to kind of refactor
  filter(0 < p_vax_18up, 1 > p_vax_18up) # drop county-weeks with no one or everyone vaccinated, will break logits

#saveRDS(clean_data_full, file = "data/fb_bias_clean_data_full_rake_aprmay.RDS")

#clean_data_full <- readRDS(file = "data/fb_bias_clean_data_full_rake_aprmay.RDS")
#########

# Try two stage approach -- generate fitted probabilities then OLS the difference
# between fitted and true to get bias correction

# Binomial GLM, p_vaccinated | county, week -- shrinkage of probability towards large sample sizes
fit.prob <- glmer(cbind(n_vax, n_response - n_vax) ~ 1 + (1|fips) + poly(week_rank, 2), 
                  family = binomial(link  = 'logit'), 
                  data = clean_data_full,
                  verbose = 2, # spit out more info during fitting process on how sampler is doing
                  nAGQ = 3) # kind of like tolerance of sampler, Gaussian quadriture method, higher is more precise
# laplace approximation: things that are exponential distributions converge to normal distributions with more observations
#saveRDS(fit.prob, file = "data/fb_bias_vax_ests_rake_aprmay.RDS")

# run the following instead of rerunning model
#fit.prob <- readRDS(file = "data/fb_bias_vax_ests_rake_aprmay.RDS")
clean.fitted <- clean_data_full %>% 
  mutate(fitted_p_fb = fitted(fit.prob), # extracts fitted values from model
         logit_vax = qlogis(p_vax_18up),
         logit_fitted_fb = qlogis(fitted_p_fb),
         #logit_diff = logit_vax - logit_fitted)
         logit_diff = logit_fitted_fb - logit_vax)

# fit a LMM to the difference in logits, random intercepts on county shrinking
# the bias correction (no sample size or SE correction here)
fit.bias.shrinkage <- lmer(logit_diff ~ 1+1|fips, data = clean.fitted)
#saveRDS(fit.bias.shrinkage, file = "data/fb_bias_vax_lmm_rake_aprmay.RDS")
# just for April and May
# apr_may <- clean.fitted %>% filter(week %within% interval(start = "2021-04-01", end = "2021-05-31"))
# fit.bias.shrinkage.aprmay <- lmer(logit_diff ~ 1 + 1|fips, data = apr_may)
# saveRDS(fit.bias.shrinkage.aprmay, file = "data/fb_bias_vax_lmm_aprmay_rake.RDS")

# run the following instead of rerunning model
#fit.bias.shrinkage <- readRDS(file = "data/fb_bias_vax_lmm_rake_aprmay.RDS")
# fit.bias.shrinkage.aprmay <- readRDS(file = "data/fb_bias_vax_lmm_aprmay_rake.RDS")

###########
# Save results 

df.bias <- clean.fitted %>% 
  select(p_vax_18up, p_fb, fips, week) %>% 
  mutate(bias_correction_logit = fitted(fit.bias.shrinkage)) %>% 
  select(fips, bias_correction_logit) %>% 
  unique()

df.bias <- df.bias %>% mutate(fips = as.numeric(levels(fips))[as.integer(fips)])

write_csv(df.bias, 'data/bias_correction_from_vaccination_08_23_rake_age_sex_educ_aprmay.csv')

###########
# Generate plots
clean.fitted <- clean.fitted %>% mutate(pred = fitted(fit.bias.shrinkage))

clean.fitted %>% 
  slice_sample(n=1000) %>% 
  filter(0 < p_vax_18up, 1 > p_vax_18up) %>% 
  ggplot(aes(p_vax_18up, plogis(qlogis(p_fb) - pred), color = log(n_response)))+ # + pred if logit_vax - logit_fitted
  geom_point()+
  geom_abline()+
  geom_smooth(se = F)+
  labs(x = 'True probability',
       y = 'Observed probabilty + bias correction',
       color = 'Sample size (log)')

clean.fitted %>% 
  slice_sample(n=1000) %>% 
  filter(0 < p_vax_18up, 1 > p_vax_18up) %>% 
  ggplot(aes(p_vax_18up, plogis(qlogis(fitted_p_fb) - pred), color = log(n_response)))+ # + pred if logit_vax - logit_fitted
  geom_point()+
  geom_abline()+
  geom_smooth(se = F)+
  labs(x = 'True probability',
       y = 'Fitted probabilty + bias correction',
       color = 'Sample size (log)')
