library(tidyverse)
library(brms)
library(cmdstanr)
library(bayesplot)
options(brms.backend = "cmdstanr")
library(tidybayes)
#library(tidylog)

### read in data
# from here: https://covid19.census.gov/datasets/USCensus::average-household-size-and-population-density-county/explore?location=4.945434%2C0.315550%2C1.99&showTable=true
pop_density <- read_csv("data/Average_Household_Size_and_Population_Density_-_County.csv",
                        col_types = "ddiddccddddddiddccccddd") %>% 
  select(GEOID, B01001_001E, B01001_001M, B01001_calc_PopDensity, B01001_calc_PopDensityM) %>% 
  rename(fips = GEOID, population = B01001_001E, population_moe = B01001_001M, 
         density = B01001_calc_PopDensity, density_moe = B01001_calc_PopDensityM)

onm <- read_csv("data/estimates/onm_processed.csv", col_types = "iDddiiiddd") %>% 
  left_join(pop_density, by = "fips") %>% 
  filter(! is.na(density)) %>% 
  filter(total_obs != 0) %>% 
  select(-c(mask_grocery_some_only_prop, mask_grocery_notso_only_prop, 
            mask_grocery_notatall_only_prop)) %>% 
  mutate(z_log10_density = (log10(density) - mean(log10(density)))/sd(log10(density)))

# pop_density NA removes fips codes from Wade Hampton, AK; Shannon, SD; Virgin Islands, American Samoa
# I fixed Wade Hampton, AK and Shannon, SD in mask_data.r, but not dealing with Virgin Islands or Samoa
# some county months have a total observation count of 0, this is because the county contained zip 
#   codes with observed successes, but those zip codes were primarily in a different county, so when
#   allocating successes, the county month of interest got < 1 successes which rounded to 0 so 
#   the county month is left in the df from before rounding the nonzero observations, but after 
#   rounding they have 0 observations
# we remove these county months from the model (there are 2682 of them)


# FORMULA
m1_formula <- bf(mask_grocery_very_suc | trials(total_obs) ~ 1 + z_log10_density,
                 family = binomial())

### prior
# set sample_prior = "only" so that data are not used
# think we can just test once because it doesn't use the data so should not vary across months
m1p <- brm(formula = m1_formula,
           data = onm %>% filter(month == "2020-10-01"),
           prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                     set_prior("normal(1, 1)", class = "b", coef = "z_log10_density")),
           sample_prior = "only",  # only for prior predictive dist
           warmup = 1000,
           iter = 2000,
           chains = 4,
           cores = 4)

# prior predictive checks
summary(m1p)
conditional_effects(m1p)
plot(m1p)
prior_summary(m1p)
pp_check(m1p, ndraws = 100) +  scale_x_log10()

# function to pull out predicted values
make_output_dfs <- function(model, onm_data){
  lin_predictions <- model %>% linpred_draws(newdata = onm_data) %>% 
    mutate(onm_est = plogis(.linpred)) %>% 
    mean_qi() %>% 
    ungroup() %>% 
    mutate(intercept_coef = fixef(model)[1],
           z_log10_dens_coef = fixef(model)[2],
           xbeta = intercept_coef + (z_log10_dens_coef * z_log10_density)) %>% 
    select(-c(.point, .interval, .row, population_moe, density_moe, 
              mask_grocery_some_prop, mask_grocery_some_suc))
  return(lin_predictions)
}

### posteriors
# formula is same across months, so only need to run different model for each month
months <- unique(onm$month)
onm_all_estimates <- data.frame()
for(i in 1:length(months)){
  print(months[i])
  # run model
  model <- brm(formula = m1_formula,
               data = (onm %>%  filter(month == months[i])),
               prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                         set_prior("normal(1, 1)", class = "b")),
               sample_prior = "no", 
               warmup = 1000, iter = 2000, chains = 4, cores = 4,
               file = paste0("data/fits/onm_", gsub("-", "_", substr(months[i], 3, 7)),
                             "fit_binomreg"),
               file_refit = "on_change")
  
  # doing posterior predictive checks
  summary(model)
  # y is observed responses, y_rep are predicted responses
  # saving posterior predictive check to make fig outside loop
  assign(paste0("p", i), pp_check(model, ndraws = 100) + scale_x_log10())
  
  model <- add_criterion(model, "loo")
  model_loo <- loo(model, moment_match = TRUE)
  model_loo # check that all values are not bad or very bad
  plot(model_loo)
  assign(paste0("m", gsub("-", "_", substr(months[i], 3, 7)), "_pareto_ks"), 
         model_loo$pointwise[,5])
  # sometimes want this to identify influential fips
  
  # save predictions
  if(i == 1){
    onm_all_estimates <- make_output_dfs(model, (onm %>% filter(month == months[i])))
  }else{
    onm_all_estimates <- onm_all_estimates %>% 
      bind_rows(make_output_dfs(model, (onm %>% filter(month == months[i]))))
  }
  
}

write_csv(onm_all_estimates, "data/estimates/onm_estimates_binomreg_FINAL.csv")

