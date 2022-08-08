library(tidyverse)
library(brms)
library(cmdstanr)
library(bayesplot)
library(insight)
library(tidybayes)
#library(tidylog)
library(ggpubr)
options(brms.backend = "cmdstanr")

### read in data
# from here: https://covid19.census.gov/datasets/USCensus::average-household-size-and-population-density-county/explore?location=4.945434%2C0.315550%2C1.99&showTable=true
pop_density <- read_csv("data/Average_Household_Size_and_Population_Density_-_County.csv",
                        col_types = "ddiddccddddddiddccccddd") %>% 
  select(GEOID, B01001_001E, B01001_calc_PopDensity) %>% 
  rename(fips = GEOID, population = B01001_001E, density = B01001_calc_PopDensity)

resampled_df <- read_csv("data/fb_resampled_data_from_raking_age_sex_educ_FINAL_corrected.csv",
                         col_types = "iDiid")

to_model <- resampled_df %>% 
  left_join(pop_density, by = "fips") %>% 
  ungroup() %>% 
  mutate(z_log10_density = (log10(density) - mean(log10(density)))/sd(log10(density)))


# FORMULA
m1_formula <- bf(successes | trials(trials) ~ 1 + z_log10_density,
                 family = binomial())

### prior
# set sample_prior = "only" so that data are not used
# think we can just test once because it doesn't use the data so should not vary across months
m1p <- brm(formula = m1_formula,
           data = to_model %>% filter(month == "2020-10-01"),
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
make_output_dfs <- function(model, fb_data){
  lin_predictions <- model %>% linpred_draws(newdata = fb_data) %>% 
    mutate(p_est = plogis(.linpred)) %>% 
    mean_qi() %>% 
    ungroup() %>% 
    mutate(intercept_coef = fixef(model)[1],
           z_log10_dens_coef = fixef(model)[2],
           xbeta = intercept_coef + (z_log10_dens_coef * z_log10_density)) %>% 
    select(-c(.point, .interval, .row))
  
  return(lin_predictions)
}

### posteriors
### RAKING WITH POOLING NO BIAS  -----------------------------------------------------------------------
months <- unique(to_model$month)
fb_all_estimates <- data.frame()
for(i in 1:length(months)){
  print(months[i])
  # run model
  model <- brm(formula = m1_formula,
               data = (to_model %>% filter(month == months[i])),
               prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                         set_prior("normal(1, 1)", class = "b")),
               sample_prior = "no", 
               warmup = 1000, iter = 2000, chains = 4, cores = 4,
               file = paste0("data/fits/fb_", gsub("-", "_", substr(months[i], 3, 7)),
                             "fit_binomreg_rake"),
               file_refit = "on_change")
  
  # doing posterior predictive checks
  print(summary(model))
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
    fb_all_estimates <- make_output_dfs(model, (to_model %>% filter(month == months[i])))
  }else{
    fb_all_estimates <- fb_all_estimates %>% 
      bind_rows(make_output_dfs(model, (to_model %>% filter(month == months[i]))))
  }
}

pdf("figures/for-paper/supplement/posterior-predictive-checks-binomreg-rake.pdf", height = 8, width = 12)
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3, nrow = 3, 
          labels = c("  A", "  B", "  C", "  D", "  E", "  F", "  G", "  H", "  I"))
dev.off()

write_csv(fb_all_estimates, "data/estimates/fb_estimates_binomreg_rake_FINAL.csv")

###
