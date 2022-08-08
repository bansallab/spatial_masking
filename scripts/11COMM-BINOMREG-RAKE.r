library(brms)
library(cmdstanr)
library(bayesplot)
library(insight)
library(tidybayes)
library(tidyverse)
#library(tidylog)
library(ggpubr)
options(brms.backend = "cmdstanr")

### read in data
# from here: https://covid19.census.gov/datasets/USCensus::average-household-size-and-population-density-county/explore?location=4.945434%2C0.315550%2C1.99&showTable=true
pop_density <- read_csv("data/Average_Household_Size_and_Population_Density_-_County.csv",
                        col_types = "ddiddccddddddiddccccddd") %>% 
  select(GEOID, B01001_001E, B01001_calc_PopDensity) %>% 
  rename(fips = GEOID, population = B01001_001E, density = B01001_calc_PopDensity)

resampled_df <- read_csv("data/fb_resampled_comm_data_from_raking_age_sex_educ_FINAL_corrected.csv", 
                         col_types = "iDiid")

# if excluding fips with pareto k's >= 0.7
# resampled_df <- resampled_df %>% filter(! fips %in% c(4019, 6037, 6071, 12071,
#                                                       12103, 36005, 40143, 41039,
#                                                       48113, 48201, 48439, 53033))
# added 48113, lost 45045 on 8/5/22

# remove this mutate and do earlier when rerun
to_model <- resampled_df %>% 
  left_join(pop_density, by = "fips") %>% 
  ungroup() %>% 
  filter(month != "2020-11-01") %>% 
  mutate(z_log10_density = (log10(density) - mean(log10(density)))/sd(log10(density)))


# FORMULA
m1_formula <- bf(successes | trials(trials) ~ 1 + z_log10_density,
                 family = binomial())

### prior
# set sample_prior = "only" so that data are not used
# think we can just test once because it doesn't use the data so should not vary across months
m1p <- brm(formula = m1_formula,
           data = to_model %>% filter(month == "2021-02-01"),
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
    mutate(comm_p_est = plogis(.linpred)) %>%
    mean_qi() %>%
    ungroup() %>%
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
                             "fit_comm_binomreg_rake"),
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
    fb_all_estimates <- make_output_dfs(model, (to_model %>% filter(month == months[i])))
  }else{
    fb_all_estimates <- fb_all_estimates %>% 
      bind_rows(make_output_dfs(model, (to_model %>% filter(month == months[i]))))
  }
}

pdf("figures/for-paper/supplement/posterior-predictive-checks-comm.pdf", height = 8, width = 12)
ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 3, 
          labels = c("  A", "  B", "  C", "  D", "  E", "  F"))
dev.off()

write_csv(fb_all_estimates, "data/estimates/fb_comm_estimates_binomreg_rake_FINAL.csv")

### PARETO Ks ### 

make_output_dfs_short <- function(model, fb_data){
  predictions <- as.data.frame(predict(model)) %>% cbind(fips = fb_data$fips)
  predictions_p <- data.frame(logit_p_est = colMeans(posterior_linpred(model))) %>%
    cbind(fips = fb_data$fips)
  out <- predictions %>% left_join(predictions_p, by = "fips") %>%
    left_join(fb_data, by = "fips") %>%
    mutate(p_est = plogis(logit_p_est))
  
  return(out)
}

months <- unique(to_model$month)
for(i in 1:length(months)){
  pks <- make_output_dfs_short(readRDS(paste0("data/fits/fb_", gsub("-", "_", substr(months[i], 3, 7)),
                                              "fit_comm_binomreg_rake.rds")),
                               to_model %>% filter(month == months[i])) %>% 
    mutate(pareto_ks = eval(parse(text = paste0("m", gsub("-", "_", substr(months[i], 3, 7)),
                                                "_pareto_ks"))))
  assign(paste0("est_", gsub("-", "_", substr(months[i], 3, 7))), pks)
  
}

all_ests <- est_20_12 %>% bind_rows(est_21_01, est_21_02, est_21_03, 
                                    est_21_04, est_21_05) %>% 
  group_by(month) %>% mutate(x = row_number()) %>% ungroup()

high_vals <- all_ests %>% filter(pareto_ks >= 0.7) %>% select(fips) %>% distinct() %>% pull(fips)


# old model testing
# library(viridis)
# all_ests %>% 
#   filter(month == "2021-02-01") %>% 
#   ggplot(aes(x = x, y = z_log10_density)) +
#   geom_point(aes(col = fips), alpha = 0.75) +
#   scale_color_viridis(direction = -1)
# all_ests %>% filter(pareto_ks > 0.01) %>% 
#   ggplot(aes(x = x, y = prop_mask_comm)) +
#   geom_point(aes(col = pareto_ks>0.5), alpha = 0.75) +
#   facet_wrap(~month, ncol = 1) +
#   scale_color_viridis(direction = -1, discrete = TRUE)
# all_ests %>% filter(pareto_ks > 0.01) %>% 
#   ggplot(aes(x = x, y = p_est)) +
#   geom_point(aes(col = pareto_ks>0.5), alpha = 0.75) +
#   facet_wrap(~month, ncol = 1) +
#   scale_color_viridis(direction = -1, discrete = TRUE)
# all_ests %>% filter(pareto_ks > 0.01) %>% 
#   ggplot(aes(x = x, y = prop_mask_comm - p_est)) +
#   geom_point(aes(col = pareto_ks>0.5), alpha = 0.75) +
#   facet_wrap(~month, ncol = 1) 
# all_ests %>% filter(pareto_ks > 0.01) %>% 
#   ggplot(aes(x = x, y = num_resp)) +
#   geom_point(aes(col =  pareto_ks>0.5), alpha = 0.75) +
#   facet_wrap(~month, ncol = 1)
# # HERE IS THE PROBLEM
# library(readxl)
# library(viridis)
# urb_rur_codes <- read_excel("data/NCHSURCodes2013.xlsx") %>% rename(fips = `FIPS code`) %>% 
#   mutate(`CBSA 2012 pop` = as.integer(`CBSA 2012 pop`),
#          `County 2012 pop` = as.integer(`County 2012 pop`)) # can ignore warnings
# all_ests <- all_ests %>% left_join(urb_rur_codes, by = "fips") %>% 
#   mutate(ur_code = as.factor(`2013 code`))
# pdf(file = "figures/for-paper/supplement/comm-masking-outlier-fips.pdf", height = 5, width = 8)
# all_ests %>% ggplot(aes(x = num_resp, y = (prop_mask_comm - p_est), fill = ur_code, col = pareto_ks > 0.5)) + 
#   #, alpha = pareto_ks > 0.5)) + 
#   geom_point(pch = 21) + 
#   labs(x = "Number of responses in county month", y = "Residual proportion", fill = "Urb/Rur Class", col = "Bad pareto k") +
#   scale_fill_viridis(discrete = TRUE) +
#   scale_color_manual(values = c(rgb(0, 0, 0, alpha=0), "red")) +
#   xlim(0, 20000) +
#   ylim(-1, 1)
# dev.off()
# 
# all_ests %>% filter(`State Abr.` %in% c("FL", "AZ", "WA", "NY", "OK", "TX", "OR", "NV", "SC")) %>% 
#   ggplot(aes(x = num_resp, y = population, fill = ur_code, col = pareto_ks > 0.5)) + 
#   geom_point(pch = 21) + 
#   labs(x = "Number of responses in county month", fill = "Urb/Rur Class", col = "Bad pareto k") +
#   scale_fill_viridis(discrete = TRUE) +
#   scale_color_manual(values = c(rgb(0, 0, 0, alpha=0), "red")) +
#   facet_wrap(~`State Abr.`)
# 
# # other diagnostics of model ------------------------------------------------------------
# # observed v predicted
# pdf(file = "figures/for-paper/supplement/diagnostic-comm.pdf", height = 6, width = 8)
# fb_all_estimates %>% ggplot(aes(x = prop_mask_comm, y = p_est)) +
#   geom_point(col = "dimgrey", alpha = 0.5) + 
#   geom_abline() + 
#   facet_wrap(~month) +
#   labs(x = "Observed", y = "Predicted")
# dev.off()
# 
# # observed v residual
# pdf(file = "figures/for-paper/supplement/diagnostic-comm-resid.pdf", height = 6, width = 8)
# fb_all_estimates %>% mutate(resid = prop_mask_comm - p_est) %>% 
#   ggplot(aes(x = prop_mask_comm, y = resid)) +
#   geom_point(col = "dimgrey", alpha = 0.5) + 
#   facet_wrap(~month) +
#   labs(x = "Observed", y = "Residual") +
#   ylim(-1, 1)
# dev.off()
# 
# 
# all_ests_test <- all_ests %>% mutate(foi = as.factor(ifelse(fips %in% c(4019, 6037, 6065, 6071, 6085, 12071, 12103, 
#                                                                         12031, 32031, 36005, 40143, 41039, 45045, 
#                                                                         48113, 48201, 48439, 53033), 1, 0)))
# 
# all_ests_test %>% ggplot(aes(x = prop_mask_comm, z_log10_density)) +
#   geom_point(aes(col = foi, alpha = foi)) +
#   facet_wrap(~month) +
#   scale_alpha_discrete(range = c(0.1, 0.9))
# 
# all_ests_test %>% ggplot(aes(x = num_resp, z_log10_density)) +
#   geom_point(aes(col = foi, alpha = foi)) +
#   facet_wrap(~month) +
#   scale_alpha_discrete(range = c(0.5, 0.9)) +
#   scale_x_log10()
# 
# all_ests_test %>% ggplot(aes(x = Estimate, z_log10_density)) +
#   geom_point(aes(col = foi, alpha = foi)) +
#   facet_wrap(~month) +
#   scale_alpha_discrete(range = c(0.5, 0.9)) +
#   scale_x_log10()
# 
