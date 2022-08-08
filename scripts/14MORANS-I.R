library(tidyverse)
library(ape)
library(lubridate)
library(ggpubr)

# load data
nbrs <- read_csv("data/county_neighbors_fips.txt", col_names = FALSE, col_types = "ii") %>% 
  rename(fips1 = X1, fips2 = X2) %>% 
  mutate(fips1 = ifelse(fips1 == 2270, 2158, ifelse(fips1 == 46113, 46102, fips1)),
         fips2 = ifelse(fips2 == 2270, 2158, ifelse(fips2 == 46113, 46102, fips2)))
# nbrs is bidirectional, but need to switch Kusivlak and Oglala counties to match fb data
# nbrs_unique <- nbrs %>% select(fips1) %>% distinct()
# do we need to exclude Alaska and Hawaii? nope already excluded in the txt file

fb_pp <- read_csv("data/estimates/fb_estimates_binomreg_FINAL.csv",
                  col_types = "iDdiiidddddddddddd") %>% 
  rename(p_est_fb = p_est) %>% 
  select(fips, month, mask_prop_most, mask_n_most, count, p_est_fb, .linpred, 
         population, density)
obs_prop_col <- fb_pp %>% select(fips, month, mask_prop_most)
# fb_pp_estimates has bias without raking
fb_pp_rake <- read_csv("data/estimates/fb_estimates_binomreg_rake_FINAL.csv", 
                       col_types = "iDiididddddddddddd") %>% 
  select(fips, month, successes, trials, mask_prop_most, p_est, population, density) %>% 
  rename(resamp_mask_prop_most = mask_prop_most) %>% 
  left_join(obs_prop_col, by = c("fips", "month"))
fb_pp_rake_debias <- read_csv("data/estimates/fb_estimates_binomreg_rake_debias_FINAL.csv", 
                              col_types = "iDiididdddddddddddddd") %>% 
  select(fips, month, successes, trials, mask_prop_most, biased_p, unbiased_p, 
         bias_correction, bias_ratio, population, density) %>% 
  rename(resamp_mask_prop_most = mask_prop_most) %>% 
  left_join(obs_prop_col, by = c("fips", "month"))
fb_comm_pp <- read_csv("data/estimates/fb_comm_estimates_binomreg_rake_FINAL.csv",
                       col_types = "iDiididdddddddd") 

onm_pp <- read_csv("data/estimates/onm_estimates_binomreg_FINAL.csv", 
                   col_types = "iDdiiidddddddddddd") %>% 
  filter(fips < 72000) %>%
  rename(onm_obs_prop = mask_grocery_very_prop) %>% 
  select(fips, month, onm_obs_prop, mask_grocery_very_suc, total_obs, onm_est, population, density)

### this is package method, no package method may give slightly different results
calculate_morans_i <- function(df, value_name){
  months <- unique(df$month)
  MIs <- data.frame(month = c(), mi = c(), sd = c())
  for(j in 1:length(months)){
    # rank needs to be made once fips not in the data have been removed and fips not in neighbor list have been removed
    fips_from_data <- df %>% filter(month == months[j]) %>% select(fips)
    fips_from_nbrs <- nbrs %>% select(fips1) %>% rename(fips = fips1) %>% distinct()
    fips_to_consider <- fips_from_data %>% inner_join(fips_from_nbrs, by = "fips")
    
    # make natural number ordering
    fips_list <- fips_to_consider[order(fips_to_consider$fips),]
    fips_list <- fips_list %>% mutate(fips_num = row_number()) %>% rename(fips1 = fips)
    
    # join with neighbors to make neighbors into natural number ordering
    nbr_df <- nbrs %>% inner_join(fips_list, by = "fips1") %>% rename(fips_num1 = fips_num)
    fips_list <- fips_list %>% rename(fips2 = fips1)
    nbr_df <- nbr_df %>% inner_join(fips_list, by = "fips2") %>% rename(fips_num2 = fips_num) %>% 
      select(fips_num1, fips_num2)
    
    # create matrix of weights with 1s and 0s
    num_fips <- fips_to_consider %>% nrow()
    nbr_mat <- as.matrix(nbr_df)
    mat <- matrix(0, num_fips, num_fips)
    mat[nbr_mat] <- 1 # uses fips_num as matrix row and column indices
    
    month_data <- df %>% filter(month == months[j]) %>% select(fips, eval(value_name))
    values <- fips_to_consider %>% inner_join(month_data, by = "fips") %>% 
      pull(eval(value_name))
    MI <- Moran.I(values, mat)
    MIs <- MIs %>% bind_rows(data.frame(month = months[j], mi = MI[1], sd = MI[3]))
  }
  return(MIs)
}

plot_morans_i <- function(df, lower, upper){
  out <- df %>% ggplot(aes(x = month, y = observed, group = 1)) +
    geom_line() +
    geom_point() + 
    labs(x = "Month", y = "Moran's I") +
    ylim(lower, upper) +
    theme_classic() +
    geom_ribbon(aes(ymax = observed + sd, ymin = observed - sd), alpha = 0.25)
  return(out)
}

# fb_pp_mi <- calculate_morans_i(fb_pp, "p_est_fb")
# p1 <- plot_morans_i(fb_pp_mi, 0.6, 0.8)
# 
# fb_pp_rake_mi <- calculate_morans_i(fb_pp_rake, "p_est")
# p2 <- plot_morans_i(fb_pp_rake_mi, 0.6, 0.8)

fb_pp_rake_debias_mi <- calculate_morans_i(fb_pp_rake_debias, "unbiased_p")
p3 <- plot_morans_i(fb_pp_rake_debias_mi, 0.6, 0.8)

pdf("figures/supplement/fb-morans-i.pdf", height = 4, width = 4)
ggarrange(p3, ncol = 1)
dev.off()

 # MI ON FB RESIDUALS ---------------------------------------------------------------------------
fb_pp_resid <- fb_pp %>% mutate(diff = mask_prop_most - p_est_fb)
fb_pp_rake_resid <- fb_pp_rake %>% mutate(diff = mask_prop_most - p_est)
fb_pp_rake_debias_resid <- fb_pp_rake_debias %>% mutate(diff = mask_prop_most - unbiased_p)

fb_pp_resid_mi <- calculate_morans_i(fb_pp_resid, "diff")
p4 <- plot_morans_i(fb_pp_resid_mi, 0.2, 0.6)

fb_pp_rake_resid_mi <- calculate_morans_i(fb_pp_rake_resid, "diff")
p5 <- plot_morans_i(fb_pp_rake_resid_mi, 0.2, 0.6)

fb_pp_rake_debias_resid_mi <- calculate_morans_i(fb_pp_rake_debias_resid, "diff")
p6 <- plot_morans_i(fb_pp_rake_debias_resid_mi, 0.2, 0.6)

pdf("figures/supplement/fb-morans-i-residuals.pdf", height = 4, width = 12)
ggarrange(p4, p5, p6, labels = "AUTO", ncol = 3)
dev.off()

# ONM --------------------------------------------------------------------------------
onm_pp_mi <- calculate_morans_i(onm_pp, "onm_est")
onm_plot <- plot_morans_i(onm_pp_mi, lower = 0.6, upper = 0.8)

pdf("figures/supplement/onm-morans-i.pdf", height = 4, width = 4)
onm_plot
dev.off()

onm_pp_resid <- onm_pp %>% mutate(diff = mask_grocery_very_suc - onm_est)
onm_pp_resid_mi <- calculate_morans_i(onm_pp_resid, "diff")
onm_resid_plot <- plot_morans_i(onm_pp_resid_mi, lower = 0.2, upper = 0.6)

pdf("figures/supplement/onm-morans-i-residuals.pdf", height = 4, width = 4)
onm_resid_plot
dev.off()

### no package method
calculate_morans_i_np <- function(data, var_name){
  months <- unique(data$month)
  morans_i <- rep(NA, length(months))
  for(j in 1:length(months)){
    df <- data %>% filter(month == months[j]) %>% select(fips, eval(var_name)) %>% rename(fips1 = fips)
    
    # merge neighbors and fb data to get list of neighboring masking values
    df_neigh_mask <- nbrs %>% inner_join(df, by = "fips1") %>% rename(p1 = eval(var_name))
    df <- df %>% rename(fips2 = fips1)
    df_neigh_mask <- df_neigh_mask %>% inner_join(df, by = "fips2") %>% rename(p2 = eval(var_name))
    
    # calculate spatial autocorrelation
    num_counties <- df_neigh_mask %>% select(fips1) %>% distinct() %>% nrow() # N
    mean_mask <- df %>% pull(eval(var_name)) %>% mean() # global mean
    total_weight <- df_neigh_mask %>% nrow() # W, number of county pairs
    numer <- df_neigh_mask %>% mutate(diff_1 = p1 - mean_mask,
                                      diff_2 = p2 - mean_mask,
                                      prod = diff_1 * diff_2) %>% # sum_i sum_j w_{ij} (x_i - x_bar) (x_j - x_bar)
      select(prod) %>% sum() # weights are binary, if neighbor or not, so we aren't calculating product when weight is 0, not neighbors
    df$orig_prop <- df %>% pull(eval(var_name))
    denom <- df %>% mutate(diff = orig_prop - mean_mask,
                           square = diff^2) %>% 
      select(square) %>% sum()
    morans_i[j] <- (num_counties/total_weight) * (numer/denom) # normalize
  }
  
  MI <- data.frame(month = months, mi = morans_i)
  return(MI)
}

plot_morans_i_np <- function(data, upper, lower){
  data %>% ggplot(aes(x = month, y = mi, group = 1)) + 
    geom_point() +
    geom_line() +
    ylim(upper, lower) +
    labs(x = "Month", y = "Moran's I")
}

fb_pp_mi_np <- calculate_morans_i_np(fb_pp, "p_est_fb")
p7 <- plot_morans_i_np(fb_pp_mi_np, 0.5, 0.7)

fb_pp_rake_mi_np <- calculate_morans_i_np(fb_pp_rake, "p_est")
p8 <- plot_morans_i_np(fb_pp_rake_mi_np, 0.5, 0.7)

fb_pp_rake_debias_mi_np <- calculate_morans_i_np(fb_pp_rake_debias, "unbiased_p")
p9 <- plot_morans_i_np(fb_pp_rake_debias_mi_np, 0.5, 0.7)

pdf("figures/supplement/fb-morans-i-no-package.pdf", height = 4, width = 12)
ggarrange(p7, p8, p9, ncol = 3)
dev.off()



