# Figures for paper
# Juliana Taube
# last run on 07/05/22
library(tidyverse)
library(choroplethr)
library(choroplethrMaps)
library(lubridate)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(readxl)
library(Hmisc)
library(vroom)
library(grid)
library(patchwork)

# read in data
fb_pp <- read_csv("data/estimates/fb_estimates_binomreg_FINAL.csv",
                         col_types = "iDdiiidddddddddddd") %>% 
  rename(p_est_fb = p_est) %>% 
  select(fips, month, mask_prop_most, mask_n_most, count, p_est_fb, .linpred, 
         population, density) %>% 
  filter(month %within% interval(start = ymd("2020-09-01"), end = ymd("2021-05-01")))
obs_prop_col <- fb_pp %>% select(fips, month, mask_prop_most) # so can get raw mask_prop_most into other dfs
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
onm_pp <- read_csv("data/estimates/onm_estimates_binomreg_FINAL.csv", col_types = "iDdiiidddddddddddd") %>% 
  filter(fips < 72000) %>%
  rename(onm_obs_prop = mask_grocery_very_prop) %>% 
  select(fips, month, onm_obs_prop, mask_grocery_very_suc, total_obs, onm_est, population, density) %>% 
  filter(month %within% interval(start = ymd("2020-09-01"), end = ymd("2021-05-01")))
urb_rur_codes <- read_excel("data/NCHSURCodes2013.xlsx") %>% rename(fips = `FIPS code`) %>% 
  mutate(`CBSA 2012 pop` = as.integer(`CBSA 2012 pop`),
         `County 2012 pop` = as.integer(`County 2012 pop`)) %>% 
  select(fips, `State Abr.`, `County name`, `CBSA title`, `2013 code`) %>% 
  rename(state = `State Abr.`, county = `County name`, area = `CBSA title`, ur_code = `2013 code`) %>%  # can ignore warnings
  mutate(ur_code = as.factor(ur_code))
  
### DATA PROCESSING AND BIAS MEASUREMENT -------------------------------------------------------------------
#(1) original data, (2) binomial regression, (3) binomial regression + raking, (4) binomial regression + raking + debias
map_cutoffs <- c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

# make original data choropleth
d1 <- fb_pp %>% mutate(value = cut(mask_prop_most, breaks = map_cutoffs, include.lowest = T),
                              region = ifelse(fips == 2158, 2270, ifelse(fips == 46102, 46113, fips)))
m1 <- CountyChoropleth$new(filter(d1, month == "2021-02-01"))
m1$ggplot_polygon <- geom_polygon(aes(fill = value), color = NA, size = 0.01)
m1$ggplot_scale <- scale_fill_brewer(name = "Raw mask prop.", palette = "YlGnBu", drop = FALSE)
p1 <- m1$render()
pdf("figures/supplement/observed-fb-masking-02-21.pdf", height = 4, width = 8)
p1
dev.off()

# getting all three datasets into one df
all_data_to_compare <- fb_pp %>% select(fips, month, mask_prop_most, p_est_fb) %>% 
  left_join(fb_pp_rake, by = c("fips", "month", "mask_prop_most")) %>% 
  rename(pp_prop = p_est_fb, rake_pp_prop = p_est) %>% 
  left_join(fb_pp_rake_debias, by = c("fips", "month", "mask_prop_most", "successes",
                                      "trials", "resamp_mask_prop_most", "population", 
                                      "density")) %>% 
  select(fips, month, trials, successes, mask_prop_most, resamp_mask_prop_most, 
         pp_prop, rake_pp_prop, biased_p, unbiased_p, bias_correction, 
         population, density) %>% 
  rename(rake_pp_uncorrected_prop = biased_p, rake_pp_debiased_prop = unbiased_p, 
         raw_prop = mask_prop_most)

# make new variables (differences and residuals) for plotting
all_data_to_compare2 <- all_data_to_compare %>% mutate(pp_diff = pp_prop - raw_prop,
                                                       pp_rake_diff = rake_pp_prop - raw_prop,
                                                       pp_rake_debias_diff = rake_pp_debiased_prop - raw_prop,
                                                       pp_resid = raw_prop - pp_prop,
                                                       pp_rake_resid = raw_prop - rake_pp_prop,
                                                       pp_rake_debias_resid = raw_prop - rake_pp_debiased_prop)

diff_cutoffs <- c(-1, -0.75, -0.5, -0.25, -0.05, 0.05, 0.25, 0.5, 0.75, 1)

# make binomial regression only choropleth
d2 <- all_data_to_compare2 %>% 
  mutate(value = cut(pp_resid, breaks = diff_cutoffs, include.lowest = T),
         region = as.integer(ifelse(fips == 2158, 2270, ifelse(fips == 46102, 46113, fips))))
m2 <- CountyChoropleth$new(filter(d2, month == "2021-02-01"))
m2$ggplot_polygon <- geom_polygon(aes(fill = value), color = NA, size = 0.01)
m2$ggplot_scale <- scale_fill_brewer(name = "Mask prop.\ndifference", 
                                     palette = "RdYlBu", drop = FALSE, 
                                     guide_legend("none"))
m2 <- m2$render()
p2 <- m2 + labs(caption = "Binomial regression model") + 
  theme(plot.caption = element_text(hjust = 0.5, face = "bold", size = 30, vjust = 0))

# make binomial regression + raking choropleth
d3 <- all_data_to_compare2 %>% 
  mutate(value = cut(pp_rake_resid, breaks = diff_cutoffs, include.lowest = T),
         region = as.integer(ifelse(fips == 2158, 2270, ifelse(fips == 46102, 46113, fips))))
m3 <- CountyChoropleth$new(filter(d3, month == "2021-02-01"))
m3$ggplot_polygon <- geom_polygon(aes(fill = value), color = NA, size = 0.01)
m3$ggplot_scale <- scale_fill_brewer(name = "Mask prop.\ndifference", 
                                     palette = "RdYlBu", drop = FALSE, guide_legend("none"))
m3 <- m3$render()
p3 <- m3 + labs(caption = "with raking") + 
  theme(plot.caption = element_text(hjust = 0.5, face = "bold", size = 30, vjust = 0))

# make binomial regression + raking + debiasing choropleth
d4 <- all_data_to_compare2 %>% 
  mutate(value = cut(pp_rake_debias_resid, breaks = diff_cutoffs, include.lowest = T),
         region = as.integer(ifelse(fips == 2158, 2270, ifelse(fips == 46102, 46113, fips))))
m4 <- CountyChoropleth$new(filter(d4, month == "2021-02-01"))
m4$ggplot_polygon <- geom_polygon(aes(fill = value), color = NA, size = 0.01)
m4$ggplot_scale <- scale_fill_brewer(name = "Residual\nmask prop.", 
                                     palette = "RdYlBu", drop = FALSE, guide_legend("none"))
m4 <- m4$render()
p4 <- m4 + labs(caption = "with raking and debiasing") + 
  theme(plot.caption = element_text(hjust = 0.5, face = "bold", size = 30, vjust = 0))

# need to get the legend of an object I'm not plotting
mL <- CountyChoropleth$new(filter(d4, month == "2021-02-01"))
mL$ggplot_polygon <- geom_polygon(aes(fill = value), color = NA, size = 0.01)
mL$ggplot_scale <- scale_fill_brewer(name = "Residual masking proportion", 
                                     palette = "RdYlBu", drop = FALSE)
mL <- mL$render()
pL <- mL + theme(legend.position = "bottom",
                 legend.text = element_text(size = 22),
                 legend.title = element_text(hjust = 0.5, size = 28),
                 legend.key = element_rect(color = "dimgrey")) +
  guides(fill = guide_legend(nrow = 1, title.position = "bottom"))
legend <- cowplot::get_legend(pL)

pdf("figures/fig1.pdf", height = 6, width = 20, onefile = F)
# setEPS()
# postscript("figures/fig1.eps", height = 8, width = 12)
ggarrange(ggarrange(p2, p3, p4, ncol = 3, widths = c(4, 4, 4), labels = "AUTO",
                    font.label = list(size = 30)),
          legend, nrow = 2, heights = c(4, 1))
          #common.legend = T, legend.grob = get_legend(pL), legend = "bottom",
          #font.label = list(size = 30))
dev.off()

### SPATIAL HETEROGENEITY IN FB DATA -------------------------------------------------------------------
unbiased_cutoffs <- c(0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
fb_est_map <- fb_pp_rake_debias %>% mutate(region = as.integer(ifelse(fips == 2158, 2270, ifelse(fips == 46102, 46113, fips))),
                                           value = cut(unbiased_p, breaks = unbiased_cutoffs, include.lowest = TRUE))
# make choropleth showing spatial heterogeneity for single month
fb_oct <- CountyChoropleth$new(filter(fb_est_map, month == "2020-10-01"))
fb_oct$ggplot_polygon <- geom_polygon(aes(fill = value), color = NA, size = 0.01)
fb_oct$ggplot_scale <- scale_fill_brewer(name = "Masking prop.",
                                         palette = "YlGnBu", drop = FALSE, guide_legend("none"))
p3a <- fb_oct$render()
# make fake plot so can extract legend
fake <- CountyChoropleth$new(filter(fb_est_map, month == "2020-10-01"))
fake$ggplot_polygon <- geom_polygon(aes(fill = value), color = NA, size = 0.01)
fake$ggplot_scale <- scale_fill_brewer(name = "Mask\nprop.",
                                         palette = "YlGnBu", drop = FALSE)
fake <- fake$render() + theme(legend.position = "right",
                            legend.text = element_text(size = 16),
                            legend.title = element_text(hjust = 0.5, vjust = 1, size = 20),
                            legend.key = element_rect(color = "dimgrey"),
                            #legend.key.size = unit(1, "lines"),
                            legend.spacing.y = unit(0.1, 'cm')) +
  guides(fill = guide_legend(ncol = 1, title.position = "top"))
#fake_legend <- cowplot::get_legend(fake)
p3amod <- p3a + inset_element(cowplot::get_legend(fake), left = 1.8, bottom = 0.3, right = 0, top = 0.3)

# urban/rural differences box/violin plot
fb_pp_rake_debias %>% left_join(urb_rur_codes, by = "fips") %>% 
  ggplot(aes(x = ur_code, y = unbiased_p, fill = ur_code)) +
  geom_violin(alpha = 0.5) + 
  geom_boxplot(alpha = 0.5, width = 0.35) +
  scale_fill_viridis(direction = 1, end = 0.95, discrete = TRUE) +
  annotate(geom = "text", x = 3.5, y = 0, 
           label = expression("Urban      "   %<->%   "     Rural"), size = 7) + # makes arrow between
  xlab("NCHS Class") +
  ylab("Masking proportion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(vjust = 1),
        plot.margin = margin(2,0.1,0.8,.1, "cm")) -> p3b

# urban/rural statistics - STAT TEST 1
# to use a one-way anova, first need to examine variances/sds in each group
ur_data <- fb_pp_rake_debias %>% left_join(urb_rur_codes, by = "fips") 
ur_data %>% 
  group_by(ur_code) %>% 
  summarise(std_dev = sd(unbiased_p))
# aov_test <- aov(unbiased_p ~ ur_code, data = ur_data)
# # given that significant differences, use multiple comparisons to tease out
# TukeyHSD(aov_test)
# library(multcomp) # an alternative way
# summary(glht(aov_test, linfct = mcp(ur_code = "Tukey")))
# pairwise.t.test(ur_data$unbiased_p, ur_data$ur_code, # additional alternative
#                 p.adjust.method = "BH")
# # check homogeneity of variance
# plot(aov_test, 1)
# library(car)
# leveneTest(unbiased_p ~ ur_code, data = ur_data) # significant so variance is different
# # anova with no assumption of equal variances
# oneway.test(unbiased_p ~ ur_code, data = ur_data)
# pairwise.t.test(ur_data$unbiased_p, ur_data$ur_code,
#                 p.adjust.method = "BH", pool.sd = FALSE)
# # check normality assumption
# plot(aov_test, 2)
# aov_residuals <- residuals(object = aov_test)
# shapiro.test(x = aov_residuals[0:5000]) # Shapiro-Wilk test, too many data points, but significant on 1st 5000

# use kruskal wallis since assumptions not met
kruskal.test(unbiased_p ~ ur_code, data = ur_data)
pairwise.wilcox.test(ur_data$unbiased_p, ur_data$ur_code, p.adjust.method = "BH")

# figure 2, I think can ignore warning
pdf(file = "figures/fig2.pdf", height = 6, width = 14) #height = 8, width = 10)
ggarrange(p3amod, p3b, ncol = 2, widths = c(10, 4),
          labels = c("A", "B"), font.label = list(size = 28))
# ggarrange(p3a, fake_legend, p3b, nrow = 1, ncol = 3, widths = c(10, 1, 4.3), 
#           labels = c("A", "", "B"), font.label = list(size = 28))
invisible(dev.off())

# additional supplemental plot showing urban-rural trends by month
fb_pp_rake_debias %>% left_join(urb_rur_codes, by = "fips") %>%
  ggplot(aes(x = ur_code, y = unbiased_p, fill = ur_code)) +
  geom_violin(alpha = 0.5) + 
  geom_boxplot(alpha = 0.5, width = 0.35) +
  facet_wrap(~month) +
  scale_fill_viridis(direction = 1, end = 0.95, discrete = TRUE) +
  labs(x = "NCHS Urban-Rural Class", y = "Debiased raked mean prop. mask\nmost of the time or more") +
  theme_classic() +
  theme(legend.position = "none") -> by_month
# pdf("figures/supplement/fb-ur-code-by-month.pdf", height = 8, width = 8)
# by_month
# dev.off()


### TEMPORAL HETEROGENEITY IN FB DATA -------------------------------------------------------------------
# generate time series data
time_series_mean <- fb_pp_rake_debias %>%
  group_by(fips) %>% mutate(avg = mean(unbiased_p), 
                            avg_raw = mean(mask_prop_most),
                            sd = sd(unbiased_p),
                            sd_raw = sd(mask_prop_most)) %>% 
  ungroup() %>% 
  mutate(avg_diff = unbiased_p - avg,
         avg_diff_raw = mask_prop_most - avg_raw, 
         zscore = (unbiased_p - avg)/sd,
         zscore_raw = (mask_prop_most - avg_raw)/sd_raw,
         centered_date = as.Date(ifelse(month == "2020-09-01", "2020-09-15", ifelse(month == "2020-10-01", "2020-10-15",
                                 ifelse(month == "2020-11-01", "2020-11-15", ifelse(month == "2020-12-01", "2020-12-15",
                                 ifelse(month == "2021-01-01", "2021-01-15", ifelse(month == "2021-02-01", "2021-02-15",
                                 ifelse(month == "2021-03-01", "2021-03-15", ifelse(month == "2021-04-01", "2021-04-15",
                                                                                    "2021-05-15"))))))))))
        # centers date so plotted mid-month not at beginning of month
p1 <- time_series_mean %>% ggplot(aes(x = centered_date, y = zscore, group = fips, col = avg)) + 
  geom_line(alpha = 0.2) +
  theme_classic() + 
  labs(x = "Month", y = "Masking proportion (z-score)", col = "Mask prop.\ncounty mean") +
  scale_x_date(breaks = seq(as.Date("2020-09-15"), as.Date("2021-05-15"), by = "1 month"), 
               #date_labels ="%B",
               labels = c("Sep 2020", "Oct", "Nov", "Dec", "Jan 2021", "Feb", "Mar", "Apr", "May"),
               expand = c(0.08, 1)) + 
  #scale_y_continuous(expand = c(0.5, 0.01)) +
  scale_y_continuous(limits = c(-4, 1.75), breaks = c(-2, -1, 0, 1)) + # can change this
  scale_color_viridis(option = "G", direction = -1) +
  theme(legend.position = "right",
        legend.text = element_text(size = 16),
        legend.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(vjust = 0),
        axis.title.x = element_text(vjust = 0)) +
  guides(col = guide_colourbar(barwidth = 2, barheight = 10, title.hjust = 0))
#scale_color_distiller(palette = "PuBuGn", direction = 1)

# nyt data to plot below: https://github.com/nytimes/covid-19-data/blob/master/rolling-averages/us.csv
nyt_case_data <- read_csv("data/nyt-us-rolling-avg.csv", col_types = "Dciddidd") %>% 
  filter(date %within% interval(start = "2020-09-01", end = "2021-05-31")) # we really want this thru 5/31 but then xaxes don't line up
# vax_path <- "archived_vaccination_data/data_master_county_20210823.csv"
# vax_data <- read_csv(vax_path) %>%
#   mutate(week = floor_date(DATE, unit = 'week')) %>%
#   filter(CASE_TYPE == "Partial Coverage" | CASE_TYPE == 'Partial Coverage') %>% # includes completely vaxxed
#   group_by(week, COUNTY, CASE_TYPE) %>%
#   summarise(p_vax = mean(CASES*.01)) %>% # cases is a percentage, we want proportion
#   mutate(fips = as.numeric(COUNTY)) %>%
#   ungroup() %>%
#   select(week, fips, p_vax) %>%
#   group_by(week) %>%
#   summarise(cum_mean_vax = mean(p_vax, na.rm = TRUE)) %>%
#   filter(week %within% interval(start = "2020-09-01", end = "2021-05-31"))
# saveRDS(vax_data, "data/vax_data_for_time_series.RDS")
vax_data <- readRDS("data/vax_data_for_time_series.RDS")
vax_data %>% ggplot(aes(x = week, y = cum_mean_vax)) +
  geom_line(size = 1) +
  theme_classic() +
  labs(x = "Date", y = "Cumulative mean national proportion vaccinated\n(1 or more doses)") -> vaxs

# backfills vax data so have line across whole plot
test <- vax_data %>% complete(week = seq.Date(as.Date("2020-09-01"), as.Date("2021-01-03"), by="week"))
test <- test %>% mutate(cum_mean_vax = ifelse(is.na(cum_mean_vax), 0, cum_mean_vax))
test %>% ggplot(aes(x = week, y = cum_mean_vax)) +
  geom_line(size = 1) +
  theme_classic() +
  labs(x = "Date", y = "Cum. mean national\nprop. vaccinated\n(1 or more doses)") -> vaxs

### worry data
# G1 and C9
# path <- dir('INSERT FACEBOOK DATA PATH')[6:14]
# # 2020-09 through 2021-05
# # only read in masking column, date, sample weight, fips code
# cols <- cols_only(C9 = col_integer(),
#                   G1 = col_integer(),
#                   StartDatetime = col_datetime(),
#                   weight = col_double(),
#                   fips = col_integer(), 
#                   wave = col_integer())
# load FB data
# worry_data <- paste0('INSERT FACEBOOK DATA PATH', path) %>%
#   map(vroom, 
#       delim = ',',
#       col_types = cols) %>% 
#   bind_rows() 
# worry_df <- worry_data %>% filter(!(is.na(C9) & is.na(G1))) %>% 
#   filter(!is.na(fips)) %>% 
#   mutate(date = ymd_hms(StartDatetime),
#          day = ymd(floor_date(date, unit = "day")),
#          week = ymd(floor_date(date, unit = "week")),
#          very_worried = ifelse(C9 == 1, 1, 0),
#          somewhat_worried = ifelse(C9 <= 2, 1, 0)) %>% 
#   group_by(week) %>% 
#   summarise(very_worried_avg = mean(very_worried),
#             somewhat_worried_avg = mean(somewhat_worried))
# 
# saveRDS(worry_df, "data/worry_df_for_time_series.RDS")
worry_df <- readRDS("data/worry_df_for_time_series.RDS")

cases <- nyt_case_data %>% ggplot(aes(x = date, y = cases_avg/1000)) + # cases_avg is rolling average, smooths drops in weekend/holidays
  geom_line(size = 1, col = "firebrick") +
  geom_line(data = test, aes(x = week, y = cum_mean_vax * 641), size = 1, col = "dodgerblue") + # 625 was from max of nyt/1000/max of vax
  geom_line(data = worry_df, aes(x = day, y = somewhat_worried_avg * 200), size = 1, col = "violet") +
  theme_classic() +
  scale_x_date(expand = c(0,0), breaks = seq(as.Date("2020-09-15"), as.Date("2021-05-15"), by = "1 month"), date_labels ="%B") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

# inset plot with cases, vax, and worry
nyt_case_data %>% ggplot(aes(x = date, y = (cases_avg-mean(cases_avg))/sd(cases_avg))) + # cases_avg is rolling average, smooths drops in weekend/holidays
  geom_line(size = 1, col = "#1B9E77") +
  geom_line(data = test, aes(x = week, y = (cum_mean_vax-mean(cum_mean_vax))/sd(cum_mean_vax)), size = 1, col = "#D95F02") +
  geom_line(data = worry_df, aes(x = week, # weekly is smoother but can do daily, too
                                 y = (somewhat_worried_avg-mean(somewhat_worried_avg, na.rm = TRUE))/sd(somewhat_worried_avg, na.rm = TRUE)), 
            size = 1, col = "#7570B3") +
  theme_classic() +
  scale_x_date(expand = c(0,0), breaks = seq(as.Date("2020-09-15"), as.Date("2021-05-15"), by = "1 month"), date_labels ="%B") +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) -> cases

vp <- viewport(width = 0.73, height = 0.3, x = 0.435, y = 0.27)

pdf("figures/fig3.pdf", height = 6, width = 10)
print(p1)
print(cases, vp = vp)
dev.off()

time_series_compare <- time_series_mean %>% 
  select(fips, month, unbiased_p, avg, avg_raw, avg_diff, avg_diff_raw, zscore, zscore_raw) %>% 
  pivot_longer(c(zscore, zscore_raw), names_to = "est_type", values_to = "value") %>% 
  mutate(est_type = as.factor(est_type))
levels(time_series_compare$est_type) <- c("Modeled", "Raw")
time_series_compare %>% 
  ggplot(aes(x = month, y = value, group = fips)) + 
  geom_line(alpha = 0.2, col = "darkslateblue") +
  facet_wrap(~est_type) +
  theme_classic() + 
  labs(x = "Month", y = "Masking prop. (z-score)") -> ps
# pdf("figures/supplement/fb-raw-vs-model-time-series.pdf", height = 4, width = 8)
# ps
# dev.off()

### COMMUNITY MASKING -------------------------------------------------------------------
comm_data <- fb_comm_pp %>%
  rename(successes_comm = successes, trials_comm = trials) %>% 
  select(fips, month, successes_comm, trials_comm, comm_mask_prop_most, comm_p_est) %>% 
  left_join(fb_pp_rake_debias, by = c("fips", "month")) %>% 
  mutate(diff = comm_p_est - unbiased_p) %>% 
  left_join(urb_rur_codes, by = "fips")

comm_data$monthf <- factor(comm_data$month, levels = c("2020-12-01", "2021-01-01",
                                                        "2021-02-01", "2021-03-01",
                                                        "2021-04-01", "2021-05-01"))
levels(comm_data$monthf) <- c("Dec 2020", "Jan 2021", "Feb 2021", "Mar 2021",
                              "Apr 2021", "May 2021")

library(viridis)
pdf(file = "figures/fig4.pdf", height = 8, width = 14)
comm_data %>% 
  ggplot(aes(x = unbiased_p, y = comm_p_est, col = ur_code)) +
  geom_point(alpha = 0.5, size = 4) +
  geom_abline() +
  facet_wrap(~monthf, nrow = 2, ncol = 3) +
  labs(x = "Bias-corrected masking proportion", y = "Community-reported masking proportion", 
       col = "NCHS Class") +
  scale_color_viridis(discrete = TRUE, labels = c("1 (Most urban)", "2", "3", 
                                                  "4", "5", "6 (Most rural)")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 20),
        legend.title = element_text(hjust = 0.5, size = 22),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(vjust = 0),
        strip.text.x = element_text(size = 20), # makes facet labels bigger
        axis.title.x = element_text(vjust = 0)) 
# warning is for county months without unbiased p??
dev.off()

### exclude influential fips codes
fb_comm_pp_excl <- read_csv("data/estimates/fb_comm_estimates_binomreg_rake_FINAL_excl.csv",
                       col_types = "iDiididdddddddd") 
comm_data_excl <- fb_comm_pp_excl %>%
  rename(successes_comm = successes, trials_comm = trials) %>% 
  select(fips, month, successes_comm, trials_comm, comm_mask_prop_most, comm_p_est) %>% 
  left_join(fb_pp_rake_debias, by = c("fips", "month")) %>% 
  mutate(diff = comm_p_est - unbiased_p) %>% 
  left_join(urb_rur_codes, by = "fips")

comm_data_excl$monthf <- factor(comm_data_excl$month, 
                                levels = c("2020-12-01", "2021-01-01", "2021-02-01",
                                           "2021-03-01", "2021-04-01", "2021-05-01"))
levels(comm_data_excl$monthf) <- c("Dec 2020", "Jan 2021", "Feb 2021", "Mar 2021",
                              "Apr 2021", "May 2021")

library(viridis)
pdf(file = "figures/supplement/fig4-excl.pdf", height = 8, width = 14)
comm_data_excl %>% 
  ggplot(aes(x = unbiased_p, y = comm_p_est, col = ur_code)) +
  geom_point(alpha = 0.5, size = 4) +
  geom_abline() +
  facet_wrap(~monthf, nrow = 2, ncol = 3) +
  labs(x = "Bias-corrected masking proportion", y = "Community-reported masking proportion", 
       col = "NCHS Class") +
  scale_color_viridis(discrete = TRUE, labels = c("1 (Most urban)", "2", "3", 
                                                  "4", "5", "6 (Most rural)")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 20),
        legend.title = element_text(hjust = 0.5, size = 22),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(vjust = 0),
        strip.text.x = element_text(size = 20), # makes facet labels bigger
        axis.title.x = element_text(vjust = 0)) 
# warning is for county months without unbiased p??
dev.off()


### OTHER SUPPLEMENTARY RESULTS ---------------------------------------------------------

### ONM vs FB comparison ----------------------------------------------------------------
# (A) map of onm:fb ratio
fips_w_both <- fb_pp %>% full_join(onm_pp, by = c("fips", "month")) %>% # inner
  mutate(has_both = (!is.na(p_est_fb) & !is.na(onm_est))) %>% 
  group_by(fips) %>% 
  summarise(both_data = sum(has_both),
            total_data = n()) %>% ungroup() %>% 
  filter(both_data >= 5)  # county must have 5/9 observations in each dataset

ratio <- fb_pp %>% full_join(onm_pp, by = c("fips", "month")) %>% 
  inner_join(fips_w_both, by = "fips") %>% 
  mutate(ratio_diff = onm_est/p_est_fb) %>% 
  group_by(fips) %>% 
  summarise(avg_ratio = mean(ratio_diff, na.rm = TRUE))

cutoffs <- c(0.67, 0.71, 0.77, 0.83, 0.91, 1.1, 1.2, 1.3, 1.4, 1.5)
ratio <- ratio %>% mutate(region = ifelse(fips == 2158, 2270, ifelse(fips == 46102, 46113, fips)),
                          value = cut(avg_ratio, breaks = cutoffs, include.lowest = TRUE))
comp_map <- CountyChoropleth$new(ratio)
comp_map$ggplot_polygon <- geom_polygon(aes(fill = value), color = "gray", size = 0.01)  #color = NA, size = 1
comp_map$ggplot_scale <- scale_fill_brewer(palette = "RdBu", name = "ONM:FB", na.value = "dimgrey", drop = FALSE)
p1a <- comp_map$render()

# (B) time series of ONM and FB for a couple counties
ts <- fb_pp %>% inner_join(onm_pp, by = c("fips", "month")) %>% left_join(urb_rur_codes, by = "fips") %>% 
  select(fips, month, p_est_fb, onm_est, state, county, ur_code) %>% 
  filter(fips %in% c(48113, 23003, 53073)) %>%
  pivot_longer(c(onm_est, p_est_fb), values_to = "prop", names_to = "dataset") %>% 
  mutate(county_factor = ifelse(county == "Aroostook County", "Aroostook County, ME\nnon-core",
                                ifelse(county == "Dallas County", "Dallas County, TX\nlarge central metro",
                                       "Whatcom County, WA\nsmall metro")))
ts$dataset <- factor(ts$dataset, levels = c("p_est_fb", "onm_est")) # change order
levels(ts$dataset) <- c("FB", "ONM") # then change name
ts$county_factor <- factor(ts$county_factor, levels = c("Dallas County, TX\nlarge central metro",
                                                        "Whatcom County, WA\nsmall metro", "Aroostook County, ME\nnon-core")) # change order
levels(ts$county_factor) <- c("Dallas, TX\nlarge central metro", "Whatcom, WA\nsmall metro", "Aroostook, ME\nnon-core") # change name

ts %>% ggplot(aes(x = month, y = prop, col = county_factor, group = interaction(fips, dataset))) +
  geom_line(aes(lty = dataset), size = 1) +
  labs(x = "Month", y = "Masking proportion", lty = "Dataset", col = "County") +
  scale_color_viridis(direction = 1, discrete = TRUE, end = 0.9) +
  theme_classic() +
  theme(legend.position = c(0.52, 0.18),
        legend.box = "horizontal",
        legend.background = element_rect(fill="transparent")) -> p1b

# figure 1
pdf("figures/supplement/onm-vs-fb.pdf", height = 4, width = 12)
ggarrange(p1a, p1b, ncol = 2, nrow = 1, widths = c(6, 3.5), labels = "AUTO")
dev.off()

# accompanying supplemental figures
# time series different by urban/rural classification
ts_diff <- fb_pp %>% inner_join(onm_pp, by = c("fips", "month")) %>% left_join(urb_rur_codes, by = "fips") %>%
  select(fips, month, p_est_fb, onm_est, state, county, ur_code) %>%
  mutate(diff = onm_est - p_est_fb)
pdf("figures/supplement/onm-fb-diff-over-time-by-ur-code-nobias.pdf", height = 5, width = 5)
ts_diff %>% ggplot(aes(x = month, y = diff, group = fips)) +
  geom_line(col = "dimgrey", alpha = 0.5) +
  facet_wrap(~ur_code) +
  labs(x = "Month", y = "ONM - FB masking proportion") +
  theme_classic()
invisible(dev.off())

# time series correlations
corr_df <- fb_pp %>% select(fips, month, p_est_fb) %>% 
  full_join(onm_pp, by = c("fips", "month")) %>% 
  inner_join(fips_w_both, by = "fips")
# to calculate spearman CI from https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation
spearman_CI <- function(x, y, alpha = 0.05){
  rs <- cor(x, y, method = "spearman", use = "complete.obs")
  n <- sum(complete.cases(x, y))
  sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
}

# correlation across space
corrT <- corr_df %>% ungroup() %>% group_by(fips) %>% 
  summarise(corS = cor(onm_est, p_est_fb, method = "spearman", use = "complete.obs"),
            upper = spearman_CI(onm_est, p_est_fb)[2],
            lower = spearman_CI(onm_est, p_est_fb)[1],
            diff = upper - lower,
            corP = cor.test(onm_est, p_est_fb, method = "pearson", use = "complete.obs")$estimate,
            lowerP = cor.test(onm_est, p_est_fb, method = "pearson", use = "complete.obs")$conf.int[1],
            upperP = cor.test(onm_est, p_est_fb, method = "pearson", use = "complete.obs")$conf.int[2],
            diffP = upperP - lowerP) %>% 
  ungroup() %>% 
  rename(region = fips) %>% mutate(region = ifelse(region == 2158, 2270, ifelse(region == 46102, 46113, region)))

cutoffs <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
corrT_cut <- corrT %>% mutate(corS_cut = cut(corS, breaks = cutoffs, include.lowest = TRUE))
spearman <- CountyChoropleth$new(rename(.data = corrT_cut, value = corS_cut))
spearman$ggplot_polygon <- geom_polygon(aes(fill = value), color = "gray", size = 0.01)  #color = NA, size = 1
spearman$ggplot_scale <- scale_fill_brewer(name = "Time series\ncorrelation", palette = "RdYlBu", drop = FALSE)
#scale_fill_brewer(name = "Time series\ncorrelation", palette = "Greens", drop = FALSE)
spearman <- spearman$render()
p1c <- spearman + theme(plot.title = element_text(face = "bold", hjust = 0.5))

# correlation over time
# alternative option is package DescTools and function SpearmanRho with conf.level = 0.95
corr <- corr_df %>% group_by(month) %>% summarise(corS = cor(onm_est, p_est_fb, method = "spearman", use = "complete.obs"),
                                                  upper = spearman_CI(onm_est, p_est_fb)[2],
                                                  lower = spearman_CI(onm_est, p_est_fb)[1],
                                                  num_counties = n())
# complete.obs uses only complete observations, throws out NAs which is fine bc we should have limited how many

corr %>% ggplot(aes(x = month, y = corS, group = 1)) + 
  geom_line() +
  geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.3) + # decide if error_bar is better
  geom_point() +
  labs(y = "Spearman correlation", x = "Month") +
  theme(legend.position = "bottom") +
  theme_classic() +
  ylim(0.75, 1) -> p1d


pdf("figures/supplement/time-series-correlations.pdf", height = 4, width = 12)
ggarrange(p1c, p1d, ncol = 2, nrow = 1, widths = c(6, 3.5), labels = "AUTO")
dev.off()

# MODEL DIAGNOSTICS ------------------------------------------------------------------
# observed v predicted - > this is a residual so no need to make maps since they already make this comparison
# FB biased & unbiased
pdf(file = "figures/residuals/diagnostic-fb-nobias.pdf", height = 6, width = 8)
fb_pp %>% ggplot(aes(x = p_est_fb, y = mask_prop_most)) +
  geom_point(col = "dimgrey", alpha = 0.4) + 
  geom_abline() + 
  facet_wrap(~month) +
  labs(x = "Predicted (binomial model)", y = "Observed") +
  geom_smooth(method = "lm") #+
  # scale_x_sqrt() +
  # scale_y_sqrt()
dev.off()

# pdf(file = "figures/residuals/raking-residual.pdf", height = 6, width = 8)
# fb_pp_rake %>% ggplot(aes(x = mask_prop_most, y = resamp_mask_prop_most)) +
#   geom_point(col = "dimgrey", alpha = 0.5) + 
#   geom_abline() + 
#   facet_wrap(~month) +
#   labs(x = "Observed", y = "Resampled") +
#   geom_smooth(method = "lm") 
# dev.off()

pdf(file = "figures/residuals/diagnostic-fb-rake.pdf", height = 6, width = 8)
fb_pp_rake %>% ggplot(aes(x = p_est, y = mask_prop_most)) +
  geom_point(col = "dimgrey", alpha = 0.5) + 
  geom_abline() + 
  facet_wrap(~month) +
  labs(x = "Predicted (binomial model + raking)", y = "Observed") +
  geom_smooth(method = "lm") #+
  # scale_x_sqrt() +
  # scale_y_sqrt()
dev.off()
# FB de bias
pdf(file = "figures/residuals/diagnostic-fb-rake-debias.pdf", height = 6, width = 8)
fb_pp_rake_debias %>% ggplot(aes(x = unbiased_p, y = mask_prop_most)) +
  geom_point(col = "dimgrey", alpha = 0.5) + 
  geom_abline() + 
  facet_wrap(~month) +
  labs(x = "Predicted (binomial model + raking + debiasing)", y = "Observed") +
  geom_smooth(method = "lm") #+
  # scale_x_sqrt() +
  # scale_y_sqrt()
dev.off()

# ONM no bias
pdf(file = "figures/residuals/diagnostic-onm.pdf", height = 6, width = 8)
onm_pp %>% ggplot(aes(x = qlogis(onm_est), y = qlogis(onm_obs_prop))) +
  geom_point(col = "dimgrey", alpha = 0.5) + 
  geom_abline() + 
  facet_wrap(~month) +
  labs(x = "Predicted (binomial model)", y = "Observed") +
  geom_smooth(method = "lm") #+
  #ylim(0, 1) + xlim(0, 1)
dev.off()

# observed v residual
pdf(file = "figures/residuals/diagnostic-fb-nobias-resid.pdf", height = 6, width = 8)
fb_pp %>% mutate(resid = mask_prop_most - p_est_fb,
                 logit_resid = qlogis(mask_prop_most - p_est_fb),
                 log_resid = qlogis(mask_prop_most) - qlogis(p_est_fb)) %>% 
  ggplot(aes(x = p_est_fb, y = log_resid)) +
  geom_point(col = "dimgrey", alpha = 0.5) + 
  facet_wrap(~month) +
  labs(x = "Predicted", y = "Residual (difference of logits)") +
  geom_smooth(method = "lm") #+
dev.off()

pdf(file = "figures/residuals/diagnostic-fb-rake-resid.pdf", height = 6, width = 8)
fb_pp_rake %>% mutate(resid = mask_prop_most - p_est,
                      logit_resid = qlogis(mask_prop_most - p_est),
                      log_resid = qlogis(mask_prop_most) - qlogis(p_est)) %>% 
  ggplot(aes(x = p_est, y = log_resid)) +
  geom_point(col = "dimgrey", alpha = 0.5) + 
  facet_wrap(~month) +
  labs(x = "Predicted", y = "Residual (difference of logits)") +
  geom_smooth(method = "lm") #+

dev.off()
# FB no bias
pdf(file = "figures/residuals/diagnostic-fb-rake-debias-resid.pdf", height = 6, width = 8)
fb_pp_rake_debias %>% mutate(resid = mask_prop_most - unbiased_p,
                             logit_resid = qlogis(mask_prop_most - unbiased_p),
                             log_resid = qlogis(mask_prop_most) - qlogis(unbiased_p)) %>% 
  ggplot(aes(x = unbiased_p, y = log_resid)) +
  geom_point(col = "dimgrey", alpha = 0.5) + 
  facet_wrap(~month) +
  labs(x = "Predicted", y = "Residual (difference of logits)") +
  geom_smooth(method = "lm") #+

dev.off()

# ONM no bias
pdf(file = "figures/residuals/diagnostic-onm-resid.pdf", height = 6, width = 8)
onm_pp %>% mutate(resid = onm_obs_prop - onm_est,
                  logit_resid = qlogis(onm_obs_prop - onm_est),
                  log_resid = qlogis(onm_obs_prop) - qlogis(onm_est)) %>% 
  ggplot(aes(x = onm_est, y = log_resid)) +
  geom_point(col = "dimgrey", alpha = 0.5) + 
  facet_wrap(~month) +
  labs(x = "Predicted", y = "Residual (difference of logits)") +
  geom_smooth(method = "lm") 
dev.off()

# residuals versus sample size to check shrinkage
fb_pp %>% ggplot(aes(x = count, y = qlogis(mask_prop_most) - qlogis(p_est_fb))) +
  geom_point(col = "dimgrey", alpha = 0.5)
fb_pp_rake %>% ggplot(aes(x = trials, y = qlogis(mask_prop_most) - qlogis(p_est))) +
  geom_point(col = "dimgrey", alpha = 0.5)
fb_pp_rake_debias %>% ggplot(aes(x = trials, y = qlogis(mask_prop_most) - qlogis(unbiased_p))) +
  geom_point(col = "dimgrey", alpha = 0.5)

# figure out what's going on with onm
onm_pp %>% pivot_longer(cols = c(onm_est, onm_obs_prop), values_to = "prop", names_to = "type") %>% 
  ggplot(aes(x = month, y = prop, col = type)) +
  geom_jitter(alpha = 0.5)


# PLOTTING SAMPLE SIZE ---------------------------------------------------------------
#samp_cutoffs <- c(1, 10, 25, 50, 100, 500, 25000)
samp_cutoffs <- c(1, 10, 100, 1000, 10000, 100000)
samp_size <- fb_pp %>% mutate(region = ifelse(fips == 2158, 2270, ifelse(fips == 46102, 46113, fips))) %>% 
  mutate(value = cut(count, breaks = samp_cutoffs, include.lowest = TRUE))

months <- unique(fb_pp$month)
months_str <- c("September 2020", "October 2020", "November 2020", "December 2020",
                "January 2021", "February 2021", "March 2021", "April 2021", "May 2021")
for (i in seq(1:length(months))){
  pdf(file = paste0("figures/samp-size/fb-", months[i], "-sampsize.pdf"), height = 4, width = 7)
  
  sampsize <- CountyChoropleth$new(filter(samp_size, month == months[i]))
  sampsize$title <- months_str[i]
  sampsize$ggplot_polygon <- geom_polygon(aes(fill = value), color = NA, size = 0.01)  #color = NA, size = 1
  sampsize$ggplot_scale <- scale_fill_brewer(name = "Num. response", #direction =-1,
                                               palette = "Reds", drop = FALSE, na.value = "dimgrey")
  sampsize <- sampsize$render()
  plot <- sampsize + theme(plot.title = element_text(hjust=0.5, face = "bold"))
  print(plot)
  dev.off()
}

# RAW PROPORTIONS IN EACH CATEGORY OVER TIME --------------------------------------------
fb_raw <- read_csv("data/fb_processed.csv", col_types = "iDdidididdddi") %>% left_join(urb_rur_codes, by = "fips") %>% 
  select(fips, month, count, mask_prop_all, mask_prop_most_only, mask_prop_some_only, mask_prop_little_only, mask_prop_none_only,
         state, county, ur_code) %>% 
  mutate(name = paste0(county, ", ", state))
onm_raw <- read_csv("data/onm_processed.csv", col_types = "iDddiiiddd") %>% left_join(urb_rur_codes, by = "fips") %>% 
  select(fips, month, total_obs, mask_grocery_very_prop, mask_grocery_some_only_prop, 
         mask_grocery_notso_only_prop, mask_grocery_notatall_only_prop,
         state, county, ur_code) %>% 
  mutate(name = paste0(county, ", ", state)) 

to_plot <- fb_raw %>% 
  filter(fips %in% c(8031, 24031, 36067, 1003, 41019, 30031, 48499, 21097, 23021)) %>%
  pivot_longer(c(mask_prop_all, mask_prop_most_only, mask_prop_some_only, mask_prop_little_only, mask_prop_none_only),
               values_to = "mask_prop", names_to = "mask_level")
to_plot$mask_level <- factor(to_plot$mask_level, levels = c("mask_prop_all", "mask_prop_most_only", "mask_prop_some_only", 
                                                            "mask_prop_little_only", "mask_prop_none_only")) # change order
levels(to_plot$mask_level) <- c("All", "Most", "Some", "A little", "None") # then change name
pdf("figures/supplement/fb-raw-mask-props-over-time.pdf", height = 9, width = 18)
to_plot %>% ggplot(aes(x = month)) +
  geom_col(aes(y = mask_prop, fill = mask_level), position = position_fill(reverse = TRUE), width = 35) + 
  labs(x = "Month",  y = "Raw proportion masking",  fill = "Amount of the\ntime masking") +
  scale_fill_viridis(discrete = TRUE) +
  geom_text(aes(y = 1.01, label=count), vjust=0, size = 3, col = "dimgrey") +
  facet_wrap(~as.factor(name))
invisible(dev.off())

to_plot2 <- onm_raw %>% filter(fips %in% c(8031, 24031, 36067, 1003, 41019, 30031, 48499, 21097, 23021)) %>%
  pivot_longer(c(mask_grocery_very_prop, mask_grocery_some_only_prop, mask_grocery_notso_only_prop, mask_grocery_notatall_only_prop),
               values_to = "mask_prop", names_to = "mask_level")
to_plot2$mask_level <- factor(to_plot2$mask_level, levels = c("mask_grocery_very_prop", "mask_grocery_some_only_prop", 
                                                              "mask_grocery_notso_only_prop", "mask_grocery_notatall_only_prop")) # change order
levels(to_plot2$mask_level) <- c("Very", "Somewhat", "Not so", "Not at all") # then change name
pdf("figures/supplement/onm-raw-mask-props-over-time.pdf", height = 9, width = 18)
to_plot2 %>% ggplot(aes(x = month, y = mask_prop, fill = mask_level)) +
  geom_col(position = position_fill(reverse = TRUE), width = 35) + 
  labs(x = "Month",  y = "Raw proportion masking",  fill = "How likely\nto mask\nwhile grocery\nshopping") +
  scale_fill_viridis(discrete = TRUE) +
  geom_text(aes(y = 1.01, label=total_obs), vjust=0, size = 3, col = "dimgrey") +
  facet_wrap(~name)
dev.off()

# MOTIVATE BIAS CORRECTION - MAP OF TRUE VERSUS FB VAX DATA -------------------------------------
clean_data_full <- readRDS(file = "data/fb_bias_clean_data_full_rake.RDS") %>% 
  mutate(fips = as.numeric(levels(fips))[as.integer(fips)]) %>% 
  left_join(urb_rur_codes, by = "fips")
# do we want to use the modeled FB estimates or the raw for this comparison? modeled - so then presumably raked is ok and now fixed vax denominator
clean_data_full %>% mutate(ratio = p_fb/p_vax_18up) %>% 
  ggplot(aes(x = week, y = ratio, group = fips, col = state)) +
  geom_line() +
  facet_wrap(~state, scales = "free_y") +
  geom_hline(aes(yintercept = 0), lty = "dashed") +
  theme(legend.position = "none") +
  labs(x = "Week", y = "FB raw vax prop:true vax prop")
  

# run the following instead of rerunning model
fit.prob <- readRDS(file = "data/fb_bias_vax_ests_rake.RDS")
clean.fitted <- clean_data_full %>% 
  mutate(fitted_p_fb = fitted(fit.prob), # extracts fitted values from model
         logit_vax = qlogis(p_vax_18up),
         logit_fitted_fb = qlogis(fitted_p_fb),
         #logit_diff = logit_vax - logit_fitted)
         logit_diff = logit_fitted_fb - logit_vax,
         p_ratio = fitted_p_fb/p_vax_18up)

pdf(file = "figures/supplement/bias-correction-vax-motivation-rake-denom-fixed.pdf", height = 10, width = 12)
clean.fitted %>% # fitted_p_fb - p_vax
  ggplot(aes(x = week, y = p_ratio, group = fips, col = state)) +
  geom_line() +
  facet_wrap(~state, scales = "free_y") +
  geom_hline(aes(yintercept = 1), lty = "dashed") +
  theme(legend.position = "none") +
  labs(x = "Week", y = "FB modeled vax prop:true vax prop")
dev.off()



