library(tidyverse)
library(readr)
library(ggplot2)
library(lubridate)

# from here: https://data.cdc.gov/Policy-Surveillance/U-S-State-and-Territorial-Public-Mask-Mandates-Fro/62d6-pm5i
county_data <- read_csv("data/U.S._State_and_Territorial_Public_Mask_Mandates_From_April_10__2020_through_August_15__2021_by_County_by_Day.csv") 
state_data <- read_csv("data/U.S._State_and_Territorial_Public_Mask_Mandates_From_April_8__2020_through_August_15__2021_by_State_by_Day.csv")
# looks like order_code = 2 denotes that mandate applies in specific counties, which are given by name, 1 is no mandate, 3 is state mandate

clean_data <- county_data %>% 
  mutate(Date = mdy(date),
         fips = as.integer(paste0(FIPS_State, ifelse(nchar(FIPS_County) == 3, FIPS_County, 
                                                     ifelse(nchar(FIPS_County) == 2, paste0("0", FIPS_County),
                                                            paste0("00", FIPS_County)))))) %>% 
  rename(public_mask = Face_Masks_Required_in_Public,
         source = Source_of_Action,
         state = State_Tribe_Territory) %>% 
  select(fips, state, Date, order_code, public_mask, source, URL, Citation) %>% 
  #filter(! is.na(public_mask)) %>% 
  filter(fips < 60000) %>% 
  mutate(is_mandate = ifelse(order_code == 1, 1, 0)) # 1 is yes mandate, 2 is no

clean_state_data <- state_data %>%
  mutate(Date = mdy(date)) %>%
  rename(public_mask = Face_Masks_Required_in_Public,
         source = Source_of_Action,
         state = State_Tribe_Territory) %>%
  select(state, Date, order_code, public_mask, Specific_counties, source, URL, Citation) %>%
  mutate(is_state_mandate = ifelse(order_code == 3, 1, 0), # 1 is no mandate, 2 is some counties, 3 is yes
         is_county_spec_mandate = ifelse(order_code == 2, 1, 0))

# 1 means yes and 2 means no masking in public with order_code
# with is_mandate variable 1 means mandate, 0 means no mandate

# states that seem to drop mandates in april/may
# AL, AR, CO, IN, KS, NH, UT, WI
# a little later but still may, LA, MN

# states that change, MS excluded all over the place
clean_data %>% filter(state %in% c("AL", "AR", "CO", "DC", "DE", "IA", "IL", "IN", "KS", "KY", "LA", "MA", "MD",
                                   "ME", "MI", "MN", "MT", "NC", "NH", "NJ", "NY", "OH", "OR", "PA", "TX", "UT",
                                   "VA", "VT", "WI", "WV", "WY")) %>% 
  ggplot(aes(x = Date, y = is_mandate, group = fips, col = state)) +
  geom_line() + 
  scale_x_date(date_labels = "%B", date_breaks = "1 month", limits = c(as.Date("2021-02-01"), as.Date("2021-07-01"))) +
  scale_y_continuous(limits=c(0, 1), breaks = c(0, 1))

pdf("figures/supplement/mask-mandates-hist.pdf", height = 4, width = 12)
clean_data %>% 
  filter(Date %within% interval(start = "2021-02-01", end = "2021-08-15")) %>% 
  group_by(fips) %>% 
  slice(which.max(is_mandate == 0)) %>%  # returns date of first TRUE, so when mandate is lifted
  ungroup() %>% filter(! is.na(public_mask)) %>% 
  filter(! state %in% c("AK", "AZ", "FL", "GA", "ID", "MO", "NE", "OK", "SC", "SD", "TN")) %>%  # filter out states that never change 
  ggplot() + 
  geom_histogram(aes(x = Date, fill = state), alpha = 0.8, col = "white") +
  geom_vline(xintercept = as.Date("2021-06-01"), lty = "dashed") +
  scale_x_date(date_labels = "%m-%Y", breaks = "1 month") +
  #ggtitle("AK, AZ, FL, GA, ID, MO, NE, OK, SC, SD, TN never had mask mandates and are excluded") +
  guides(fill=guide_legend(ncol=4)) +
  labs(y = "Number of counties lifting mask mandates", fill = "State")
dev.off()

# calculating proportion of mandates lifted by a certain time
tot_counties <- clean_data %>% select(fips) %>% distinct() %>% nrow()
tot_counties_w_mandate_ever <- clean_data %>% 
  filter(! state %in% c("AK", "AZ", "FL", "GA", "ID", "MO", "NE", "OK", "SC", "SD", "TN")) %>% select(fips) %>% 
  distinct() %>% nrow()
tot_counties_w_mandate_lifted_before_aug21 <- clean_data %>% 
  filter(! state %in% c("AK", "AZ", "FL", "GA", "ID", "MO", "NE", "OK", "SC", "SD", "TN")) %>%
  filter(! fips %in% c(28055, 28125)) %>%  # MS counties that dropped before 1/1/21
  filter(Date %within% interval(start = "2021-01-10", end = "2021-08-15")) %>% group_by(fips) %>% 
  slice(which.max(is_mandate == 0)) %>% ungroup() %>% filter(is_mandate == 0) %>% 
  select(fips) %>% distinct() %>% nrow()
num_lifted <- clean_data %>% 
  filter(Date %within% interval(start = "2021-01-10", end = "2021-08-15")) %>% 
  filter(! state %in% c("AK", "AZ", "FL", "GA", "ID", "MO", "NE", "OK", "SC", "SD", "TN")) %>% # I think need to filter out these states that never had mandates again
  filter(! fips %in% c(28055, 28125)) %>%  # MS counties that dropped before 1/1/21
  group_by(fips) %>% 
  slice(which.max(is_mandate == 0)) %>%  # returns date of first TRUE, so when mandate is lifted
  ungroup() %>% 
  filter(is_mandate == 0) %>% 
  filter(Date %within% interval(start = "2021-01-10", end = "2021-04-30")) %>% 
  mutate(lift_date = if_else(fips == 28021 | fips == 28143, as_date("2021-03-03"), Date)) %>%  
  # this is for fips 28021, 28143 that impose mandates on 1/15 and drop later, have to use if_else not ifelse
  # how many lifted before may 1
  select(fips) %>% distinct() %>% nrow()

frac_lifted_all <- num_lifted / tot_counties
frac_lifted_w_mandate_ever <- num_lifted / tot_counties_w_mandate_ever
frac_lifted_w_mandate_lifted_before_aug21 <- num_lifted / tot_counties_w_mandate_lifted_before_aug21
