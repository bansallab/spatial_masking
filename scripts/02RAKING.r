### Performs raking on age, gender, and education for use in masking,
###   vaccination, and community masking questions
### Juliana Taube
### Last updated 6/27/22

### GOAL: come up with new county-month weights which we will use to resample the data to feed into the brms partial pooling model
### APPROACH: to come up with these weights, use raking

### start with data from underlying population on a few demographic characteristics: age, sex, education (we don't have race)
### compare these distributions to the distributions in the data (county-month aggregated) --> either use survey package or anesrake

# first get fb data on age, sex, education, and group to county-month in desired duration I think, 
#   then we calc these distributions for county month
library(tidyr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(readr)
library(vroom)

path <- dir('INSERT PATH TO FACEBOOK DATA (including race/ethnicity)')[1:9]
# 2020-09 through 2021-05

# only read in mask columns, date, sample weight, fips code, vaccination
cols <- cols_only(D1 = col_integer(), # gender
                  D2 = col_integer(), # age
                  D8 = col_integer(), # education
                  raceethnicity = col_character(), # race and ethnicity
                  V1 = col_integer(),
                  V2 = col_integer(),
                  C14 = col_integer(), # masking 5 days
                  C14a = col_integer(), # masking 7 days
                  C16 = col_double(), # how many ppl in public wore masks
                  H2 = col_double(), # replaces c16 at end of may
                  StartDatetime = col_datetime(),
                  weight = col_double(),
                  fips = col_integer())

# load FB data
data <- paste0('INSERT PATH TO FACEBOOK DATA (including race/ethnicity)', path) %>%
  map(vroom, # maps function to a vector, vroom is faster version of read_csv
      delim = ',',
      col_types = cols) %>% # map passes back a list of each of the separate paths
  bind_rows() %>% 
  mutate(response_id = row_number())


clean_data <- data %>% 
  mutate(week = ymd(floor_date(StartDatetime, unit = 'week')),
         month = ymd(floor_date(StartDatetime, unit = 'month')),
         self_response = ifelse(is.na(C14), C14a, C14),
         comm_response = case_when(H2 == 1 ~ 5, # H2 is opposite order of C16
                                   H2 == 2 ~ 4,
                                   H2 == 3 ~ 3,
                                   H2 == 4 ~ 2,
                                   H2 == 5 ~ 1,
                                   is.na(H2) ~ C16)) %>% 
  filter(!is.na(fips), !is.na(D1), !is.na(D2), !is.na(D8)) %>%
  filter(fips < 60000) %>% # not including outside 50 + DC
  filter(D1 <= 2) %>% # ACS can't compare non-binary, prefer not to answer, or other genders; assuming gender == sex
  mutate(education = ifelse(D8 >= 6, 6, D8)) %>% # put professional degree, masters, doctorate in one education level
  rename(sex = D1, age = D2, vax = V1, vax_doses = V2) %>% 
  select(-C14, -C14a, -C16, -H2, -StartDatetime, -D8) 

# check correlation of race and education variables
# since both are categorical variables, I think we need to use a chi-squared test
# code from: https://datascience.stackexchange.com/questions/893/how-to-get-correlation-between-two-categorical-variable-and-a-categorical-variab
my_tbl <- clean_data %>% filter(!is.na(raceethnicity)) %>% 
  group_by(education, raceethnicity) %>% 
  summarise(num_obs = n()) %>% 
  pivot_wider(names_from = raceethnicity, values_from = num_obs)
total_obs <- clean_data %>% filter(!is.na(raceethnicity)) %>% nrow()

chi2 <- chisq.test(my_tbl, simulate.p.value = T) #, correct=F)
library(rcompanion)
cramerV(as.matrix(my_tbl), ci = T) # vector memory exhausted

v <- sqrt((chi2$statistic[[1]]/total_obs)/min(nrow(my_tbl) - 1, ncol(my_tbl) - 1))


# don't aggregate because anesrake wants to give a weight to each response
fb_data_to_rake <- clean_data %>% mutate(age_cat = as.factor(case_when(age == 1 ~ "age18to24",
                                                                       age == 2 ~ "age25to34",
                                                                       age == 3 ~ "age35to44",
                                                                       age == 4 ~ "age45to54",
                                                                       age == 5 ~ "age55to64",
                                                                       age == 6 ~ "age65to74",
                                                                       age == 7 ~ "age75up")),
                                         sex_cat = as.factor(case_when(sex == 1 ~ "male",
                                                                       sex == 2 ~ "female")),
                                         edu_cat = as.factor(case_when(education == 1 ~ "ltHS",
                                                                       education == 2 ~ "HS",
                                                                       education == 3 ~ "ltColl",
                                                                       education == 4 ~ "2yrColl",
                                                                       education == 5 ~ "4yrColl",
                                                                       education == 6 ~ "postColl")),
                                         # raceeth_cat = as.factor(case_when(raceethnicity == "NonHispanicAmericanIndianAlaskaNative" ~ "NHAIANNHPI",
                                         #                                   raceethnicity == "NonHispanicNativeHawaiianPacificIslander" ~ "NHAIANNHPI",
                                         #                                   raceethnicity == "NonHispanicMultipleOther" ~ "NHMultOther",
                                         #                                   raceethnicity == "NonHispanicWhite" ~ "NHWhite",
                                         #                                   raceethnicity == "NonHispanicAsian" ~ "NHAsian",
                                         #                                   raceethnicity == "NonHispanicBlackAfricanAmerican" ~ "NHBlack",
                                         #                                   raceethnicity == "Hispanic" ~ "HAll")),
                                         caseid = row_number()) %>% # may need to further combine Multiple other with the Alaskan and Hawaiian groups
  rename(given_weight = weight) %>% 
  mutate(fips = ifelse(fips == 2270, 2158, ifelse(fips == 46113, 46102, fips)))

write_csv(fb_data_to_rake, "data/fb_data_to_rake_age_sex_educ.csv")

fb_data_to_rake <- read_csv("data/fb_data_to_rake_age_sex_educ.csv", col_types = "iicidiiiDDiiifffi")

# we are using acs data for this, plain census data that Shweta sent didn't have education
library(tidycensus)
# INSERT CENSUS API KEY

# acs example
v18 <- load_variables(2018, "acs5", cache = TRUE)
View(v18)

### EXTRACTING DESIRED ACS DATA ### 
# need age, race, and sex
# denominator to calculate prop of county residents in each category should be population over 18
# FB categories c("18-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75+")
men_by_age <- get_acs(geography = "county", 
                      variables = c(male_18to19 = "B01001_007", male_20 = "B01001_008", male_21 = "B01001_009", 
                                    male_22to24 = "B01001_010", male_25to29 = "B01001_011", male_30to34 = "B01001_012", 
                                    male_35to39 = "B01001_013", male_40to44 = "B01001_014", male_45to49 = "B01001_015",
                                    male_50to54 = "B01001_016", male_55to59 = "B01001_017", male_60and61 = "B01001_018",
                                    male_62to64 = "B01001_019", male_65and66 = "B01001_020", male_67to69 = "B01001_021", 
                                    male_70to74 = "B01001_022", male_75to79 = "B01001_023", male_80to84 = "B01001_024", 
                                    male_85plus = "B01001_025"),
                      output = "wide",
                      year = 2018)

men_by_age_agg <- men_by_age %>% 
  mutate(male_18to24 = male_18to19E + male_20E+ male_21E + male_22to24E,
         male_25to34 =  male_25to29E + male_30to34E,
         male_35to44 = male_35to39E + male_40to44E,
         male_45to54 = male_45to49E + male_50to54E, 
         male_55to64 = male_55to59E + male_60and61E + male_62to64E,
         male_65to74 = male_65and66E + male_67to69E+ male_70to74E,
         male_75up = male_75to79E + male_80to84E + male_85plusE,
         male_18up = male_18to19E + male_20E+ male_21E + male_22to24E + male_25to29E + male_30to34E + male_35to39E + male_40to44E +
           male_45to49E + male_50to54E + male_55to59E + male_60and61E + male_62to64E + male_65and66E + male_67to69E+ male_70to74E +
           male_75to79E + male_80to84E + male_85plusE,
         .keep = "unused") %>% 
  select(-ends_with("M")) #selecting out margin of error estimates for this df

#Women by age
women_by_age <- get_acs(geography = "county", 
                        variables = c(female_18to19 = "B01001_031", female_20 = "B01001_032", female_21 = "B01001_033", 
                                      female_22to24 = "B01001_034", female_25to29 = "B01001_035", female_30to34 = "B01001_036", 
                                      female_35to39 = "B01001_037", female_40to44 = "B01001_038", female_45to49 = "B01001_039", 
                                      female_50to54 = "B01001_040", female_55to59 = "B01001_041", female_60and61 = "B01001_042", 
                                      female_62to64 = "B01001_043", female_65and66 = "B01001_044", female_67to69 = "B01001_045", 
                                      female_70to74 = "B01001_046", female_75to79 = "B01001_047", female_80to84 = "B01001_048", 
                                      female_85plus = "B01001_049"),
                        output = "wide",
                        year = 2018)

women_by_age_agg <- women_by_age %>% 
  mutate(female_18to24 = female_18to19E + female_20E+ female_21E + female_22to24E,
         female_25to34 =  female_25to29E + female_30to34E,
         female_35to44 = female_35to39E + female_40to44E,
         female_45to54 = female_45to49E + female_50to54E, 
         female_55to64 = female_55to59E + female_60and61E + female_62to64E,
         female_65to74 = female_65and66E + female_67to69E+ female_70to74E,
         female_75up = female_75to79E + female_80to84E + female_85plusE,
         female_18up = female_18to19E + female_20E+ female_21E + female_22to24E + female_25to29E + female_30to34E + female_35to39E + female_40to44E +
           female_45to49E + female_50to54E + female_55to59E + female_60and61E + female_62to64E + female_65and66E + female_67to69E+ female_70to74E +
           female_75to79E + female_80to84E + female_85plusE,
         .keep = "unused") %>% 
  select(-ends_with("M")) #selecting out margin of error estimates for this df

#Male and female separated age distributions
age_by_sex <- men_by_age_agg %>%
  left_join(women_by_age_agg)

#Total in each age group, percentage of population in that age group
age_dist <- age_by_sex %>% 
  mutate(all_18to24 = female_18to24 + male_18to24,
         all_25to34 = female_25to34 + male_25to34,
         all_35to44 = female_35to44 + male_35to44,
         all_45to54 = female_45to54 + male_45to54, 
         all_55to64 = female_55to64 + male_55to64,
         all_65to74 = female_65to74 + male_65to74,
         all_75up = female_75up + male_75up,
         all_18up = female_18up + male_18up,
         prop_18to24 = all_18to24 / all_18up,
         prop_25to34 = all_25to34 / all_18up,
         prop_35to44 = all_35to44 / all_18up,
         prop_45to54 = all_45to54 / all_18up, 
         prop_55to64 = all_55to64 / all_18up,
         prop_65to74 = all_65to74 / all_18up,
         prop_75up = all_75up / all_18up,
         .keep = "unused")

sex_dist <- age_by_sex %>% 
  mutate(male = male_18to24 + male_25to34 + male_35to44 + male_45to54 + male_55to64 + male_65to74 + male_75up,
         female = female_18to24 + female_25to34 + female_35to44 + female_45to54 + female_55to64 + female_65to74 + female_75up,
         all_18up = female_18up + male_18up,
         prop_male = male/all_18up,
         prop_female = female/all_18up,
         .keep = "unused")


# education attainment for pop 25+ or sex by educational attainment for pop 25+ available but i think we need 18+ to match fb which has sex and age
education_by_age_by_sex <- get_acs(geography = "county",
                                   variables = c(m_18to24_lt_9th = "B15001_004", m_18to24_lt_dip = "B15001_005",
                                                 m_18to24_hs = "B15001_006", m_18to24_lt_coll = "B15001_007",
                                                 m_18to24_2yrcoll = "B15001_008", m_18to24_4yrcoll = "B15001_009",
                                                 m_18to24_gt_coll = "B15001_010", m_25to34_lt_9th = "B15001_012",
                                                 m_25to34_lt_dip = "B15001_013", m_25to34_hs = "B15001_014",
                                                 m_25to34_lt_coll = "B15001_015", m_25to34_2yrcoll = "B15001_016",
                                                 m_25to34_4yrcoll = "B15001_017", m_25to34_gt_coll = "B15001_018",
                                                 m_35to44_lt_9th = "B15001_020", m_35to44_lt_dip = "B15001_021",
                                                 m_35to44_hs = "B15001_022", m_35to44_lt_coll = "B15001_023",
                                                 m_35to44_2yrcoll = "B15001_024", m_35to44_4yrcoll = "B15001_025",
                                                 m_35to44_gt_coll = "B15001_026", m_45to64_lt_9th = "B15001_028",
                                                 m_45to64_lt_dip = "B15001_029", m_45to64_hs = "B15001_030",
                                                 m_45to64_lt_coll = "B15001_031", m_45to64_2yrcoll = "B15001_032",
                                                 m_45to64_4yrcoll = "B15001_033", m_45to64_gt_coll = "B15001_034",
                                                 m_65up_lt_9th = "B15001_036", m_65up_lt_dip = "B15001_037",
                                                 m_65up_hs = "B15001_038", m_65up_lt_coll = "B15001_039",
                                                 m_65up_2yrcoll = "B15001_040", m_65up_4yrcoll = "B15001_041",
                                                 m_65up_gt_coll = "B15001_042", f_18to24_lt_9th = "B15001_045",
                                                 f_18to24_lt_dip = "B15001_046", f_18to24_hs = "B15001_047",
                                                 f_18to24_lt_coll = "B15001_048", f_18to24_2yrcoll = "B15001_049",
                                                 f_18to24_4yrcoll = "B15001_050", f_18to24_gt_coll = "B15001_051",
                                                 f_25to34_lt_9th = "B15001_053", f_25to34_lt_dip = "B15001_054",
                                                 f_25to34_hs = "B15001_055", f_25to34_lt_coll = "B15001_056",
                                                 f_25to34_2yrcoll = "B15001_057", f_25to34_4yrcoll = "B15001_058",
                                                 f_25to34_gt_coll = "B15001_059", f_35to44_lt_9th = "B15001_061",
                                                 f_35to44_lt_dip = "B15001_062", f_35to44_hs = "B15001_063",
                                                 f_35to44_lt_coll = "B15001_064", f_35to44_2yrcoll = "B15001_065",
                                                 f_35to44_4yrcoll = "B15001_066", f_35to44_gt_coll = "B15001_067",
                                                 f_45to64_lt_9th = "B15001_069", f_45to64_lt_dip = "B15001_070",
                                                 f_45to64_hs = "B15001_071", f_45to64_lt_coll = "B15001_072",
                                                 f_45to64_2yrcoll = "B15001_073", f_45to64_4yrcoll = "B15001_074",
                                                 f_45to64_gt_coll = "B15001_075", f_65up_lt_9th = "B15001_077",
                                                 f_65up_lt_dip = "B15001_078", f_65up_hs = "B15001_079",
                                                 f_65up_lt_coll = "B15001_080", f_65up_2yrcoll = "B15001_081",
                                                 f_65up_4yrcoll = "B15001_082", f_65up_gt_coll = "B15001_083"),
                                   output = "wide",
                                   year = 2018)

education_agg <- education_by_age_by_sex %>% 
  mutate(edu_lt_hs = m_18to24_lt_9thE + m_18to24_lt_dipE + m_25to34_lt_9thE + m_25to34_lt_dipE + m_35to44_lt_9thE + m_35to44_lt_dipE +
           m_45to64_lt_9thE + m_45to64_lt_dipE + m_65up_lt_9thE + m_65up_lt_dipE + f_18to24_lt_9thE + f_18to24_lt_dipE +
           f_25to34_lt_9thE + f_25to34_lt_dipE + f_35to44_lt_9thE + f_35to44_lt_dipE + f_45to64_lt_9thE + f_45to64_lt_dipE + 
           f_65up_lt_9thE + f_65up_lt_dipE,
         edu_hs = m_18to24_hsE + m_25to34_hsE + m_35to44_hsE + m_45to64_hsE + m_65up_hsE +
           f_18to24_hsE + f_25to34_hsE + f_35to44_hsE + f_45to64_hsE + f_65up_hsE,
         edu_lt_coll = m_18to24_lt_collE + m_25to34_lt_collE + m_35to44_lt_collE + m_45to64_lt_collE + m_65up_lt_collE + 
           f_18to24_lt_collE + f_25to34_lt_collE + f_35to44_lt_collE + f_45to64_lt_collE + f_65up_lt_collE,
         edu_2yrcoll = m_18to24_2yrcollE + m_25to34_2yrcollE + m_35to44_2yrcollE + m_45to64_2yrcollE + m_65up_2yrcollE +
           f_18to24_2yrcollE + f_25to34_2yrcollE + f_35to44_2yrcollE + f_45to64_2yrcollE + f_65up_2yrcollE,
         edu_4yrcoll = m_18to24_4yrcollE + m_25to34_4yrcollE + m_35to44_4yrcollE + m_45to64_4yrcollE + m_65up_4yrcollE +
           f_18to24_4yrcollE + f_25to34_4yrcollE + f_35to44_4yrcollE + f_45to64_4yrcollE + f_65up_4yrcollE,
         edu_gtcoll = m_18to24_gt_collE + m_25to34_gt_collE + m_35to44_gt_collE + m_45to64_gt_collE + m_65up_gt_collE +
           f_18to24_gt_collE + f_25to34_gt_collE + f_35to44_gt_collE + f_45to64_gt_collE + f_65up_gt_collE,
         .keep = "unused") %>% 
  select(-ends_with('M'))


# now need to get 18up denominator into this dataset
denom <- sex_dist %>% select(GEOID, all_18up)
education_dist <- education_agg %>% left_join(denom, by = c("GEOID")) %>% 
  mutate(prop_lt_hs = edu_lt_hs/all_18up,
         prop_hs = edu_hs/all_18up,
         prop_lt_coll = edu_lt_coll/all_18up,
         prop_2yrcoll = edu_2yrcoll/all_18up,
         prop_4yrcoll = edu_4yrcoll/all_18up,
         prop_gt_coll = edu_gtcoll/all_18up)

all_acs <- age_dist %>% left_join(sex_dist) %>% left_join(education_dist) %>% 
  left_join(race_by_hisp_dist) %>% 
  mutate(fips = as.integer(GEOID)) %>% 
  select(fips, NAME, all_18up, contains("prop"))

write_csv(all_acs, "data/acs_target_data_for_fb_raking.csv")
all_acs <- read_csv("data/acs_target_data_for_fb_raking.csv", col_types = "icidddddddddddddddddddddd")

# NOW WE HAVE FACEBOOK RESPONSES AND ACS DATA WITH PROPORTIONS FOR GENDER/SEX, AGE, AND EDUCATION ATTAINMENT
# GET THINGS INTO RAKING FORMAT AND DO THE WEIGHTING!

# inputter is a list of target values, list contains vectors 
# what's tricky is that we have different targets for each county
# so do we need to apply anesrake to each county-month using a for loop? yes I think so

# for a given county month we need a list - this will be constant across months for target county but dataframe will change 
# now need data frame that matches targets

# need columns of this facebook data to be my desired attributes
# try collapsing age categories
fb_data_to_rake <- fb_data_to_rake %>% mutate(age_cat_col = as.factor(case_when(age_cat == "age18to24" ~ "age18to34",
                                                                                age_cat == "age25to34" ~ "age18to34", #"age25to34",
                                                                                age_cat == "age35to44" ~ "age35to54",
                                                                                age_cat == "age45to54" ~ "age35to54",
                                                                                age_cat == "age55to64" ~ "age55up",
                                                                                age_cat == "age65to74" ~ "age55up",
                                                                                age_cat == "age75up" ~ "age55up")),
                                              edu_cat_col = as.factor(case_when(edu_cat == "ltHS" ~ "HSOrLess",
                                                                                edu_cat == "HS" ~ "HSOrLess",
                                                                                edu_cat == "ltColl" ~ "SomeAllColl",
                                                                                edu_cat == "2yrColl" ~ "SomeAllColl",
                                                                                edu_cat == "4yrColl" ~ "SomeAllColl",
                                                                                edu_cat == "postColl" ~ "PostColl"))) #,

all_acs <- all_acs %>% mutate(prop_18to34 = prop_18to24 + prop_25to34,
                              prop_35to54 = prop_35to44 + prop_45to54,
                              prop_55up = prop_55to64 + prop_65to74 + prop_75up,
                              prop_HSOrLess = prop_lt_hs + prop_hs,
                              prop_SomeAllColl = prop_lt_coll + prop_2yrcoll + prop_4yrcoll,
                              prop_PostColl = prop_gt_coll, 
                              prop_allother = prop_notraking + prop_asian)

library(anesrake)
my_log <- file("data/raking_output.txt") # File name of output log
sink(my_log, append = TRUE, type = "output") # Writing console output to log file
sink(my_log, append = TRUE, type = "message")

all_acs_test <- all_acs
all_acs_test[all_acs_test == 0] <- 0.0000000001

counties <- fb_data_to_rake %>% pull(fips) %>% unique() %>% sort()
#counties <- counties[counties != 48301]
months <- fb_data_to_rake %>% pull(month) %>% unique()
all_weights <- data.frame(caseid = c(), weight = c())
non_convergent <- data.frame(month = c(), fips = c())
# don't rerun the above if doing months separately
for(i in 1:length(counties)){
  print(counties[i])
  # targets just need once per county
  acsoi <- all_acs_test %>% filter(fips == counties[i])
  # age
  agetarg <- c(acsoi$prop_18to34, acsoi$prop_35to54, acsoi$prop_55up)
  names(agetarg) <- c("age18to34", "age35to54", "age55up")
  # sex/gender
  sextarg <- c(acsoi$prop_male, acsoi$prop_female)
  names(sextarg) <- c("male", "female")
  # education
  edutarg <- c(acsoi$prop_HSOrLess, acsoi$prop_SomeAllColl, acsoi$prop_PostColl)
  names(edutarg) <- c("HSOrLess", "SomeAllColl", "PostColl")
  
  targets <- list(agetarg, sextarg, edutarg)
  names(targets) <- c("age_cat_col", "sex_cat", "edu_cat_col")
  for(j in 1:length(months)){
    print(months[j])
    if(months[j] == "2020-09-01" & counties[i] == 48311){next} # problems with this county
    fbdoi <- fb_data_to_rake %>% filter(month == months[j], fips == counties[i]) 
    if(nrow(fbdoi) >= 3){
      # do raking
      out <- anesrake(inputter = targets, dataframe = as.data.frame(fbdoi), caseid = fbdoi$caseid, 
                      type = "nolim", verbose = F, cap = 10, maxit = 200) 
      # nolim means include all variables so I don't get error that data are too close to target 
      outweights <- data.frame(caseid = out$caseid, weight = out$weightvec)
      all_weights <- all_weights %>% bind_rows(outweights)
      #print(summary(out))
      print(out$converge)
      if(!"Complete convergence was achieved" %in% out$converge){
        non_convergent <- non_convergent %>% bind_rows(data.frame(month = months[j], fips = counties[i]))
      }
      # use caseid to match up with original response
    }
    
  }
}
# can try upping the max it back higher than 100
write_csv(non_convergent, "data/fb_raking_nonconvergence_age_sex_educ_FINAL.csv")
non_convergent <- read_csv("data/fb_raking_nonconvergence_age_sex_educ_FINAL.csv") %>% 
  mutate(did_not_converge = 1)

closeAllConnections()
# are the counties with nas smaller sample size?
fb_weighted <- fb_data_to_rake %>% left_join(all_weights, by = "caseid")
write_csv(fb_weighted, "data/fb_raking_weights_age_sex_educ_FINAL.csv")

# correct nonconvergent weights to all be one
fb_corrected <- fb_weighted %>% left_join(non_convergent) %>% ungroup() %>% 
  rowwise() %>% 
  mutate(updated_weight = ifelse(is.na(did_not_converge), weight, 1))
write_csv(fb_corrected, "data/fb_raking_weights_age_sex_educ_FINAL_corrected.csv")

# checking if weights change over time
# pdf("figures/for-paper/supplement/weights-hist-by-month.pdf", height = 8, width = 12)
# fb_weighted %>% ggplot(aes(x = weight)) +
#   geom_histogram() +
#   facet_wrap(~month)
# dev.off()
# 
# pdf("figures/for-paper/supplement/weights-hist-by-month-converge.pdf", height = 8, width = 12)
# fb_weighted %>% left_join(non_convergent) %>% 
#   filter(is.na(did_not_converge)) %>% 
#   ggplot(aes(x = weight)) +
#   geom_histogram() +
#   facet_wrap(~month)
# dev.off()
# # 440 NAs are counties with <3 obs plus fips 48311 in Sept
# 
# summarized_data <- clean_data %>% 
#   mutate(fips = ifelse(fips == 2270, 2158, ifelse(fips == 46113, 46102, fips))) %>% 
#   group_by(month, fips) %>% 
#   summarise(num_obs = n())
# 
# non_convergent_samp_size <- non_convergent %>% left_join(summarized_data, by = c("month", "fips"))
# 
# # num county months included
# num_cnty_mos <- fb_data_to_rake %>% group_by(fips, month) %>% 
#   summarise(num_obs = n()) %>% 
#   mutate(is_incl = ifelse(num_obs >= 3, 1, 0),
#          is_incl = ifelse(month == "2020-09-01" & fips == 48311, 0, is_incl))
# sum(num_cnty_mos$is_incl)


