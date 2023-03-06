# Data and code for "Spatiotemporal trends in self-reported mask-wearing behavior in the United States: Analysis of a large cross-sectional survey

This repository provides the data and source code for the following study: Juliana C Taube, Zachary Susswein, Shweta Bansal. "Spatiotemporal trends in self-reported mask-wearing behavior in the United States: Analysis of a large cross-sectional survey." JMIR Public Health and Surveillance. https://doi.org/10.2196/42128. 

## How to use this resource
To rerun the models and reproduce the figures, start by opening `covid_masking.Rproj`. From here, run the files in numerical order, starting with `06GENERATE-BIAS.r` if you don't have access to the individual data (see Individual Data below). Code will not run correctly if files are sourced out of order. All necessary input data files should be in the repository, with the exception of a large CDC table with daily county mask mandate data which can be downloaded by users (see `data/`). 

## Estimates (`data/estimates/`)
Raw survey responses aggregated to the county-month level and model estimates at the county-month level for both the CTIS (`fb`) and ONM (`onm`) surveys.
* `fb_estimates_binomreg_FINAL.csv` contains self-reported mask-wearing estimates (`p_est`) from binomial regression model only
* `fb_estimates_binomreg_rake_FINAL.csv` contains self-reported mask-wearing estimates (`p_est`) from binomial regression model with raked and resampled individual responses 
* `fb_estimates_binomreg_rake_debias_FINAL.csv` contains self-reported mask-wearing estimates (`unbiased_p`) from binomial regression model with raked and resampled individual responses and an offset for bias
* `fb_estimates_binomreg_rake_debias_FINAL_excl.csv` contains self-reported mask-wearing estimates (`unbiased_p`) from binomial regression model with raked and resampled individual responses and an offset for bias where influential fips (pareto k $\ge$ 0.7 in the initial model) are excluded
* `fb_comm_estimates_binomreg_rake_FINAL.csv` contains community-reported mask-wearing estimates (`comm_p_est`) from binomial regression model with raked and resampled individual responses
* `fb_comm_estimates_binomreg_rake_FINAL_excl.csv` contains community-reported mask-wearing estimates (`comm_p_est`) from binomial regression model with raked and resampled individual responses where influential fips (pareto k $\ge$ 0.7 in the initial model) are excluded
* `onm_estimates_binomreg_FINAL.csv` contains self-reported mask-wearing estimates (`onm_est`) for the grocery store setting in the Outbreaks Near Me survey from the binomial regression model only
* `fb_processed.csv` contains raw self-reported mask-wearing data (`mask_prop_most`) aggregated to the county-month level
* `onm_processed.csv` contains raw self-reported mask-wearing data from the ONM survet (`mask_grocery_very_prop`) aggregated to the county-month level

## Data (`data/`)
Reference files that may be required to run the code, including fips and zipcode crosswalk files, urban/rural classifications, etc.
* `acs_target_data_for_fb_raking.csv` contains age, sex, education, and some race/ethnicity targets for raking from the 2019 ACS
* `archived_vaccination_data/data_master_county_20210823.csv` contains weekly vaccination coverage data for each fips code as recorded at https://vaccinetracking.us
* `Average_Household_Size_and_Population_Density_-_County.csv` contains population density for each fips code (from https://covid19.census.gov/datasets/USCensus::average-household-size-and-population-density-county/explore?location=4.945434%2C0.315550%2C1.99&showTable=true)
* `county_neighbors_fips.txt` contains a fips code neighbor list
* `COUNTY_ZIP_062021.xlsx` is the file for crosswalking fips code to zip code (from https://www.huduser.gov/portal/datasets/usps_crosswalk.html)
* `ZIP_COUNTY_062021.xlsx` is the file for crosswalking zip code to fips code (from https://www.huduser.gov/portal/datasets/usps_crosswalk.html)
* `fb_resampled_comm_data_from_raking_age_sex_educ_FINAL_corrected.csv` contains county-month level resampled data on community-reported masking following raking on age, sex, and education, and with correction of weights in nonconvergent county-months
* `fb_resampled_data_from_raking_age_sex_educ_FINAL_corrected.csv` contains county-month level resampled data on self-reported masking following raking on age, sex, and education, and with correction of weights in nonconvergent county-months
* `fb_resampled_vax_data_from_raking_age_sex_educ_FINAL_corrected_aprmay.csv` contains county-month level resampled data on vaccination following raking on age, sex, and education, and with correction of weights in nonconvergent county-months
* `NCHSURCodes2013.xlsx` contains urban-rural classifications for each fips code (from https://www.cdc.gov/nchs/data_access/urban_rural.htm)
* `nyt-us-rolling-avg.csv` contains COVID-19 case data from the New York Times (from https://github.com/nytimes/covid-19-data/blob/master/rolling-averages/us.csv)
* `state_fips.csv` contains each state's fips code
* `U.S._State_and_Territorial_Public_Mask_Mandates_From_April_8__2020_through_August_15__2021_by_State_by_Day.csv` contains mask mandate data at the state-level (from https://data.cdc.gov/Policy-Surveillance/U-S-State-and-Territorial-Public-Mask-Mandates-Fro/62d6-pm5i) 
* `U.S._State_and_Territorial_Public_Mask_Mandates_From_April_10__2020_through_August_15__2021_by_County_by_Day.csv` contains mask mandate data at the county-level (from https://data.cdc.gov/Policy-Surveillance/U-S-State-and-Territorial-Public-Mask-Mandates-Fro/62d6-pm5i)
  * This file is large and will be uploaded in the future. In the meantime, users can download the csv from the above link.
* `vax_data_for_time_series.RDS` contains cumulative proportions of all U.S. residents vaccinated by week (from `vaccinetracking.us`, code to produce these estimates is commented out in `13FIGURES.csv`)
* `worry_df_for_time_series.RDS` contains aggregated weekly data on proportion of respondents very or somewhat worried about severe COVID-19 from the CTIS 

## Code (`scripts/`)
Scripts to clean data, rake survey responses, resampled weighted survey responses, run binomial regression models, and reproduce figures. Scripts for analyzing individual responses are provided for reproducibility but will not run without the original individual-level data (see Individual Data section below). File names briefly describe the purpose of each script (where `COMM` stands for community-reported and `ONM` for the Outbreaks Near Me dataset).

## Individual Data
Individual CTIS survey responses cannot be shared by the authors, but researchers can visit https://cmu-delphi.github.io/delphi-epidata/symptom-survey/data-access.html if they would like to enter an agreement for data usage with CMU Delphi. Individual ONM responses also cannot be shared by the authors, but researchers can contact the OutbreaksNearMe team at Boston Children's Hospital and Momentive to inquire about access to these data. 

If users have access to the individual data, they will be able to run files 01 through 05.
