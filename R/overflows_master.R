###############################################################################
## Master Overflows analysis script ###########################################
###############################################################################


# NOTE: for functions called below (time series creation, forecast and regression)
# only the parameter combinations below have been tested

source(file = "src/packages.R")
source(here::here("src", "file_paths.R"))
source(here::here("src", "prep_functions.R"))

source(file = "src/cohort_identification.R")
source(file = "src/transitions_coding.R")
source(file = "src/time_series_creator.R")
source(file = "src/time_series_forecast.R")
source(file = "src/regression_analyses.R")

# Parameters --------------------------------------------------------------

HES_years <- c("2015", "2016", "2017", "2018", "2019", "2020")
HES_years_condensed <- c("2015", "2016", "2017", "2018", "201920")

# Number of cores to use for cohort allocation and time series creation
core_nums <- 12


# Prepare HES data  -------------------------------------------------------
# functions in data_prep.R

# Connect to db with clean HES data and extract an process one year at a time
db <- DBI::dbConnect(RSQLite::SQLite(), database_path_VM)

HES_years %>% 
  walk(prepare_HES_year, db = db, output_folder = "processed_data")
     
dbDisconnect(db)
    

# Combine 2019 and 2020 data ----------------------------------------------
# 2020 isn't a full year

APC_CC_2019 <-  readRDS(here::here("processed_data", "APC_CC_2019.Rds"))
APC_CC_2020 <-  readRDS(here::here("processed_data", "APC_CC_2020.Rds"))

APC_CC_201920 <-  APC_CC_2019 %>% 
  bind_rows(APC_CC_2020)

saveRDS(APC_CC_201920, here::here("processed_data", "APC_CC_201920.Rds"))

rm(APC_CC_2019)
rm(APC_CC_2020)
rm(APC_CC_201920)

# Create cohorts and transitions ------------------------------------------------------
# iterate over HES years

hes_data_file_vec <- str_c(here::here("processed_data", "APC_CC_"), HES_years_condensed, ".Rds")

forecast_start_date_vec <- str_c(as.numeric(HES_years[1:5]), "-03-01")

transitions_data_loc <- str_c(here::here("transitions_data", "APC_CC_"), HES_years_condensed, "_transitions.csv")


cohort_transitions_pipeline <- function(hes_data_file, transitions_data_loc, core_nums){
  
  hes_data <- readRDS(hes_data_file) 
  
  cohort_allocation <- cohort_set_up(num_cores = core_nums,
                                     data_loc = hes_data)
  
  transitions_data <- hes_transitions(cohort_allocation)
  
  vroom_write(transitions_data,
              path = transitions_data_loc,
              delim = ",")
}

walk2(hes_data_file_vec, transitions_data_loc,
     ~cohort_transitions_pipeline(hes_data_file = .x, transitions_data_loc = .y, core_nums = core_nums))


# Combining transitions data ----------------------------------------------
# not the most elegant way to do this

transitions_data <- read_csv(transitions_data_loc[1])

transitions_data_2016 <- read_csv(transitions_data_loc[2])
transitions_data <- bind_rows(transitions_data, transitions_data_2016)
rm(transitions_data_2016)

transitions_data_2017 <- read_csv(transitions_data_loc[3])
transitions_data <- bind_rows(transitions_data, transitions_data_2017)
rm(transitions_data_2017)

transitions_data_2018 <- read_csv(transitions_data_loc[4])
transitions_data <- bind_rows(transitions_data, transitions_data_2018)
rm(transitions_data_2018)

transitions_data_201920 <- read_csv(transitions_data_loc[5])
transitions_data <- bind_rows(transitions_data, transitions_data_201920)
rm(transitions_data_201920)

# Ironing out some last bits
transitions_data <- transitions_data %>% 
  filter(is.na(WaitingTime) | WaitingTime >= 0)

transitions_data <- transitions_data %>% 
  filter(!is.na(agegrp_v3))

transitions_data <- transitions_data %>% 
  mutate(rttstart_week = as.numeric(rttstart_week),
         epistart_week = as.numeric(epistart_week)) 

saveRDS(transitions_data, here::here("transitions_data", "transitions_combined.Rds"))

# Time series creation ----------------------------------------------------

time_series_data <- time_series_creator(hes_data = transitions_data, 
                                        num_cores = 8, 
                                        forecast_date = "2019-03-05",
                                        emergency_run = TRUE, 
                                        elective_ts = TRUE, 
                                        forecast_cutoff = "2019-03-05",
                                        waiting_pool = TRUE)

saveRDS(time_series_data, here::here("time_series_forecasts", "time_series_data.Rds"))

# Save outputs: patients in hospital and waiting listts
write_csv(time_series_data[[3]], here::here("time_series_forecasts", "waiting_patients.csv"))

time_series_data[[6]] %>%  
  separate(p, into = c("ICD", "age"), sep = "_") %>% 
  group_by(a, ICD, s) %>% 
  summarise(one = sum(one, na.rm = TRUE)) %>% 
  write_csv(., here::here("time_series_forecasts", "in_hosp_patients.csv"))


# Time series forecast ----------------------------------------------------

cut_off_date_1 <- "2015-01-01"
cut_off_date_2 <- "2020-02-29"
cutoff_dates <- c(cut_off_date_1, cut_off_date_2)


times_series_forecasts <- running_forecasts(total_cohort_data = time_series_data, 
                                            train_date = as.Date("2020-03-02"),
                                            forecast_period = 53,  
                                            base_dir = "time_series_forecasts/",
                                            run_admis = TRUE, forecast_admis = TRUE, 
                                            forecast_WT = TRUE, forecast_frail = TRUE,
                                            cutoff_dates = cutoff_dates)

saveRDS(times_series_forecasts, here::here("time_series_forecasts", "time_series_forecasts.Rds"))



# -------------------------------------------------------------------------

transitions_data <- transitions_data %>% 
  rename(GA_LoS = "GA_LOS",
         cc_start_flg = "cc_start_flag")

elective_no_trends_half <- regression_cluster_set_up(patient_group = "elective",
                                                     forecast_start = "2019-03-02",
                                                     start_date = cut_off_date_1,
                                                     hes_data = transitions_data, week_num = seq(1,5),
                                                     failure_function_run = FALSE, half_week = TRUE)

elective_no_trends_whole <- regression_cluster_set_up(patient_group = "elective",
                                                      forecast_start = "2019-03-02",
                                                      start_date = cut_off_date_1,
                                                      hes_data = transitions_data, week_num = seq(1,5),
                                                      failure_function_run = FALSE, half_week = FALSE)

elective_dat <- list(dplyr::bind_rows(elective_no_trends_half[[1]], elective_no_trends_whole[[1]]),
                     dplyr::bind_rows(elective_no_trends_half[[2]], elective_no_trends_whole[[2]]))

write_csv(elective_dat[[1]], here::here("regression_results", "elective_dat_table1.csv"))
write_csv(elective_dat[[2]], here::here("regression_results", "elective_dat_table2.csv"))


emergency_no_trends_half <- regression_cluster_set_up(patient_group = "emergency",
                                                      forecast_start = "2019-03-02",
                                                      start_date = cut_off_date_1,
                                                      hes_data = transitions_data, week_num = seq(1,5),
                                                      failure_function_run = FALSE, half_week = TRUE)

emergency_no_trends_whole <- regression_cluster_set_up(patient_group = "emergency",
                                                       forecast_start = "2019-03-02",
                                                       start_date = cut_off_date_1,
                                                       hes_data = transitions_data, week_num = seq(1,5),
                                                       failure_function_run = FALSE, half_week = FALSE)


emergency_dat <- list(dplyr::bind_rows(emergency_no_trends_half[[1]], emergency_no_trends_whole[[1]]),
                      dplyr::bind_rows(emergency_no_trends_half[[2]], emergency_no_trends_whole[[2]]))

write_csv(emergency_dat[[1]], here::here("regression_results", "emergency_dat_table1.csv"))
write_csv(emergency_dat[[2]], here::here("regression_results", "emergency_dat_table2.csv"))

failure_func <- regression_cluster_set_up(patient_group = "failure",
                                          hes_data = transitions_data, forecast_length = 52,
                                          forecast_start = "2019-03-02",
                                          month_trend = FALSE, time_trend = FALSE, week_num = seq(1,5),
                                          failure_function_run = TRUE, half_week = TRUE)

write_csv(failure_func[[3]], here::here("regression_results", "elective_to_emergencies.csv"))

regression_results <- list(elective_dat, emergency_dat, failure_func[[3]])

