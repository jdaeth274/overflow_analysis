###############################################################################
## Master Overflows analysis script ###########################################
###############################################################################

source(file = "src/packages.R")
source(here::here("src", "file_paths.R"))
source(here::here("src", "prep_functions.R"))

source(file = "src/cohort_identification.R")
source(file = "src/transitions_coding.R")
source(file = "src/time_series_creator.R")
source(file = "src/time_series_forecast.R")
source(file = "src/regression_analyses.R")
source(file = "src/costs_merging_HF.R")

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

#forecast_start_date_vec <- str_c(as.numeric(HES_years[1:5]), "-03-01")

transitions_data_loc <- str_c(here::here("transitions_data", "APC_CC_"), HES_years_condensed, "_transitions.csv")


cohort_transitions_pipeline <- function(hes_data_file, transitions_data_loc, core_nums){
  
  hes_data <- readRDS(hes_data_file) 
  
  hes_data <- hes_data %>% 
    rename(GA_LoS = "GA_LOS",
           cc_start_flg = "cc_start_flag",
           cc_dis_flg = "cc_dis_flag")
  
  hes_data <- hes_data %>% 
    mutate(cc_start_flg = replace_na(cc_start_flg, 0),
           cc_dis_flg = replace_na(cc_dis_flg, 0))
  
  hes_data <- hes_data %>% 
    mutate(GA_LoS = ifelse(!is.na(cc_LoS), Total_LoS - cc_LoS, Total_LoS))
  
  
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
transitions_data_2016 <- transitions_data_2016 %>% 
  mutate(gaccga_death = as.numeric(gaccga_death))
transitions_data <- bind_rows(transitions_data, transitions_data_2016)

rm(transitions_data_2016)

transitions_data_2017 <- read_csv(transitions_data_loc[3])
transitions_data_2017 <- transitions_data_2017 %>% 
  mutate(gaccga_death = as.numeric(gaccga_death))
transitions_data <- bind_rows(transitions_data, transitions_data_2017)
rm(transitions_data_2017)

transitions_data_2018 <- read_csv(transitions_data_loc[4])
transitions_data_2018 <- transitions_data_2018 %>% 
  mutate(gaccga_death = as.numeric(gaccga_death))
transitions_data <- bind_rows(transitions_data, transitions_data_2018)
rm(transitions_data_2018)

transitions_data_201920 <- read_csv(transitions_data_loc[5])
transitions_data_201920 <- transitions_data_201920 %>% 
  mutate(gaccga_death = as.numeric(gaccga_death))
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

transitions_data <- transitions_data %>% 
  rename(SUSHRGep = "SUSHRG") 

#  floor_date is used to make sure this works for Feb 29
transitions_data <- transitions_data %>% 
  mutate(EpiEnd_FY = as.numeric(if_else(month(epiend) %in% c(4:12), 
                                        paste0(format(epiend, format = "%y"),
                                               format(floor_date(epiend, "month") + years(1), format = "%y" )),
                                        paste0(format(floor_date(epiend, "month") - years(1), format = "%y" ),
                                               format(epiend, format = "%y")))))


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

# Exclude electives after March 2019
time_series_data[[2]] <- time_series_data[[2]] %>%
  tidylog::filter(date <= ymd("2019-03-31"))

time_series_data[[4]] <- time_series_data[[4]] %>%
  tidylog::filter(!(rttstart_YYYY == 2020 | (rttstart_YYYY == 2019 & rttstart_week > 13)))

times_series_forecasts <- running_forecasts(total_cohort_data = time_series_data,
                                            train_date = as.Date("2020-03-02"),
                                            forecast_period = 53,
                                            base_dir = "time_series_forecasts/",
                                            run_admis = TRUE, forecast_admis = TRUE,
                                            forecast_WT = TRUE, forecast_frail = TRUE,
                                            cutoff_dates = cutoff_dates)

saveRDS(times_series_forecasts, here::here("time_series_forecasts", "time_series_forecasts.Rds"))



# -------------------------------------------------------------------------

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


elective_dat_table1 <- dplyr::bind_rows(elective_no_trends_half[[1]], elective_no_trends_whole[[1]])

elective_dat_table2 <- dplyr::bind_rows(elective_no_trends_half[[2]], elective_no_trends_whole[[2]])

saveRDS(elective_dat_table1, here::here("regression_results", "elective_dat_table1.Rds"))
saveRDS(elective_dat_table2, here::here("regression_results", "elective_dat_table2.Rds"))


elective_dat_table1_censored <- elective_dat_table1 %>% 
  tidylog::filter(!(group_size %in% c(1:5)))

write_csv(elective_dat_table1_censored, 
          here::here("regression_results", "elective_dat_table1.csv"))

elective_dat_table2_censored <- elective_dat_table2 %>% 
  mutate(p = paste0("ICD", str_pad(ICD, width = 2, side = "left" ,pad = "0"), "_AGE", age)) %>% 
  semi_join(elective_dat_table1_censored, by = c("week", "p"))

write_csv(elective_dat_table2_censored, here::here("regression_results", "elective_dat_table2.csv"))


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


emergency_dat_table1 <- dplyr::bind_rows(emergency_no_trends_half[[1]], emergency_no_trends_whole[[1]])

saveRDS(emergency_dat_table1, here::here("regression_results", "emergency_dat_table1.Rds"))


emergency_dat_table1_censored <- emergency_dat_table1 %>% 
  tidylog::filter(!(group_size %in% c(1:5))) %>% 
  group_by(a, p, s, week) %>% 
  mutate(group_size_prop = group_size * pi_y) %>% 
  mutate(pi_y_censored = case_when(all(group_size_prop == 0 | group_size_prop >= 10) ~ pi_y, 
                                   group_size_prop %in% c(1:9) ~ 0),
         group_size_prop_censored = case_when(all(group_size_prop == 0 | group_size_prop >= 10) ~ group_size_prop, 
                                              group_size_prop %in% c(1:9) ~ 0),
         censored = case_when(all(group_size_prop == 0 | group_size_prop >= 10) ~ 0,
                              group_size_prop %in% c(1:9) ~ 1),
         rank = rank(group_size_prop, ties.method = "first"),
         second = if_else(any(group_size_prop %in% c(1:9)) &
                            rank == min(rank[group_size_prop>=10]), 1, 0),
         group_size_diff = if_else(second == 1, 
                                   sum(group_size_prop_censored[censored == 1], na.rm = TRUE) + sum(group_size_prop[censored == 1], na.rm = TRUE),
                                   0),
         group_size_prop_censored = if_else(second == 1, group_size_prop + group_size_diff, group_size_prop_censored),
         pi_y_censored = if_else(second == 1, group_size_prop_censored/group_size, pi_y_censored),
         censored = if_else(second == 1, 1, censored),
         group_size_prop_censored = if_else(is.na(group_size_prop_censored), group_size_prop, group_size_prop_censored),
         pi_y_censored = if_else(is.na(pi_y_censored), pi_y, pi_y_censored),
         censored = if_else(is.na(censored), 0, censored),
         censored_per_group = sum(censored)) %>% 
  tidylog::filter(is.na(group_size) | max(group_size_prop) >= 10)

emergency_dat_table1_censored <- emergency_dat_table1_censored %>% 
  mutate(group_size_censored = sum(group_size_prop_censored),
         pi_y_sum_censored = sum(pi_y_censored))

check_group_size <- emergency_dat_table1_censored$group_size == emergency_dat_table1_censored$group_size_censored
any(check_group_size) == FALSE

check_prop_sum <- emergency_dat_table1_censored$pi_y_sum_censored == 1
any(check_prop_sum) == FALSE

emergency_dat_table1_censored <- emergency_dat_table1_censored %>% 
  select(-group_size_diff, - group_size_prop, -pi_y, - coeff, 
         -variance, - second, - rank, - group_size_censored,
         - pi_y_sum_censored)

write_csv(emergency_dat_table1_censored, 
          here::here("regression_results", "emergency_dat_table1.csv"))



failure_func_output <- regression_cluster_set_up(patient_group = "failure",
                                          hes_data = transitions_data, forecast_length = 52,
                                          forecast_start = "2019-03-05",
                                          month_trend = FALSE, time_trend = FALSE, week_num = seq(1,5),
                                          failure_function_run = TRUE, half_week = TRUE)

write_csv(failure_func_output[[3]], here::here("survival", "elective_to_emergencies.csv"))



# Costing analysis --------------------------------------------------------



transitions_data <- transitions_data %>% 
  rename(procode3 = "Procode") 

c("1516", "1617", "1718", "1819") %>% 
  walk(~costs_function(transitions_data_whole = transitions_data, 
                costs_directory = "../Reference costs/", 
                output_dir = here::here("costing_results"), fyear = .x))

full_mergedcosts <- readRDS(here::here("costing_results", "mergedcosts_1516.Rds"))

costs_1617 <- readRDS(here::here("costing_results", "mergedcosts_1617.Rds"))
full_mergedcosts <- rbind(full_mergedcosts, costs_1617)
rm(costs_1617)

costs_1718 <- readRDS(here::here("costing_results", "mergedcosts_1718.Rds"))
full_mergedcosts <- rbind(full_mergedcosts, costs_1718)
rm(costs_1718)

costs_1819 <- readRDS(here::here("costing_results", "mergedcosts_1617.Rds"))
full_mergedcosts <- rbind(full_mergedcosts, costs_1819)
rm(costs_1819)

saveRDS(full_mergedcosts, here::here("costing_results", "full_mergedcosts.Rds"))

calculate_unitcosts(full_mergedcosts, 
                    costs_directory = "../Reference costs/", 
                    output_dir = here::here("costing_results"))
  