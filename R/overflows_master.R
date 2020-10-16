###############################################################################
## Master Overflows analysis script ###########################################
###############################################################################

usage <- "Overflows master script, outputs excel file with patient admissions and transitions probabilities"
input <- "Usage: R ./overflow_master.R <hes_data_loc> <transitions_data_write_loc> <forecast_start_date> <hes_cutoff_1> <hes_cutoff_2> <num_cores> <covid_csv_loc> <excel_loc>"
usage <- paste(usage, input, sep = "\n")

source(file = "./Rscript/R/packages.R")
source(file = "./Rscript/R/cohort_identification.R")
source(file = "./Rscript/R/transitions_coding.R")
source(file = "./Rscript/R/time_series_creator.R")
source(file = "./Rscript/R/time_series_forecast.R")
source(file = "./Rscript/R/regression_analyses.R")
source(file = "./Rscript/R/final_output_creator.R")
source(file = "./Rscript/R/costs_merging_HF.R")
## First step to load up HES dataset ##

input_args <- commandArgs(trailingOnly = TRUE)

if(length(input_args) != 10){
  cat(usage)
  cat("You have: ", length(input_args)," need: 10")
  quit(status = 1)
}

hes_data_loc <- input_args[1]
transitions_data_loc <- input_args[2]
forecast_start_date <- input_args[3]
cut_off_date_1 <- input_args[4]
cut_off_date_2 <- input_args[5]
core_nums <- input_args[6]
covid_csv_loc <- input_args[7]
excel_sheet <- input_args[8]
covid_prob_locs <- input_args[9]
costs_dir <- input_args[10]
covid_no_cc <- input_args[11]
cut_off_dates <- c(cut_off_date_1, cut_off_date_2)

cohort_allocation <- cohort_set_up(num_cores = core_nums,
                                   data_loc = input_args[1])
                                    
transitions_data <- hes_transitions(transitions_data)


## Write out transitions data ##
vroom_write(transitions_data,
            path = transitions_data_loc,
            delim = ",")

## regression analyses 


elective_no_trends_half <- regression_cluster_set_up(patient_group = "elective",
                                                hes_data = transitions_data, week_num = seq(1,5),
                                                failure_function_run = FALSE, half_week = TRUE)
elective_no_trends_whole <- regression_cluster_set_up(patient_group = "elective",
                                                     hes_data = transitions_data, week_num = seq(1,5),
                                                     failure_function_run = FALSE, half_week = FALSE)


emergency_no_trends_half <- regression_cluster_set_up(patient_group = "emergency",
                                                 hes_data = transitions_data, week_num = seq(1,5),
                                                 failure_function_run = FALSE, half_week = TRUE)
emergency_no_trends_whole <- regression_cluster_set_up(patient_group = "emergency",
                                                      hes_data = transitions_data, week_num = seq(1,5),
                                                      failure_function_run = FALSE, half_week = FALSE)

elective_dat <- list(dplyr::bind_rows(elective_no_trends_half[[1]], elective_no_trends_whole[[1]]),
                     dplyr::bind_rows(elective_no_trends_half[[2]], elective_no_trends_whole[[2]]))

emergency_dat <- list(dplyr::bind_rows(emergency_no_trends_half[[1]], emergency_no_trends_whole[[1]]),
                     dplyr::bind_rows(emergency_no_trends_half[[2]], emergency_no_trends_whole[[2]]))


failure_func <- regression_cluster_set_up(patient_group = "failure",
                                          hes_data = transitions_data, forecast_length = 52,
                                          forecast_start = forecast_start_date,
                                          month_trend = FALSE, time_trend = FALSE, week_num = seq(1,5),
                                          failure_function_run = TRUE, half_week = TRUE)
regression_results <- list(elective_dat, emergency_dat, failure_func[[3]])

## time_series ##
time_series_data <- time_series_creator(hes_data = transitions_data, num_cores = core_nums, forecast_date = forecast_start_date,
                                        emergency_run = TRUE, elective_ts = TRUE, forecast_cutoff = forecast_start_date)
source(file = "./Rscript/R/time_series_forecast.R")
times_series_forecasts <- running_forecasts(total_cohort_data = time_series_data, train_date = as.Date(forecast_start_date),
                                       forecast_period = 78,  base_dir = "./JOSH/HF_out_data/imperial_forecasts/seasonsal_prop_bundles_",
                                       run_admis = FALSE, forecast_admis = FALSE, forecast_WT = FALSE, forecast_frail = FALSE,
                                       cutoff_dates = cutoff_dates, run_diags = FALSE, forecast_cc = FALSE)

## costs ##

costs_res <- costs_function(transitions_data, costs_directory = costs_dir, FY = 1213, output_dir = "./")


## Sum up file 
covid_csv <- read.csv(file = covid_csv_loc,
                      stringsAsFactors = FALSE)
covid_probs <- read.csv(file = "D:/Dropbox/COVID19/Overflow/JOSH/6_8_data/Probabilities_covid_4_combined.csv",
                        stringsAsFactors = FALSE)
covid_no_cc <- read.csv(file = "D:/Dropbox/COVID19/Overflow/JOSH/6_8_data/Probabilities_covid_requiredCC_combined.csv",
                        stringsAsFactors = FALSE)


pi_y_df <- sum_up_function(reg_res = regression_results , time_series_data_res = time_series_data,
                time_series_forecasts = times_series_forecasts, forecast_length = 78, COVID_preds = covid_csv,
                out_sheet = excel_sheet, covid_probs = covid_probs, week_nums = seq(0.5,5,0.5), covid_no_cc = covid_no_cc)




