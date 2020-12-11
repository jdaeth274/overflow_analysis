###############################################################################
## Combining HF out data ######################################################
###############################################################################

altering_vals <- function(whole_data, new_data){
   
   p_weeks <- paste(new_data$p, new_data$week, sep = "-")
   uni_p_weeks <- unique(p_weeks)
   for( k in uni_p_weeks){
     current_p <- str_split_fixed(k,"-",2)[1]
     current_weeks <- str_split_fixed(k,"-",2)[2]
     adding_in_dat <- new_data[new_data$p == current_p & new_data$week == current_weeks,]
     present_or_not <- whole_data[whole_data$p == current_p & whole_data$week == current_weeks &
                                    whole_data$s == "C",]
     
     if(nrow(present_or_not) > 0){

       whole_data[whole_data$p == current_p & whole_data$week == current_weeks &
                    whole_data$s == "C", "pi_y" ] <- adding_in_dat$pi_y
     }else{
       whole_data <- bind_rows(whole_data, adding_in_dat)
     }
        
     
   }
    
   
  return(whole_data)
  
}

removing_20_22 <- function(forecasted_data){
  pdf_seq <- c("emergency_admissions","emergency_frailty","emergency_cc",
               "elective_pool","elective_median","elective_mean","elective_frail","elective_cc",
               "elective_bundles","emergency_bundles")
  
  
  for(k in 1:length(pdf_seq)){
    
    current_HF <- forecasted_data[[k]]
    
    if(any(grepl("_22_",current_HF$patient_group))){
      
      current_HF <- current_HF[-grep("_22_",current_HF$patient_group),]
    }
    if(any(grepl("_20_",current_HF$patient_group))){
      current_HF <- current_HF[-grep("_20_",current_HF$patient_group),]
    }
    
    current_HF <- current_HF[current_HF$date > "2020-03-01",]
    
    forecasted_data[[k]] <- current_HF
  
  }
  
  
  return(forecasted_data)  
  
}

in_hospital_maker <- function(hf_data, imp_data){

  template_out <- imp_data
  admi_groups <- c("N","E")
  s_groups <- c("G","C")
  template_out$y0 <- 0
  
  for(admissions in admi_groups){
    for(ward in s_groups){
      
      current_HF_dat <- hf_data[hf_data$a == admissions & hf_data$s == ward,]
      current_imp_dat <- imp_data[imp_data$a == admissions & imp_data$s == ward,]
      icd_groups <- unique(current_HF_dat$ICD)
      
      for(icd in icd_groups){
        current_hf_icd <- current_HF_dat[current_HF_dat$ICD == icd,]
        imp_rows <- current_imp_dat[grep(icd, current_imp_dat$p),]
        
        if(sum(imp_rows$y0) != 0){
        
          imp_split <- imp_rows$y0 / sum(imp_rows$y0)
        }else{
          imp_split <- c(0.2,0.4,0.4)
        }
        
        if(current_hf_icd$one == "<10"){
          hf_val <- 10L
        }else{
          hf_val <- as.integer(current_hf_icd$one)
        }
        
        hf_vals <- as.integer(hf_val * imp_split)
        
        ps_touse <- paste(rep(icd, 3), c("_AGE1","_AGE2","_AGE3"), sep = "")
        
        template_out[template_out$s == ward & 
                       template_out$a == admissions &
                       template_out$p %in% ps_touse, "y0"] <- hf_vals

        
      }
      
    }
  }
  
  return(template_out)
  
}

## forecasts 

hf_forecasts <- list.files("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/forecasts/THF_outputs_20200925/", full.names = TRUE)
hf_forecasts <- hf_forecasts[c(7,9,8,5,4,3,2,1,6,10)]
hf_fores <- lapply(hf_forecasts, read.csv, stringsAsFactors = FALSE)
hf_fores <- removing_20_22(hf_fores)


## regressions 
hf_piy_elec <- read.csv("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/transitions/THF_outputs_20200904_regressions/THF_outputs_20200904_regressions/elective_dat_table1.csv",
                        stringsAsFactors = FALSE)
hf_piy_emerg <- read.csv("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/transitions/THF_outputs_20200904_regressions/THF_outputs_20200904_regressions/emergency_dat_table1.csv",
                         stringsAsFactors = FALSE)

colnames(hf_piy_emerg)
colnames(hf_piy_elec)

hf_piy_elec <- hf_piy_elec[,c(1:7,9)]
hf_piy_emerg <- hf_piy_emerg[,c(1:4,6:7)]

## remove 22 & 20 ICDs

hf_piy_elec <- hf_piy_elec[-grep("ICD22", hf_piy_elec$p),]
hf_piy_elec <- hf_piy_elec[-grep("ICD20", hf_piy_elec$p),]
colnames(hf_piy_emerg)[6] <- "pi_y"

## Need to alter N_ICD12_AGE1 CC, N_ICD07_AGE3, N_ICD03_AGE3 10.5 days &
## N_ICD03_AGE1 no cc vals

## ICD12_AGE1 is C-G 1 in 1.5 in imp data 
## ICD07_AGE3 is C-C 1 in 1.5 in imp data
## ICD03_AGE3 is C-H 0, C-D 0.61, C-G 0.19, C-C 0.2 
## ICD03_AGE1 is C-G 1 for 0.5 

## ICD12 
template <- hf_piy_elec[5:8,]
icd_12 <- template
icd_07 <- template
icd_07_2 <- template
icd_03 <- template
icd_03_2 <- template
icd_031 <- template

icd_12$week <- 1.5
icd_12$p <- "ICD12_AGE1"
icd_12$pi_y <- c(0,1,0,0)




icd_07$week <- 1.5
icd_07$p <- "ICD07_AGE3"
icd_07$pi_y <- c(0,0,0,1)

icd_07_2$week <- 2.5
icd_07_2$p <- "ICD07_AGE3"
icd_07_2$pi_y <- c(0,0,1,0)

icd_03$week <- 1.5
icd_03$p <- "ICD03_AGE3"
icd_03$pi_y <- c(0,0.19,0.61,0.2)

icd_03_2$week <- 2.5
icd_03_2$p <- "ICD03_AGE3"
icd_03_2$pi_y <- c(0,0,1,0)

icd_031$week <- 0.5
icd_031$p <- "ICD03_AGE1"
icd_031$pi_y <- c(0,1,0,0)




missing_vals <- bind_rows(icd_12, icd_07, icd_07_2, icd_03, icd_03_2, icd_031)

hf_piy_elec <- altering_vals(hf_piy_elec, missing_vals)


## pools 

in_hosp_pool_df <- read.csv("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/pool_sizes/OneDrive_1_22-09-2020/in_hosp_patients.csv", 
                            stringsAsFactors = FALSE)
waiting_patients_df <- read.csv("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/pool_sizes/OneDrive_1_22-09-2020/waiting_patients.csv", 
                            stringsAsFactors = FALSE)

in_hosp_imp <- read.csv("D:/Dropbox/COVID19/Overflow/JOSH/imperial_in_hosp.csv",
                        stringsAsFactors = FALSE)

in_hosp_pool_df <- in_hosp_pool_df[-grep("ICD22",in_hosp_pool_df$ICD),]
waiting_patients_df <- waiting_patients_df[-grep(22,waiting_patients_df$icd),]
in_hosp_pool_df <- in_hosp_pool_df[-grep("ICD20",in_hosp_pool_df$ICD),]
waiting_patients_df <- waiting_patients_df[-grep(20,waiting_patients_df$icd),]
colnames(in_hosp_imp)[1] <- "a"

in_hosp_pool_df <- in_hospital_maker(in_hosp_pool_df, in_hosp_imp)


## COVID data 
covid_csv_loc <- "D:/Dropbox/COVID19/Overflow/JOSH/projections.csv"
covid_csv <- read.csv(file = covid_csv_loc,
                      stringsAsFactors = FALSE)
covid_probs <- read.csv(file = "D:/Dropbox/COVID19/Overflow/JOSH/6_8_data/Probabilities_covid_4_combined.csv",
                        stringsAsFactors = FALSE)
covid_no_cc <- read.csv(file = "D:/Dropbox/COVID19/Overflow/JOSH/6_8_data/Probabilities_covid_requiredCC_combined.csv",
                        stringsAsFactors = FALSE)

## Failure function 

failure_func_hf <- read.csv("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/transitions/elective_to_emergencies.csv",
                            stringsAsFactors = FALSE)

## lets put it all together 

reg_data <- list(hf_piy_elec, hf_piy_emerg, failure_func_hf)

forecasts_data_hf <- hf_fores


pool_data <- list(NULL,NULL,waiting_patients_df, NULL,NULL,in_hosp_pool_df)



excel_sheet <- "D:/Dropbox/COVID19/Overflow/Wolfram/HF_full_icd_28_09_2020.xlsx"

source("D:/Dropbox/COVID19/Overflow/Rscript/R/final_output_creator.R")
pi_y_df <- sum_up_function(reg_res = reg_data , time_series_data_res = pool_data,
                           time_series_forecasts = forecasts_data_hf, forecast_length = 78, COVID_preds = covid_csv,
                           out_sheet = excel_sheet, covid_probs = covid_probs, week_nums = seq(0.5,5,0.5), covid_no_cc = covid_no_cc)








