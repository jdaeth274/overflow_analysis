  ## comparing the pool outputs between HF and out HES extract ##

plotting_function <- function(hf_data, imp_data, file_name = "78_weeks", input_file){

  pdf_seq <- c("emergency_admissions","emergency_frailty","emergency_cc",
               "elective_pool","elective_median","elective_mean","elective_frail","elective_cc",
               "elective_bundles","emergency_bundles")

  
  for(k in 1:length(pdf_seq)){

    print(pdf_seq[k])
    
    
    
    current_HF <- hf_data[[k]]
    current_imp <- imp_data[[k]]
    
    current_HF$data <- "HF"
    current_imp$data <- "Imperial"
    tot_patient_groups <- unique(c(current_HF$patient_group, current_imp$patient_group))
    if(any(grepl("_22_",tot_patient_groups))){
      tot_patient_groups <- tot_patient_groups[-grep("_22_",tot_patient_groups)]
    }
    if(any(grepl("_20_",tot_patient_groups))){
      tot_patient_groups <- tot_patient_groups[-grep("_20_",tot_patient_groups)]
    }
    if(any(grepl("_15_3",tot_patient_groups))){
      tot_patient_groups <- tot_patient_groups[-grep("_15_3",tot_patient_groups)]
    }
    current_pdf <- paste(input_file,file_name ,pdf_seq[k],".pdf",sep = "")
    pdf(file = current_pdf, paper = "a4r", width = 10, height = 7)
    
    for(j in 1:length(tot_patient_groups)){
      print(j)
      
      current_title <- paste(pdf_seq[k], tot_patient_groups[j])
      
      HF_patient_group_dat <- current_HF[current_HF$patient_group == tot_patient_groups[j] &
                                           current_HF$date > "2020-03-01",]
      imp_patient_group_dat <- current_imp[current_imp$patient_group == tot_patient_groups[j],]
      HF_patient_group_dat$time <- seq(1, nrow(HF_patient_group_dat))
      imp_patient_group_dat$time <- seq(1, nrow(HF_patient_group_dat))
      
      tot_dat_pg <- dplyr::bind_rows(HF_patient_group_dat, imp_patient_group_dat)
      
      current_plot <- ggplot(data = tot_dat_pg, aes(x = time, y = median, group = data)) +
        geom_line(aes(color = data)) + geom_ribbon(aes(ymin = lower, ymax = upper, fill = data, colour = data),
                                                   alpha = 0.25) +
        theme_bw() + ggtitle(current_title) + xlab("Time") + ylab(pdf_seq[k])
      
      print(current_plot)
      
      
    }
    
    dev.off()
    
    
    
  }
  
  
  
}



hf_forecasts <- list.files("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/forecasts/THF_outputs_20200924/THF_outputs_20200924/", full.names = TRUE)
hf_forecasts <- hf_forecasts[c(7,9,8,5,4,3,2,1,6,10)]
hf_fores <- lapply(hf_forecasts, read.csv, stringsAsFactors = FALSE)
imperial_forecasts <- list.files("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/imperial_forecasts/", full.names = TRUE)
imperial_forecasts <- imperial_forecasts[c(7,9,8,5,4,3,2,1,6,10)]
imp_fores <- lapply(imperial_forecasts, read.csv, stringsAsFactors = FALSE)

plotting_function(hf_fores, imp_fores, input_file = "~/Dropbox/Overflow/JOSH/HF_out_data/25-09-2020_HF")

## checking the bundling 


elec_hf_bundles <- hf_fores[[9]]
elec_hf_bundles <- elec_hf_bundles[elec_hf_bundles$date > "2020-03-01",]

age_1 <- elec_hf_bundles[elec_hf_bundles$age == 1,]
age_1[1:78,"median"] + age_1[79:156,"median"] + age_1[157:234,"median"] + age_1[235:312,"median"] + age_1[313:390,"median"]
age_2 <- elec_hf_bundles[elec_hf_bundles$age == 2,]
age_2[1:78,"median"] + age_2[79:156,"median"] + age_2[157:234,"median"] + age_2[235:312,"median"] + age_2[313:390,"median"]
age_3 <- elec_hf_bundles[elec_hf_bundles$age == 3,]
age_3[1:78,"median"] + age_3[79:156,"median"] + age_3[157:234,"median"] + age_3[235:312,"median"]



imperial_elec_testing <- read.csv("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/imperial_forecasts/seasonsal_prop_bundles_cohort_1_prop_bundle.csv",
                                  stringsAsFactors = FALSE)


imperial_elec_testing
age_1 <- imperial_elec_testing[imperial_elec_testing$age == 1,]
age_1[1:78,"median"] + age_1[79:156,"median"] + age_1[157:234,"median"] + age_1[235:312,"median"] + age_1[313:390,"median"]
age_2 <- imperial_elec_testing[imperial_elec_testing$age == 2,]
age_2[1:78,"median"] + age_2[79:156,"median"] + age_2[157:234,"median"] + age_2[235:312,"median"] + age_2[313:390,"median"]
age_3 <- imperial_elec_testing[imperial_elec_testing$age == 3,]
age_3[1:78,"median"] + age_3[79:156,"median"] + age_3[157:234,"median"] + age_3[235:312,"median"]



