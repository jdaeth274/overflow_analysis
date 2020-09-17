## comparing the pool outputs between HF and out HES extract ##

plotting_function <- function(hf_data, imp_data){
  
  pdf_seq <- c("emergency_admissions","emergency_frailty","emergency_cc",
               "elective_pool","elective_median","elective_mean","elective_frail","elective_cc",
               "elective_bundles","emergency_bundles")
  
  for(k in 1:length(pdf_seq)){
    current_HF <- hf_data[[k]]
    current_imp <- imp_data[[k]]
    
    current_HF$data <- "HF"
    current_imp$data <- "Imperial"
    current_HF$time <- rep(seq(0, 52), length(unique(current_HF$patient_group)))
    current_imp$time <- rep(seq(0,77), length(unique(current_imp$patient_group)))
    
    tot_data <- dplyr::bind_rows(current_HF, current_imp)
    
    tot_patient_groups <- unique(tot_data$patient_groups)
    current_pdf <- paste("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/",pdf_seq[k],".pdf",sep = "")
    
    
    
    
  }
  
  
  
}



hf_forecasts <- list.files("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/forecasts/THF_outputs_20200901_2/THF_outputs_20200901_2/", full.names = TRUE)
hf_forecasts <- hf_forecasts[c(7,9,8,5,4,3,2,1,6,10)]
hf_fores <- lapply(hf_forecasts, read.csv, stringsAsFactors = FALSE)



