  ## comparing the pool outputs between HF and out HES extract ##

plotting_function <- function(hf_data, imp_data, file_name = "78_weeks"){

  pdf_seq <- c("emergency_admissions","emergency_frailty","emergency_cc",
               "elective_pool","elective_median","elective_mean","elective_frail","elective_cc",
               "elective_bundles","emergency_bundles")

  
  for(k in 1:length(pdf_seq)){

    print(pdf_seq[k])
    if(k == 9)
      next

    
    current_HF <- hf_data[[k]]
    current_imp <- imp_data[[k]]
    
    current_HF$data <- "HF"
    current_imp$data <- "Imperial"
    tot_patient_groups <- unique(c(current_HF$patient_group, current_imp$patient_group))
    if(any(grepl("_22_",tot_patient_groups))){
      tot_patient_groups <- tot_patient_groups[-grep("_22_",tot_patient_groups)]
    }
    current_pdf <- paste("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/",file_name ,pdf_seq[k],".pdf",sep = "")
    pdf(file = current_pdf, paper = "a4r", width = 10, height = 7)
    
    for(j in 1:length(tot_patient_groups)){
      print(j)
      
      current_title <- paste(pdf_seq[k], tot_patient_groups[j])
      
      HF_patient_group_dat <- current_HF[current_HF$patient_group == tot_patient_groups[j],]# &
                                           #current_HF$date > "2020-03-01",]
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



hf_forecasts <- list.files("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/forecasts/forecast_2020_09_21/forecasts/forecasts/", full.names = TRUE)
hf_forecasts <- hf_forecasts[c(7,9,8,5,4,3,2,1,6,10)]
hf_fores <- lapply(hf_forecasts, read.csv, stringsAsFactors = FALSE)
imperial_forecasts <- list.files("D:/Dropbox/COVID19/Overflow/JOSH/HF_out_data/imperial_forecasts/", full.names = TRUE)
imperial_forecasts <- imperial_forecasts[c(7,9,8,5,4,3,2,1,6,10)]
imp_fores <- lapply(imperial_forecasts, read.csv, stringsAsFactors = FALSE)

plotting_function(hf_fores, imp_fores)

