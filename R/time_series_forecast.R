library(tidyverse)
library(ggplot2)
library(lubridate)
library(KFAS)
library(ggpubr)
library(forecast)
setwd("~/Dropbox/covid19/overflow_analysis/")

## new cohort data ##

cohort_1_ts <- read.csv("~/Dropbox/covid19/overflow_analysis/time_series_data/TS_cohort1_ICDAge.csv",
                        stringsAsFactors = FALSE)
cohort_2_ts <- read.csv("~/Dropbox/covid19/overflow_analysis/time_series_data/TS_cohorts2_ICDAge.csv",
                        stringsAsFactors = FALSE)
cohort_3_ts <- read.csv("~/Dropbox/covid19/overflow_analysis/time_series_data/TS_cohorts3_ICDAge.csv",
                        stringsAsFactors = FALSE)
cohort_3_ts_admissions <- read.csv("~/Dropbox/covid19/overflow_analysis/time_series_data/TS_cohorts3_ICDAge_admidate.csv",
                                   stringsAsFactors = FALSE)

diagnostics_function <- function(train_data, metric, out, title_name){
  
  train_data$one.step.res.std<-as.numeric(rstandard(out,type='recursive'))
  
  admis <- ggplot(train_data)+aes(x=date,y=train_data[,2], colour = "blue")+geom_line()+
    ggtitle('Admissions')+
    theme_bw(base_size = 12, base_family = "Helvetica")+
    theme(legend.position="none")+ylab("")+xlab("")
  
  
  recursive_resid <- ggplot(train_data)+aes(x=date,y=one.step.res.std)+geom_line()+
    ggtitle('Recursive residuals')+
    theme_bw(base_size = 12, base_family = "Helvetica")+
    theme(legend.position="none")+ylab("")+xlab("")
  
  qq_ggplot <- ggplot(train_data, aes(sample = one.step.res.std))
  gg_qq_plot <- qq_ggplot + stat_qq() + stat_qq_line() + ggtitle("QQ plot")
  
  
  hist_data <- cbind.data.frame(train_data$one.step.res.std, dnorm(train_data$one.step.res.std, mean = 0,
                                                                   sd = 1))
  colnames(hist_data) <- c("Residuals_recursive", "normal_density")
  
  hist_plot <- ggplot(data = hist_data) + geom_histogram(aes( x = Residuals_recursive,
                                                              y = ..density..),
                                                         binwidth = 0.5, fill = "lightgrey",
                                                         color = "black") +
    geom_line(aes(x = Residuals_recursive, y = normal_density), colour = "midnightblue")
  
  gg_acf_plot <- ggAcf(train_data$one.step.res.std,na.action=na.pass, type = "correlation")
  
  train_data$state.res.std<-as.numeric(rstandard(out,type='state')[,1])
  
  smooth_residuals <- ggplot(train_data)+aes(x=date,y=state.res.std)+geom_line()+
    ggtitle('Standardised smoothed (auxiliary) residuals')+
    theme_bw(base_size = 12, base_family = "Helvetica")+
    theme(legend.position="none")+ylab("")+xlab("")
  
  train_data$trend <-as.numeric(signal(out, states = c('trend'), filtered=FALSE )$signal)
  
  trend <- ggplot(train_data)+aes(x=date,y=trend, colour = 'red')+geom_line()+
    geom_line(aes(x=date,y=train_data[,2],colour='blue'))+
    ggtitle('Trend')+
    theme_bw(base_size = 12, base_family = "Helvetica")+
    theme(legend.position="none")+ylab("")+xlab("")
  
  train_data$seasonality <-as.numeric(signal(out, states = c('season'), filtered=FALSE )$signal)
  
  seasonal <- ggplot(train_data)+aes(x=date,y=seasonality)+geom_line()+
    geom_line(aes(x=date,y=train_data[,2],colour='blue'))+
    geom_line(aes(x=date,y=trend,colour='red'))+
    ggtitle('Seasonal pattern')+
    theme_bw(base_size = 12, base_family = "Helvetica")+
    theme(legend.position="none")+ylab("")+xlab("")
  
  train_data$predicted <-as.numeric(signal(out, states = c('season','trend'), filtered=FALSE )$signal)
  predicted <- ggplot(train_data)+aes(x=date,y=predicted)+ geom_jitter()+
    geom_line(aes(x=date,y=train_data[,2],colour='Red'))+
    ggtitle('Predicted ')+
    theme_bw(base_size = 12, base_family = "Helvetica")+
    theme(legend.position="none")+ylab("")+xlab("")
  
  multi_state <- ggpubr::ggarrange(admis,hist_plot, gg_qq_plot,
                                   recursive_resid,smooth_residuals, trend,
                                     seasonal,predicted, gg_acf_plot, ncol = 3, nrow = 3)
  multi_state <- annotate_figure(multi_state,
                                 top = text_grob(title_name, color = "Midnightblue", face = "bold", size = 12))
  
  print(multi_state)
  
}



forecast_function <- function(ts_data, results_name = "results.pdf",
                              forecast_frail = FALSE, forecast_wait = FALSE,
                              forecast_admissions = TRUE,
                              admit_type = "Emergency",
                              diagnostics_only = FALSE,
                              wait_time_col = "p50_WTwkyy",
                              cohort3 = FALSE, train_date = "2012-01-01",
                              only_ICD = NULL){
  

  # join the year and week variable to a date. Use UK  convention that week starts on a Monday.
  if(cohort3 == TRUE){
    ts_data$date<-as.POSIXct( paste( 1, ts_data$rttstart_week, ts_data$rttstart_YYYY, sep = "-" ), format = "%u-%U-%Y",locale = "UK" ) 
    ts_data$date <- as.Date(ts_data$date)
    
  }else{
    ts_data$date<-as.POSIXct( paste( 1, ts_data$admidate_week, ts_data$admidate_YYYY, sep = "-" ), format = "%u-%U-%Y",locale = "UK" ) 
    ts_data$date <- as.Date(ts_data$date)
  }
  ## remove date na rows ##
  na_dates <- which(is.na(ts_data$date))
  if(length(na_dates) > 0)
    ts_data <- ts_data[-na_dates,]
  
  
  quad_y <- grep("YYYY", colnames(ts_data))
  
  
  time_start <- Sys.time()

  pdf(file = results_name, paper = "A4r", width = 10, height = 7)
  
  group_num <- length(unique(ts_data$agegrp_v3))*length(unique(ts_data$ICD))
  current_group <- 1
  
  df_preds <- data.frame(patient_group = numeric(),
                                lower=numeric(),
                                median=numeric(),
                                upper=numeric(),
                                date=as.Date(character()))
  df_wait <- data.frame(patient_group = numeric(),
                               lower=numeric(),
                               median=numeric(),
                               upper=numeric(),
                               date=as.Date(character()))
  df_frail <- data.frame(patient_group = numeric(),
                                lower=numeric(),
                                median=numeric(),
                                upper=numeric(),
                                date=as.Date(character()))
  
  cat("This number of categories:",group_num,"\n")
  ad <- admit_type
  

  
  if(length(only_ICD) > 0)
    icds_to_look_at <- ts_data$ICD[grep(only_ICD,ts_data$ICD)[1]]
  else
    icds_to_look_at <- unique(ts_data$ICD)
  
  for ( age in unique(ts_data$agegrp_v3)){
      for (d in icds_to_look_at){
          cat("\r","On group number:", current_group)
          #print(paste("On group number:", current_group))
          current_group <- current_group + 1
          temp<-ts_data[ts_data$agegrp_v3==age &
                       ts_data$ICD==d,]
          temp <- temp[temp$date > as.Date("2009-04-30") &
                         temp$date < as.Date("2013-03-01") ,]
          
          if(cohort3 == TRUE){
            
          }
          
          
          test<-sum(temp$Admissions>0)>10
          if (test){
            if(forecast_admissions == TRUE){
            #ggplot(temp)+ aes(x=date,y=Admissions,colour=as.factor(agegrp_v2))+geom_line()
            
            # large number of admissions
            # 
            if(diagnostics_only == TRUE)
              Train<-temp[,c("date","Admissions")]
            else
              Train<-temp[temp$date<as.Date(train_date),c("date","Admissions")]
            Train<-Train[order(Train$date),]
            a1 <- matrix(c(0,0),2, 1)
            P1 <- matrix(0, 2, 2)
            P1inf <- diag(2)
            

            m<-as.formula('Admissions~SSMseasonal(period = 52.18, sea.type = "trigonometric") +
                                       SSMtrend(degree = 2, Q = list(NA,0),a1 = a1, P1inf = P1inf)')
            
            model <- SSModel(m,H = matrix(NA),
                             data=Train)
            
            
            fit <- fitSSM(model, inits = c(1.,1.),method = "BFGS")
            
            
            out <- KFS(fit$model,filtering='state',smoothing=c('state','disturbance','mean'))
            Train$predicted <-as.numeric(signal(out, states = c('season','trend'), filtered=FALSE )$signal)
            
            ###################################################################
            ## Diagnostics for admissions #####################################
            ###################################################################
            
            diagnostics_function(Train, metric = "Admissions", out,
                                 title_name = paste(d,age,ad) )
            
            if(diagnostics_only == FALSE){
              
              df<-temp[,c('date',"Admissions")]
              df<-df[order(df$date),]
              df$signal<-NA
              df[df$date <= max(Train$date),'signal']<-Train$predicted
              
              
              
              n<-nrow(temp[temp$date>max(Train$date),])
              
              v<-var(residuals(out,type="response"))
              newdata<-SSModel(rep(NA,n)~SSMseasonal(period = 52.18, sea.type = "trigonometric")+
                                 SSMtrend(degree = 2,Q =as.list(fit$model$Q[1:2]) ),H = fit$model$H)
              
              pred <- predict(fit$model, newdata=newdata,  interval = c( "prediction"), level = 0.95, 
                              states =c('all'), se.fit = FALSE, nsim = 1000, 
                              maxiter = 52, filtered = TRUE)
              pred<-as.data.frame(pred)
              df$fit<-NA
              df$upr<-NA
              df$lwr<-NA
              df[df$date>max(Train$date),c('fit','lwr','upr')]<-pred
              
                pred_dates <- df$date[df$date>max(Train$date)]
                patient_group <- rep(paste(ad,d,age, sep = "_"), length(pred_dates))
                rows_to_add <- cbind.data.frame(patient_group, pred, pred_dates)
                colnames(rows_to_add) <- c("patient_group","median","lower","upper","date")
                df_preds <- rbind.data.frame(df_preds, rows_to_add)
                
              
              
              print(ggplot(data=df)+aes(x=date,y=fit, colour = "Predicted_admission")+geom_line()+
                      geom_line(aes(x=date,y=signal,colour='Fit_to_data'))+
                      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.35, colour = NA
                      ) +
                      geom_line(aes(x=date,y=Admissions,colour='Admissions_actual'))+
                      theme_bw(base_size = 12, base_family = "Helvetica")+
                      theme(legend.position="right")+ylab("")+xlab("")+ labs(colour = "Legend") +
                      scale_color_manual(values = c("Fit_to_data" = "blue",
                                                    "Admissions_actual" = "red",
                                                    "Predicted_admission" = "black"))+
                      ggtitle(paste("Admissions",ad,d, age )))
              
            }
            }
            
            if(forecast_wait == TRUE & ad == "Elective"){
            
              
              
              if(diagnostics_only == TRUE)
                Train<-temp[,c("date",wait_time_col)]
              else
                Train<-temp[temp$date<as.Date(train_date),c("date",wait_time_col)]
              
            Train<-Train[order(Train$date),]
            a1 <- matrix(c(0,0),2, 1)
            P1 <- matrix(0, 2, 2)
            P1inf <- diag(2)
            
            formula_strings <- paste(wait_time_col, '~SSMseasonal(period = 52.18, sea.type = "trigonometric") +
                                       SSMtrend(degree = 2, Q = list(NA,0),a1 = a1, P1inf = P1inf)', sep = "")
            
            m<-as.formula(formula_strings)
            
            model <- SSModel(m,H = matrix(NA),
                             data=Train)
            print("Waiting time forecast")
            
            fit <- fitSSM(model, inits = c(1.,1.),method = "BFGS")
            
            
            out <- KFS(fit$model,filtering='state',smoothing=c('state','disturbance','mean'))
            Train$predicted <-as.numeric(signal(out, states = c('season','trend'), filtered=FALSE )$signal)
            
            ###################################################################
            ## Diagnostics ####################################################
            ###################################################################
            
            diagnostics_function(Train, metric = "Wait time", out,
                                 title_name = paste(d,age,ad) )
            
            if(diagnostics_only == FALSE){
              df<-temp[,c('date',wait_time_col)]
              df<-df[order(df$date),]
              df$signal<-NA
              df[df$date <= max(Train$date),'signal']<-Train$predicted
              
              
              
              n<-nrow(temp[temp$date>max(Train$date),])
              
              v<-var(residuals(out,type="response"))
              newdata<-SSModel(rep(NA,n)~SSMseasonal(period = 52.18, sea.type = "trigonometric")+
                                 SSMtrend(degree = 2,Q =as.list(fit$model$Q[1:2]) ),H = fit$model$H)
              
              pred <- predict(fit$model, newdata=newdata,  interval = c( "prediction"), level = 0.95, 
                              states =c('all'), se.fit = FALSE, nsim = 1000, 
                              maxiter = 52, filtered = TRUE)
              pred<-as.data.frame(pred)
              df$fit<-NA
              df$upr<-NA
              df$lwr<-NA
              df[df$date>max(Train$date),c('fit','lwr','upr')]<-pred
              
              
                
              pred_dates <- df$date[df$date>max(Train$date)]
              patient_group <- rep(paste(ad,d,age, sep = "_"), length(pred_dates))
              rows_to_add <- cbind.data.frame(patient_group, pred, pred_dates)
              colnames(rows_to_add) <- c("patient_group","median","lower","upper","date")
              df_wait <- rbind.data.frame(df_wait, rows_to_add)
                
              
              #b<-as.Date("2019-07-01")
              print(ggplot(data=df)+aes(x=date,y=fit, colour = "Predicted wait time")+geom_line()+
                      geom_line(aes(x=date,y=signal,colour='Fit_to_data'))+
                      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.45,
                                  linetype = 2, colour = NA) +
                      geom_line(aes(x=date,y=df[,2],colour='Actual_WT_data'))+
                      theme_bw(base_size = 12, base_family = "Helvetica")+
                      theme(legend.position="right")+ylab("")+xlab("")+
                      scale_color_manual(values = c("Fit_to_data" = "blue",
                                                    "Actual_WT_data" = "red",
                                                    "Predicted wait time" = "black"))+
                      ggtitle(paste("Waiting time",ad,age, d )))
            }
            }
            
            if(forecast_frail == TRUE){
              
            if(!any(grepl("prop_Frail",colnames(temp)))){
              temp$prop_Frail <- temp$Frail / temp$Admissions
            }
            
  
              
            if(diagnostics_only == TRUE)
              Train<-temp[,c("date","prop_Frail")]
            else
              Train<-temp[temp$date<as.Date(train_date),c("date","prop_Frail")]
            Train<-Train[order(Train$date),]
            a1 <- matrix(c(0,0),2, 1)
            P1 <- matrix(0, 2, 2)
            P1inf <- diag(2)
            
            m<-as.formula('prop_Frail~SSMseasonal(period = 52.18, sea.type = "trigonometric") +
                                         SSMtrend(degree = 2, Q = list(NA,0),a1 = a1, P1inf = P1inf)')
            
            model <- SSModel(m,H = matrix(NA),
                             data=Train)
            print("Proportion frail")
            
            fit <- fitSSM(model, inits = c(1.,1.),method = "BFGS")
            
            
            out <- KFS(fit$model,filtering='state',smoothing=c('state','disturbance','mean'))
            Train$predicted <-as.numeric(signal(out, states = c('season','trend'), filtered=FALSE )$signal)
            
            
            diagnostics_function(Train, metric = "Proportion frail", out,
                                 title_name = paste(d,age,ad) )
            
            
            if(diagnostics_only == FALSE){
            
              df<-temp[,c('date',"prop_Frail")]
              df<-df[order(df$date),]
              df$signal<-NA
              df[df$date <= max(Train$date),'signal']<-Train$predicted
              
              
              
              n<-nrow(temp[temp$date>max(Train$date),])
              
              v<-var(residuals(out,type="response"))
              newdata<-SSModel(rep(NA,n)~SSMseasonal(period = 52.18, sea.type = "trigonometric")+
                                 SSMtrend(degree = 2,Q =as.list(fit$model$Q[1:2]) ),H = fit$model$H)
              
              pred <- predict(fit$model, newdata=newdata,  interval = c( "prediction"), level = 0.95, 
                              states =c('all'), se.fit = FALSE, nsim = 1000, 
                              maxiter = 52, filtered = TRUE)
              pred<-as.data.frame(pred)
              df$fit<-NA
              df$upr<-NA
              df$lwr<-NA
              df[df$date>max(Train$date),c('fit','lwr','upr')]<-pred
              
                pred_dates <- df$date[df$date>max(Train$date)]
                patient_group <- rep(paste(ad,d,age, sep = "_"), length(pred_dates))
                rows_to_add <- cbind.data.frame(patient_group, pred, pred_dates)
                colnames(rows_to_add) <- c("patient_group","median","lower","upper","date")
                df_frail <- rbind.data.frame(df_frail, rows_to_add)
                
              
              
              #b<-as.Date("2019-07-01")
              print(ggplot(data=df)+aes(x=date,y=fit, colour = "Predicted_prop_Frail")+geom_line()+
                      geom_line(aes(x=date,y=signal,colour='Fit_to_data'))+
                      geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.45,
                                  colour = NA) +
                      geom_line(aes(x=date,y=prop_Frail,colour='prop_Frail_data'))+
                      theme_bw(base_size = 12, base_family = "Helvetica")+
                      theme(legend.position="right")+ylab("")+xlab("")+
                      scale_color_manual(values = c("Fit_to_data" = "blue",
                                                    "prop_Frail_data" = "red",
                                                    "Predicted_prop_Frail" = "black"))+
                      ggtitle(paste("Frailty",ad,age, d )))
            }
            }
          }
        }}
  
  dev.off()
  time_end <- Sys.time()
  print((time_end - time_start))
  return(list(df_preds, df_wait, df_frail))
  
}


test_cohort_1 <- forecast_function(cohort_1_ts, results_name = "cohort_1_admissions_mar_train_frail.pdf",
                                   forecast_admissions = TRUE, diagnostics_only = FALSE, forecast_frail = TRUE,
                                   admit_type = "Emergency",
                                   train_date = "2012-02-25")

test_cohort_2 <- forecast_function(cohort_2_ts, results_name = "cohort_2_admissions.pdf",
                                   forecast_admissions = TRUE, diagnostics_only = FALSE,
                                   admit_type = "Emergency",
                                   train_date = "2012-02-01")

test_cohort_3_new_median <- forecast_function(cohort_3_ts, results_name = "cohort_3_new_patients.pdf",
                                       forecast_admissions = TRUE, forecast_frail = FALSE,
                                       forecast_wait = TRUE,
                                       admit_type = "Elective", cohort3 = TRUE,
                                       wait_time_col = "p50_WT_ICDc", train_date = "2012-02-25")

test_cohort_3_new_mean <- forecast_function(cohort_3_ts, results_name = "cohort_3_new_patientsmean.pdf",
                                            forecast_admissions = FALSE, forecast_wait = TRUE, forecast_frail = TRUE,
                                       admit_type = "Elective", cohort3 = TRUE,
                                       wait_time_col = "mean_WT_ICDc", diagnostics_only = FALSE,
                                       train_date = "2012-02-25", only_ICD = "Neoplasms")


###############################################################################
## Lets get the neoplasms out of there! #######################################
###############################################################################

cohort_1_admi <- test_cohort_1[[1]]
cohort_1_frail <- test_cohort_1[[3]]

cohort_3_admi <- test_cohort_3_new_median[[1]]
cohort_3_median_wait <- test_cohort_3_new_median[[2]]

cohort_3_mean_wait <- test_cohort_3_new_mean[[2]]
cohort_3_frail <- test_cohort_3_new_mean[[3]]



cohort_1_neos <- grep("Neoplasms", cohort_1_admi$patient_group, ignore.case = TRUE)
cohort_1_admi_neos <- cohort_1_admi[cohort_1_neos,]

patient_groupings <- str_split_fixed(cohort_1_admi_neos$patient_group, "_", 3)
cohort_1_admi_neos$cohort <- 1
cohort_1_admi_neos$age <- patient_groupings[,3]
cohort_1_admi_neos$icd_name <- patient_groupings[,2]
cohort_1_admi_neos$icd_num <- 2

cohort_1_neos_frail <- grep("Neoplasms", cohort_1_frail$patient_group, ignore.case = TRUE)
cohort_1_frail_neos <- cohort_1_frail[cohort_1_neos_frail,]
patient_groupings <- str_split_fixed(cohort_1_frail_neos$patient_group, "_", 3)
cohort_1_frail_neos$cohort <- 1
cohort_1_frail_neos$age <- patient_groupings[,3]
cohort_1_frail_neos$icd_name <- patient_groupings[,2]
cohort_1_frail_neos$icd_num <- 2



write.csv(cohort_1_admi_neos,
          file = "~/Dropbox/covid19/overflow_analysis/time_series_data/neoplasms_admissions_cohort_1.csv",
          row.names = FALSE)
write.csv(cohort_1_frail_neos,
          file = "~/Dropbox/covid19/overflow_analysis/time_series_data/neoplasms_frail_cohort_1.csv",
          row.names = FALSE)


cohort_3_neos <- grep("Neoplasms", cohort_3_admi$patient_group, ignore.case = TRUE)
cohort_3_admi_neos <- cohort_3_admi[cohort_3_neos,]
patient_groupings <- str_split_fixed(cohort_3_admi_neos$patient_group, "_", 3)
cohort_3_admi_neos$cohort <- 3
cohort_3_admi_neos$age <- patient_groupings[,3]
cohort_3_admi_neos$icd_name <- patient_groupings[,2]
cohort_3_admi_neos$icd_num <- 2



cohort_3_neos_wait <- grep("Neoplasms", cohort_3_median_wait$patient_group, ignore.case = TRUE)
cohort_3_wait_neos <- cohort_3_median_wait[cohort_3_neos_wait,]
patient_groupings <- str_split_fixed(cohort_3_wait_neos$patient_group, "_", 3)
cohort_3_wait_neos$cohort <- 3
cohort_3_wait_neos$age <- patient_groupings[,3]
cohort_3_wait_neos$icd_name <- patient_groupings[,2]
cohort_3_wait_neos$icd_num <- 2

write.csv(cohort_3_admi_neos,
          file = "~/Dropbox/covid19/overflow_analysis/time_series_data/neoplasms_new_pool_cohort_3.csv",
          row.names = FALSE)
write.csv(cohort_3_wait_neos,
          file = "~/Dropbox/covid19/overflow_analysis/time_series_data/neoplasms_median_wait_cohort_3.csv",
          row.names = FALSE)


cohort_3_meaners <- grep("Neoplasms", cohort_3_mean_wait$patient_group, ignore.case = TRUE)
cohort_3_wait_mean <- cohort_3_mean_wait[cohort_3_meaners,]
patient_groupings <- str_split_fixed(cohort_3_wait_mean$patient_group, "_", 3)
cohort_3_wait_mean$cohort <- 3
cohort_3_wait_mean$age <- patient_groupings[,3]
cohort_3_wait_mean$icd_name <- patient_groupings[,2]
cohort_3_wait_mean$icd_num <- 2



cohort_3_frail_grep <- grep("Neoplasms", cohort_3_frail$patient_group, ignore.case = TRUE)
cohort_3_frail_new <- cohort_3_frail[cohort_3_frail_grep,]
patient_groupings <- str_split_fixed(cohort_3_frail_new$patient_group, "_", 3)
cohort_3_frail_new$cohort <- 3
cohort_3_frail_new$age <- patient_groupings[,3]
cohort_3_frail_new$icd_name <- patient_groupings[,2]
cohort_3_frail_new$icd_num <- 2

write.csv(cohort_3_wait_mean,
          file = "~/Dropbox/covid19/overflow_analysis/time_series_data/neoplasms_mean_wait_cohort_3.csv",
          row.names = FALSE)
write.csv(cohort_3_frail_new,
          file = "~/Dropbox/covid19/overflow_analysis/time_series_data/neoplasms_frail_cohort_3.csv",
          row.names = FALSE)






nrow(cohort_1_admi_neos)/3



cohort_1_admi_neos$frail <- cohort_1_admi_neos$median * cohort_1_admi_neos$median



neoplasms_wt_fr <- forecast_function(data, results_name = "admissions_diagnostics.pdf",forecast_wait = FALSE,
                  forecast_admissions = TRUE, forecast_frail = FALSE)

write.csv(neoplasms_wt_fr[[2]],
          file = "~/Dropbox/covid19/overflow_analysis/neoplasms_WT_2012-13.csv",
          row.names = FALSE)
write.csv(neoplasms_wt_fr[[3]],
          file = "~/Dropbox/covid19/overflow_analysis/neoplasms_Frailty_2012-13.csv",
          row.names = FALSE)


write.csv(neoplasms,
          file = "~/Dropbox/covid19/overflow_analysis/neoplasm_predictions_2012-2013.csv",
          row.names = FALSE)


###############################################################################
## Forecasting for the different cats #########################################
###############################################################################

indiviudal_plot_function <- function(temp, ad, d, age){
  
  test<-sum(temp$Admissions>0)>10
  if (test){
    #ggplot(temp)+ aes(x=date,y=Admissions,colour=as.factor(agegrp_v2))+geom_line()
    
    # large number of admissions
    # 
    Train<-temp[temp$admidate_YYYY<2012,c("date","Admissions")]
    Train<-Train[order(Train$date),]
    a1 <- matrix(c(0,0),2, 1)
    P1 <- matrix(0, 2, 2)
    P1inf <- diag(2)
    
    m<-as.formula('Admissions~SSMseasonal(period = 52.18, sea.type = "trigonometric") +
                                     SSMtrend(degree = 2, Q = list(NA,0),a1 = a1, P1inf = P1inf)')
    
    model <- SSModel(m,H = matrix(NA),
                     data=Train)
    print(paste(ad,d, age ))
    print("Estimation model with local linear trend and trigonometric seasonality")
    fit <- fitSSM(model, inits = c(1.,1.),method = "BFGS")
    
    print("Smoothing")
    out <- KFS(fit$model,filtering='state',smoothing=c('state','disturbance','mean'))
    Train$predicted <-as.numeric(signal(out, states = c('season','trend'), filtered=FALSE )$signal)
    
    df<-temp[,c('date',"Admissions")]
    df<-df[order(df$date),]
    df$signal<-NA
    df[df$date <= max(Train$date),'signal']<-Train$predicted
    
    
    
    n<-nrow(temp[temp$date>max(Train$date),])
    
    v<-var(residuals(out,type="response"))
    newdata<-SSModel(rep(NA,n)~SSMseasonal(period = 52.18, sea.type = "trigonometric")+
                       SSMtrend(degree = 2,Q =as.list(fit$model$Q[1:2]) ),H = fit$model$H)
    
    pred <- predict(fit$model, newdata=newdata,  interval = c( "prediction"), level = 0.95, 
                    states =c('all'), se.fit = FALSE, nsim = 1000, 
                    maxiter = 52, filtered = TRUE)
    pred<-as.data.frame(pred)
    df$fit<-NA
    df$upr<-NA
    df$lwr<-NA
    df[df$date>max(Train$date),c('fit','lwr','upr')]<-pred
    
    #b<-as.Date("2019-07-01")
    print(ggplot(data=df)+aes(x=date,y=fit, colour = "Predicted_median")+geom_line()+
            geom_line(aes(x=date,y=signal,colour='Fit_to_data'))+
            geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.2) +
            geom_line(aes(x=date,y=Admissions,colour='Admissions_data'))+
            theme_bw(base_size = 12, base_family = "Helvetica")+
            theme(legend.position="right")+ylab("")+xlab("")+
            scale_color_manual(values = c("Fit_to_data" = "blue",
                                          "Admissions_data" = "red",
                                          "Predicted_median" = "black"))+
            ggtitle(paste(ad,d, age )))
  }
  
}



forecast_function_altered_grouping <- function(ts_data, results_name = "results_agg.pdf"){

  ## Run this for all data, data by admit type, data by age, admit x age, admit x ICD
  
  # join the year and week variable to a date. Use UK  convention that week starts on a Monday.
  ts_data$date<-as.POSIXct( paste( 1, ts_data$admidate_week, ts_data$admidate_YYYY, sep = "-" ), format = "%u-%U-%Y",locale = "UK" ) 
  
#  pdf(file = results_name ,paper='A4')
  
  ## All ##

  dates <- unique(ts_data$date)
  new_ts_df <- data.frame(matrix(ncol = 3, nrow = length(dates)))
  colnames(new_ts_df) <- c("date","Admissions","admidate_YYYY")
  new_ts_df$date <- dates
  
  for(k in 1:nrow(new_ts_df)){
    current_date_data <- ts_data[ts_data$date == new_ts_df$date[k],]
    current_admissions <- sum(current_date_data$Admissions)
    new_ts_df$Admissions[k] <- current_admissions
    new_ts_df$admidate_YYYY[k] <- current_date_data$admidate_YYYY[1]
    
  }
  
  new_ts_df <- new_ts_df[new_ts_df$date > as.Date("2009-04-30"),]
  indiviudal_plot_function(new_ts_df, "Emergency & Elective", "All ICDs","All Ages")
  
  ## Admit ##
  
  for(k in c("Emergency","Elective")){
    new_ts_df <- data.frame(matrix(ncol = 3, nrow = length(dates)))
    colnames(new_ts_df) <- c("date","Admissions","admidate_YYYY")
    new_ts_df$date <- dates
    
    for(j in 1:nrow(new_ts_df)){
      current_date_data <- ts_data[ts_data$date == new_ts_df$date[j] &
                                   ts_data$admimeth_C== k,]
      current_admissions <- sum(current_date_data$Admissions)
      new_ts_df$Admissions[j] <- current_admissions
      new_ts_df$admidate_YYYY[j] <- current_date_data$admidate_YYYY[1]
      
    }
    new_ts_df <- new_ts_df[new_ts_df$date > as.Date("2009-04-30"),]
    indiviudal_plot_function(new_ts_df, k, "All ICDs","All ages")
  }
  
  ## Ages ##
  
  for(age in unique(ts_data$agegrp_v3)){
    new_ts_df <- data.frame(matrix(ncol = 3, nrow = length(dates)))
    colnames(new_ts_df) <- c("date","Admissions","admidate_YYYY")
    new_ts_df$date <- dates
    
    for(j in 1:nrow(new_ts_df)){
      current_date_data <- ts_data[ts_data$date == new_ts_df$date[j] &
                                     ts_data$agegrp_v3== age,]
      current_admissions <- sum(current_date_data$Admissions)
      new_ts_df$Admissions[j] <- current_admissions
      new_ts_df$admidate_YYYY[j] <- current_date_data$admidate_YYYY[1]
      
    }
    new_ts_df <- new_ts_df[new_ts_df$date > as.Date("2009-04-30"),]
    indiviudal_plot_function(new_ts_df, "Emergency & Elective", "All ICDs",age)
  }
  
  ## Admissions by age ##
  
  for(k in c("Emergency","Elective")){
    for(age in unique(ts_data$agegrp_v3)){
      new_ts_df <- data.frame(matrix(ncol = 3, nrow = length(dates)))
      colnames(new_ts_df) <- c("date","Admissions","admidate_YYYY")
      new_ts_df$date <- dates
      
      for(j in 1:nrow(new_ts_df)){
        current_date_data <- ts_data[ts_data$date == new_ts_df$date[j] &
                                       ts_data$agegrp_v3== age &
                                       ts_data$admimeth_C== k,]
        current_admissions <- sum(current_date_data$Admissions)
        new_ts_df$Admissions[j] <- current_admissions
        new_ts_df$admidate_YYYY[j] <- current_date_data$admidate_YYYY[1]
        
      }
      new_ts_df <- new_ts_df[new_ts_df$date > as.Date("2009-04-30"),]
      indiviudal_plot_function(new_ts_df, k, "All ICDs",age)
    }
  }
  
  ## Admissions by ICD ##
  
  for (d in unique(ts_data$ICD)){
    for(k in c("Emergency","Elective")){
      new_ts_df <- data.frame(matrix(ncol = 3, nrow = length(dates)))
      colnames(new_ts_df) <- c("date","Admissions","admidate_YYYY")
      new_ts_df$date <- dates
      
      for(j in 1:nrow(new_ts_df)){
        current_date_data <- ts_data[ts_data$date == new_ts_df$date[j] &
                                       ts_data$ICD== d &
                                       ts_data$admimeth_C== k,]
        current_admissions <- sum(current_date_data$Admissions)
        new_ts_df$Admissions[j] <- current_admissions
        new_ts_df$admidate_YYYY[j] <- current_date_data$admidate_YYYY[1]
        
      }
      new_ts_df <- new_ts_df[new_ts_df$date > as.Date("2009-04-30"),]
      indiviudal_plot_function(new_ts_df, k, d,"All ages")
    }}
  
  
  
  
  dev.off()
  
}

forecast_function_altered_grouping(ts_data = data)
plot(seq(0,100,1))

