# subset only icd 15 from transitions_data

icd_15 <- subset(transitions_data, subset = ICD == 15)

# subset cohort 3
icd_15_cohort3 <- subset(icd_15, subset = cohort == 3)
CI.byagegrp <- cuminc(ftime = icd_15_cohort3$cc_LoS, fstatus = icd_15_cohort3$cc_transitions, group = icd_15_cohort3$agegrp_v3)

cif_11_time <- CI.byagegrp[["1 1"]][["time"]]
cif_11_est <- CI.byagegrp[["1 1"]][["est"]]

df_cif_11 <- cbind.data.frame(cif_11_time, cif_11_est)

df_cif_11 <- df_cif_11[order(df_cif_11$cif_11_time, -(df_cif_11$cif_11_est)),]

df_cif_11 <- df_cif_11[!duplicated(df_cif_11$cif_11_time),]

colnames(df_cif_11) <- c("time","CIF")

max_time <- max(df_cif_11$time)
multiplesofseven <- seq(7,max_time, by = 7)
seventh_day_timies <- timepoints(CI.byagegrp,multiplesofseven) 

test <- cbind.data.frame(seventh_day_timies[["est"]])

test$mean_7 <- rowMeans(test, na.rm = TRUE)
test$median_7 <- apply(test[-11], 1, FUN = median, na.rm = TRUE)

colnames(test)[1] <- "day7"

icd15export <- cbind.data.frame(test$day7, test$mean_7, test$median_7)

colnames(icd15export) <- c("day7","mean_7","median_7")

write.csv(icd15export,
          file = "D:/Dropbox/COVID19/Overflow/icd15export.csv",
          row.names = FALSE,
          quote = FALSE)
