##############################################################################################################
## Costing analysis ##########################################################################################
##############################################################################################################
# Lines where Fiona must change the file location: 26, 33-46, 222, 308-326
##############################################################################################################
# Required packages

costs_function <- function(transitions_data_whole, costs_directory, output_dir, fyear){
  
  
  ##############################################################################################################
  # Notes about merging HES and costs data:
  # ORG.CODE matched to procode3 
  # Currency.code matched to HRG or susHRG (this is for post FY1213 onwards. will need to include a change in variable name for HF pre cleaning file)
  # we need to keep the UNIT COST, MEAN COST, Actual Cost MAPPING POT and PERIOD (2009-10) EpiEnd_FY== "0910"
  # mean cost = the national mean average unit cost
  # unit cost = the average cost to the organisation of providing that activity
  # actual cost = the organisation's activity (number of cases) x the organisation unit cost
  # we are really only interested in matching the unit cost of each organisation's hrg
  # Mapping pot is matched to admimeth_C
  # Mapping to 01_EI and 02_NEI and 03_XS
  # XS bedays with the bed.days variable
  
  ##############################################################################################################
  ## Load HES data
  
  if(class(transitions_data_whole) == "character"){
    system.time(transitions_data_whole <- readRDS(here::here("time_series_forecasts", "time_series_forecasts.Rds")))  
  }
  
  # keep variables you need
  hes_data <- select(transitions_data_whole, c("procode3","admimeth_C","EpiEnd_FY","SUSHRGep", "ICD", "agegrp_v3", "admimeth_C", "WaitingTime"))
  
  # Pre-process HES data
  # Keep HES data of yoi (FY of interest)
  hes_data_yoi <- subset(hes_data, EpiEnd_FY == as.numeric(fyear))
  
  ##############################################################################################################
  ## Load costs data (all already inflated to 2020 costs)
  # Load organisational level unit costs
  if(substr(costs_directory,nchar(costs_directory),nchar(costs_directory)) != "/"){
    costs_directory <- paste(costs_directory,"/",sep = "")
  }
  
  
  costs_data_path <- paste(costs_directory,"cleaned_2015_2019_orgrcs.csv",sep = "")
  system.time(costs_dataset <- vroom(costs_data_path, delim = ","))
  
  # Load previous years reference costs to match HES admissions that are left unmatched
  #load national ref cost data for 1819 - for costs that are unmatched using organisational level cost data
  rcs1819_path <- paste(costs_directory,"national_rcs1819.csv",sep = "")
  national_rcs1819 <- read.csv(rcs1819_path, fileEncoding = "UTF-8-BOM")
  
  # load national ref cost data for 1718 - for costs that are unmatched after using national_rcs1213
  rcs1718_path <- paste(costs_directory,"national_rcs1718.csv",sep = "")
  national_rcs1718 <- read.csv(rcs1718_path,fileEncoding = "UTF-8-BOM")
  
  # load national ref cost data for 1617 - for costs that are unmatched after using national_rcs1112
  rcs1617_path <- paste(costs_directory,"national_rcs1617.csv",sep = "")
  national_rcs1617 <- read.csv(rcs1617_path,fileEncoding = "UTF-8-BOM")
  
  # load national ref cost data for 1516 - for costs that are unmatched after using national_rcs1011
  rcs1516_path <- paste(costs_directory,"national_rcs1516.csv",sep = "")
  national_rcs1516 <- read.csv(rcs1516_path,fileEncoding = "UTF-8-BOM")
  
  ##############################################################################################################
  # Pre-process Costs data
  # Split organizational reference costs into each individual year
  org_costs_1516 <- subset(costs_dataset, EpiEnd_FY == 1516)
  org_costs_1617 <- subset(costs_dataset, EpiEnd_FY == 1617)
  org_costs_1718 <- subset(costs_dataset, EpiEnd_FY == 1718)
  org_costs_1819 <- subset(costs_dataset, EpiEnd_FY == 1819)
  
  ##############################################################################################################
  # Now merge HES and reference costs
  # Since not all HES episodes will match with a cost on the first try 
  # (procode3-SUSHRGep pair does not exist in costs data but does exist in hes data), 
  # we will match in the following order:
  # 1. On organization level reference costs (1819)
  # 2. On national reference cost schedule (1819)
  # 3. On organization level reference costs (1718)
  # 4. On national reference cost schedule (1718)
  # 5. On organization level reference costs (1617)
  # 6. On national reference cost schedule (1617)
  # 7. On organization level reference costs (1516)
  # 8. On national reference cost schedule (1516)
  ############################################################################
  # 1. On organization level reference costs (1819)
  # match hes_data_yoi with org_costs_1819
  merged_orglevel_1819 <- left_join(hes_data_yoi, org_costs_1819, 
                                    by = c("procode3" = "orgcode",
                                           "admimeth_C" = "admimeth_c",
                                           "SUSHRGep" = "currencycode", "EpiEnd_FY"))
  # subset MATCHED and UNMATCHED episodes
  matched_orglevel_1819 <- merged_orglevel_1819[!is.na(merged_orglevel_1819$unitcost),]
  unmatched_orglevel_1819 <- merged_orglevel_1819[is.na(merged_orglevel_1819$unitcost),]
  
  # remove unitcost column from unmatched_orglevel_1819
  unmatched_orglevel_1819 <- select(unmatched_orglevel_1819, -unitcost)
  ############################################################################
  # 2. On national reference cost schedule (1819)
  # match unmatched_orglevel_1819 with national_rcs1819
  merged_rcslevel_1819 <- left_join(unmatched_orglevel_1819, national_rcs1819,
                                    by = c("admimeth_C" = "admimeth_c",
                                           "SUSHRGep" = "currencycode"))
  # subset MATCHED and UNMATCHED episodes
  matched_rcslevel_1819 <- merged_rcslevel_1819[!is.na(merged_rcslevel_1819$unitcost),]
  unmatched_rcslevel_1819 <- merged_rcslevel_1819[is.na(merged_rcslevel_1819$unitcost),]
  
  # remove unitcost column from unmatched_rcslevel_1819
  unmatched_rcslevel_1819 <- select(unmatched_rcslevel_1819, -unitcost, -meancost, -actualcost, -expectedcost)
  
   ############################################################################
  # 3. On organization level reference costs (1718)
  # match unmatched_rcslevel_1819 with org_costs_1718
  merged_orglevel_1718 <- left_join(unmatched_rcslevel_1819, org_costs_1718,
                                    by = c("procode3" = "orgcode",
                                           "admimeth_C" = "admimeth_c",
                                           "SUSHRGep" = "currencycode", "EpiEnd_FY"))
  # subset MATCHED and UNMATCHED episodes
  matched_orglevel_1718 <- merged_orglevel_1718[!is.na(merged_orglevel_1718$unitcost),]
  unmatched_orglevel_1718 <- merged_orglevel_1718[is.na(merged_orglevel_1718$unitcost),]
  
  # remove unitcost column from unmatched_orglevel_1718
  unmatched_orglevel_1718 <- select(unmatched_orglevel_1718, -unitcost)
  ############################################################################
  # 4. On national reference cost schedule (1718)
  # match unmatched_orglevel_1718 with national_rcs1718
  merged_rcslevel_1718 <- left_join(unmatched_orglevel_1718, national_rcs1718,
                                    by = c("admimeth_C" = "admimeth_c",
                                           "SUSHRGep" = "currencycode"))
  # subset MATCHED and UNMATCHED episodes
  matched_rcslevel_1718 <- merged_rcslevel_1718[!is.na(merged_rcslevel_1718$unitcost),]
  unmatched_rcslevel_1718 <- merged_rcslevel_1718[is.na(merged_rcslevel_1718$unitcost),]
  
  # remove unitcost column from unmatched_rcslevel_1718
  unmatched_rcslevel_1718 <- select(unmatched_rcslevel_1718, -unitcost, -meancost, -actualcost, -expectedcost)
  ############################################################################
  # 5. On organization level reference costs (1617)
  # match unmatched_rcslevel_1718 with org_costs_1617
  merged_orglevel_1617 <- left_join(unmatched_rcslevel_1718, org_costs_1617,
                                    by = c("procode3" = "orgcode",
                                           "admimeth_C" = "admimeth_c",
                                           "SUSHRGep" = "currencycode", "EpiEnd_FY"))
  # subset MATCHED and UNMATCHED episodes
  matched_orglevel_1617 <- merged_orglevel_1617[!is.na(merged_orglevel_1617$unitcost),]
  unmatched_orglevel_1617 <- merged_orglevel_1617[is.na(merged_orglevel_1617$unitcost),]
  
  # remove unitcost column from unmatched_orglevel_1617
  unmatched_orglevel_1617 <- select(unmatched_orglevel_1617, -unitcost)
  ############################################################################
  # 6. On national reference cost schedule (1617)
  # match unmatched_orglevel_1617 with national_rcs1617
  merged_rcslevel_1617 <- left_join(unmatched_orglevel_1617, national_rcs1617,
                                    by = c("admimeth_C" = "admimeth_c",
                                           "SUSHRGep" = "currencycode"))
  # subset MATCHED and UNMATCHED episodes
  matched_rcslevel_1617 <- merged_rcslevel_1617[!is.na(merged_rcslevel_1617$unitcost),]
  unmatched_rcslevel_1617 <- merged_rcslevel_1617[is.na(merged_rcslevel_1617$unitcost),]
  
  # remove unitcost column from unmatched_rcslevel_1718
  unmatched_rcslevel_1617 <- select(unmatched_rcslevel_1617, -unitcost, -meancost, -actualcost, -expectedcost)
  ############################################################################
  # 7. On organization level reference costs (1516)
  # match unmatched_rcslevel_1617 with org_costs_1516
  merged_orglevel_1516 <- left_join(unmatched_rcslevel_1617, org_costs_1516,
                                    by = c("procode3" = "orgcode",
                                           "admimeth_C" = "admimeth_c",
                                           "SUSHRGep" = "currencycode", "EpiEnd_FY"))
  # subset MATCHED and UNMATCHED episodes
  matched_orglevel_1516 <- merged_orglevel_1516[!is.na(merged_orglevel_1516$unitcost),]
  unmatched_orglevel_1516 <- merged_orglevel_1516[is.na(merged_orglevel_1516$unitcost),]
  
  # remove unitcost column from unmatched_orglevel_1516
  unmatched_orglevel_1516 <- select(unmatched_orglevel_1516, -unitcost)
  ############################################################################
  # 8. On national reference cost schedule (1516)
  # match unmatched_orglevel_1516 with national_rcs1516
  merged_rcslevel_1516 <- left_join(unmatched_orglevel_1516, national_rcs1516,
                                    by = c("admimeth_C" = "admimeth_c",
                                           "SUSHRGep" = "currencycode"))
  ############################################################################
  # append all matched dataframes to form complete hes-cost merge
  full_mergedcosts <- rbind(matched_orglevel_1819,
                            merged_rcslevel_1819,
                            matched_orglevel_1718,
                            merged_rcslevel_1718,
                            matched_orglevel_1617,
                            merged_rcslevel_1617,
                            matched_orglevel_1516,
                            merged_rcslevel_1516)
  
  # figure out how many unmatched episodes we still have after this
  print('Number of uncosted HES episodes')
  unmatched <- sum(is.na(full_mergedcosts$unitcost))
  unmatched
  
  print('Number of HES episodes')
  total <- nrow(full_mergedcosts)
  total
  
  print('Percentage of HES episodes unmatched')
  (unmatched/total)*100
  
  # drop unmatched from full_mergedcosts$unitcost
  full_mergedcosts <- full_mergedcosts[!is.na(full_mergedcosts$unitcost),]
  saveRDS(full_mergedcosts, paste0(output_dir, "/mergedcosts_", fyear, ".Rds"))
  
}  
  
calculate_unitcosts <- function(full_mergedcosts, costs_directory, output_dir){
  ##############################################################################################################
  # Now we move onto calculating unit costs
  # first, calculate average unit cost per patient group (mean, sd, min, p25, p50, p75, max)
  quantilecost <- full_mergedcosts %>% group_by(ICD, agegrp_v3, admimeth_C) %>%
    do(data.frame(t(quantile(.$unitcost, props = c(0.25, 0.50, 0.75)))))
  
  avgcost <- do.call(data.frame, aggregate(full_mergedcosts$unitcost, 
                                           list(full_mergedcosts$ICD, full_mergedcosts$agegrp_v3, full_mergedcosts$admimeth_C),
                                           function(x) c(mean = mean(x), sd = sd(x), n = length(x))))
  
  colnames(avgcost) <- c("ICD","agegrp","admimeth","avg_cost","sd_cost", "n")
  
  costs <- left_join(avgcost, quantilecost,
                     by = c("ICD" = "ICD",
                            "agegrp" = "agegrp_v3",
                            "admimeth" = "admimeth_C"))
  costs <- rename(costs, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
  
  costs$disease <- "ICD"
  
  costs$age <- paste("_AGE", costs$agegrp, sep = "")
  
  costs$p <- paste0(costs$disease, costs$ICD,costs$age)
  costs$s <- "G"
  costs$a <- ifelse(costs$admimeth == 1, "N", "E")
  
  costs_export <- select(costs, -c(ICD, agegrp, admimeth, disease, age))
  
  # Write out costs_export to csv
  if(substr(output_dir,nchar(output_dir),nchar(output_dir)) != "/")
    output_dir <- paste(output_dir,"/",sep = "")
  
  costs_export_loc <- paste(output_dir,"costs.csv")
  
  write.csv(costs_export,
            file = costs_export_loc,
            row.names = FALSE,
            quote = FALSE)
  
  ###########################################################################################
  # Now calculate average unit cost per patient group across waiting time quartiles for electives (eg. for p25 WT, what is the cost p25, p50, p75, mean, sd)
  electives_merged <- subset(full_mergedcosts, admimeth_C == 1)
  
  quantilewt <- electives_merged %>% group_by(ICD, agegrp_v3) %>%
    do(data.frame(t(quantile(.$WaitingTime, props = c(0.25, 0.50, 0.75)))))
  
  quantilewt <- rename(quantilewt, min = X0., p25 = X25., p50 = X50., p75 = X75., max = X100.)
  
  # at each quartiles, what is the avg, sd, min, 25, 50, 75, max cost?
  # identify which electives fall within each quartiles
  electives_wtquant <- left_join(electives_merged, quantilewt,
                                 by = c("ICD","agegrp_v3"))
  
  electives_wtquant$p025 <- ifelse(electives_wtquant$WaitingTime < electives_wtquant$p25, 1, 0)
  electives_wtquant$p2550 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$p25 & 
                                       electives_wtquant$WaitingTime < electives_wtquant$p50), 1, 0)
  electives_wtquant$p5075 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$p50 & 
                                       electives_wtquant$WaitingTime < electives_wtquant$p75), 1, 0)
  electives_wtquant$p75100 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$p75), 1, 0)
  
  electives_p025 <- subset(electives_wtquant, p025 == 1)
  electives_p2550 <- subset(electives_wtquant, p2550 == 1)
  electives_p5075 <- subset(electives_wtquant, p5075 == 1)
  electives_p75100 <- subset(electives_wtquant, p75100 == 1)
  
  # quantiles of costs within each WT quartiles
  quant_costs_wtp025 <- electives_p025 %>% group_by(ICD, agegrp_v3, p025) %>%
    do(data.frame(t(quantile(.$unitcost, props = c(0.25, 0.50, 0.75)))))
  
  quant_costs_wtp2550 <- electives_p2550 %>% group_by(ICD, agegrp_v3, p2550) %>%
    do(data.frame(t(quantile(.$unitcost, props = c(0.25, 0.50, 0.75)))))
  
  quant_costs_wtp5075 <- electives_p5075 %>% group_by(ICD, agegrp_v3, p5075) %>%
    do(data.frame(t(quantile(.$unitcost, props = c(0.25, 0.50, 0.75)))))
  
  quant_costs_wtp75100 <- electives_p75100 %>% group_by(ICD, agegrp_v3, p75100) %>%
    do(data.frame(t(quantile(.$unitcost, props = c(0.25, 0.50, 0.75)))))
  
  # mean and sd cost within each WT quartiles
  avg_costs_wtp025 <- do.call(data.frame, aggregate(electives_p025$unitcost, 
                                                    list(electives_p025$ICD, electives_p025$agegrp_v3),
                                                    function(x) c(mean = mean(x), sd = sd(x), n = length(x))))
  colnames(avg_costs_wtp025) <- c("ICD","agegrp_v3","avg_cost","sd_cost", "n")
  
  avg_costs_wtp2550 <- do.call(data.frame, aggregate(electives_p2550$unitcost, 
                                                     list(electives_p2550$ICD, electives_p2550$agegrp_v3),
                                                     function(x) c(mean = mean(x), sd = sd(x), n = length(x))))
  colnames(avg_costs_wtp2550) <- c("ICD","agegrp_v3","avg_cost","sd_cost", "n")
  
  avg_costs_wtp5075 <- do.call(data.frame, aggregate(electives_p5075$unitcost, 
                                                     list(electives_p5075$ICD, electives_p5075$agegrp_v3),
                                                     function(x) c(mean = mean(x), sd = sd(x), n = length(x))))
  colnames(avg_costs_wtp5075) <- c("ICD","agegrp_v3","avg_cost","sd_cost", "n")
  
  avg_costs_wtp75100 <- do.call(data.frame, aggregate(electives_p75100$unitcost, 
                                                      list(electives_p75100$ICD, electives_p75100$agegrp_v3),
                                                      function(x) c(mean = mean(x), sd = sd(x), n = length(x))))
  colnames(avg_costs_wtp75100) <- c("ICD","agegrp_v3","avg_cost","sd_cost", "n")
  
  # combine mean, sd, and quantile costs at each WT quartiles
  costs_wtp025 <- left_join(avg_costs_wtp025, quant_costs_wtp025,
                            by = c("ICD" = "ICD","agegrp_v3" = "agegrp_v3"))
  costs_wtp025 <- rename(costs_wtp025, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
  costs_wtp025 <- subset(costs_wtp025, select = -p025)
  
  costs_wtp2550 <- left_join(avg_costs_wtp2550, quant_costs_wtp2550,
                             by = c("ICD" = "ICD","agegrp_v3" = "agegrp_v3"))
  costs_wtp2550 <- rename(costs_wtp2550, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
  costs_wtp2550 <- subset(costs_wtp2550, select = -p2550)
  
  costs_wtp5075 <- left_join(avg_costs_wtp5075, quant_costs_wtp5075,
                             by = c("ICD" = "ICD", "agegrp_v3" = "agegrp_v3"))
  costs_wtp5075 <- rename(costs_wtp5075, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
  costs_wtp5075 <- subset(costs_wtp5075, select = -p5075)
  
  costs_wtp75100 <- left_join(avg_costs_wtp75100, quant_costs_wtp75100,
                              by = c("ICD" = "ICD", "agegrp_v3" = "agegrp_v3"))
  costs_wtp75100 <- rename(costs_wtp75100, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)
  costs_wtp75100 <- subset(costs_wtp75100, select = -p75100)
  
  # write out to csv
  wtp025_loc <- paste(output_dir,"costs_p025.csv")
  
  write.csv(costs_wtp025,
            file = wtp025_loc, 
            row.names = FALSE,
            quote = FALSE)
  
  wtp2550_loc <- paste(output_dir,"costs_p2550.csv")
  write.csv(costs_wtp2550,
            file = wtp2550_loc,
            row.names = FALSE,
            quote = FALSE)
  wtp5075_loc <- paste(output_dir,"costs_p5075.csv")
  write.csv(costs_wtp5075,
            file = wtp5075_loc,
            row.names = FALSE,
            quote = FALSE)
  wtp75100_loc <- paste(output_dir,"costs_p75100.csv")
  write.csv(costs_wtp75100,
            file = wtp75100_loc,
            row.names = FALSE,
            quote = FALSE)
  
  
  return(list(costs_export, costs_wtp025, costs_wtp2550, costs_wtp5075, costs_wtp75100))
  
  
  
}
