###############################################################################
## Costing analysis ###########################################################
###############################################################################
require(stringr)
require(dplyr)
require(vroom)
require(dtplyr)
require(data.table)

# Organisational level unit costs
system.time(costs_dataset <- vroom("D:/NHS_Costs/merged_2009_2019_orgrcs.csv", delim = ","))

# HES data
system.time(transitions_data_whole <- fread("E:/HES/COVID/HES_APC_CC_0913_transitions_all_ICD.csv",
                                            ))
#load ref cost national data for 1213 - for costs that are unmatched using organisational level cost data
national_rcs1213=read.csv("D:/Dropbox/COVID19/Overflow/Data/national_rcs1213.csv",fileEncoding = "UTF-8-BOM")

# load national ref cost data for 1112 - for costs that are unmatched after using national_rcs1213
national_rcs1112 = read.csv("D:/Dropbox/COVID19/Overflow/Data/national_rcs1112.csv",fileEncoding = "UTF-8-BOM")

# load national ref cost data for 1011 - for costs that are unmatched after using national_rcs1112
national_rcs1011 = read.csv("D:/Dropbox/COVID19/Overflow/Data/national_rcs1011.csv",fileEncoding = "UTF-8-BOM")

# load national ref cost data for 0910 - for costs that are unmatched after using national_rcs1011
national_rcs0910 = read.csv("D:/Dropbox/COVID19/Overflow/Data/national_rcs0910.csv",fileEncoding = "UTF-8-BOM")

#ORG.CODE matched to procode3 
#Currency.code matched to HRG or susHRG (thisi s for post FY1213 onwards. will need to include a change in variable name for HF pre cleaning file)
#we need to keep the UNIT COST, MEAN COST, Actual Cost MAPPING POT and PERIOD (2009-10) EpiEnd_FY== "0910"
#Mean cost = the national mean average unit cost
#unit cost = the average cost to the organisation of providing that activity
#actual cost = the organisation's activity (number of cases) x the organisation unit cost
#we are really only interested in matching the unit cost of each organisation's hrg
#Mapping pot is matched to admimeth_C
#Mapping to 01_EI and 02_NEI and 03_XS
#XS bedays with the bed.days variable


###############################################################################
## Data is loaded up, lets merge these costs! #################################
###############################################################################

# first, only keep FY1213 in transitions_data_whole as it is the only year with SUSHRGep code (cost dataset has HRG as HRG4.0, imperial hes has hrg 3.5)
transitions_data_1213 <- subset(transitions_data_whole, EpiEnd_FY == "1213")

## alter costs data to include epiEnd_FY version of PERIOd & admimeth_version of Mapping Pot - create EpiEnd_FY as an integer not charc to match the HES

costs_period <- costs_dataset$PERIOD
costs_sub <- gsub("20","",costs_period)
costs_sub <- gsub("-","",costs_sub)

costs_dataset$EpiEnd_FY <- costs_sub

## admimeth version 


costs_dataset$admimeth_c <- 0

elective_rows <- which(costs_dataset$MAPPING.POT == "01_EI")
emerg_rows <- which(costs_dataset$MAPPING.POT == "02_NEI")
extra_days <- which(costs_dataset$MAPPING.POT == "03_XS")

costs_dataset$admimeth_c[elective_rows] <- 1
costs_dataset$admimeth_c[emerg_rows] <- 2


###############################################################################
## Mapping the codes across ###################################################
###############################################################################

system.time(narrowed_costs <- costs_dataset[which(costs_dataset$admimeth_c != 0),], )

narrowed_costs <- as.data.frame(narrowed_costs)

## limit to current epiendfy 

fy_epiends <- as.integer(plyr::count(transitions_data_1213$EpiEnd_FY)[,1])

narrowed_costs <- narrowed_costs[narrowed_costs$EpiEnd_FY %in% fy_epiends,] #check this maybe replace with row 66

narrowed_costs$EpiEnd_FY <- as.integer(narrowed_costs$EpiEnd_FY)

###############################################################################
## Loop through the narrowed costs and match into the HES data ################
###############################################################################

only_int_cols <- narrowed_costs[,colnames(narrowed_costs) %in% c("ORG.CODE","CURRENCY.CODE","UNIT.COST","MEAN",
                                                                 "ACTUAL.COST","admimeth_c","EpiEnd_FY",
                                                                 "MAPPING.POT")]


rm(costs_dataset)
rm(narrowed_costs)
rm(emerg_rows)
rm(costs_sub)
rm(elective_rows)
rm(extra_days)
rm(costs_period)
rm(fy_epiends)
gc()



hacky_merging <- function(groupings_data,hes_data, costs_data){
  
  
  if(is.data.table(hes_data) == FALSE)
    hes_data <- as.data.table(hes_data)
  
  

  current_admimeth <- str_split_fixed(groupings_data, "-",3)[1]
  current_year <- str_split_fixed(groupings_data, "-",3)[2]
  current_half <- as.integer(str_split_fixed(groupings_data, "-",3)[3])
  
  hes_data <- hes_data[hes_data$admimeth_C == as.integer(current_admimeth) &
                         hes_data$EpiEnd_FY == current_year]
  cost_row <- which(costs_data$EpiEnd_FY == current_year &
                      costs_data$admimeth_c == as.integer(current_admimeth))
  costs_data <- costs_data[cost_row,]
  
  procode_currency_grp <- unique(paste(hes_data$procode3, hes_data$SUSHRGep, sep = "-")) #need to check if HRG_35 is also used, the SUSHRGep is only for 1213
  
  split_point <- length(procode_currency_grp) / 2
  
  if(current_half == 1)
    procode_split <- procode_currency_grp[1:split_point]
  else
    procode_split <- procode_currency_grp[(split_point + 1):length(procode_currency_grp)]
  
  
  print("Starting loop now")
  for(k in 1:length(procode_split)){
    current_group <- procode_split[k]
    current_split <- str_split_fixed(current_group, "-", 2)
    current_procode <- current_split[1]
    current_SUSHRG <- current_split[2]
    
    cost_row <- which(costs_data$ORG.CODE == current_procode &
                        costs_data$CURRENCY.CODE == current_SUSHRG)
                        
    if(length(cost_row) >= 1){
      hes_rows <- hes_data[hes_data$procode3 == current_procode &
                             hes_data$SUSHRGep == current_SUSHRG]
      
      hes_rows <- hes_rows$index_2
      
      hes_data$unit_cost[hes_data$index_2 == hes_rows] <- costs_data$UNIT.COST[cost_row[1]]
      hes_data$mean_cost[hes_data$index_2 == hes_rows] <- costs_data$MEAN[cost_row[1]]
      hes_data$actual_cost[hes_data$index_2 == hes_rows] <- costs_data$ACTUAL.COST[cost_row[1]]
      multi_row <- ifelse(length(cost_row) > 1, "Yes_multi", "Yes")
      hes_data$altered[hes_data$index_2 == hes_rows] <- multi_row
      
      
    }
    
    
  }
  
  return(hes_data)
  
}


hacky_merging_cluster <- function(hes_data, costs_data, num_core){
  
  ## This function sets up the data to run in parrallel 
  
  require(snow)
  require(tictoc)
  print("adding extra cols to hes_data")
  tic("adding extra cols")
  hes_data$unit_cost <- NA
  hes_data$mean_cost <- NA
  hes_data$actual_cost <- NA
  hes_data$altered <- "No"
  hes_data$index_2 <- seq(1, nrow(hes_data))
  toc()
  print("Narrowing HES data")
  tic("narrowing data")
  hes_data <- hes_data[,c("procode3","admimeth_C","SUSHRGep","EpiEnd_FY",
                          "index_2","unit_cost","mean_cost","actual_cost",
                          "altered")]
  toc()
  
  print("Getting unique combos")
  tic("year_admi combos")
  year_admi_groupings <- rep(unique(paste(hes_data$admimeth_C, hes_data$EpiEnd_FY, sep = "-")), 2)
  halves <- rep(c("1","2"), each = length(year_admi_groupings)/2)
  year_admi_halves <- paste(year_admi_groupings, halves, sep = "-")
  toc()
  
  
  ### So we run each of the admimeth & Financial years on the 
  
  # electives_data <- hes_data[hes_data$admimeth_C == 1,]
  # groupings_elec <- paste(electives_data$procode3, electives_data$admimeth_C, 
  #                         electives_data$SUSHRGep, electives_data$EpiEnd_FY,
  #                    sep = "-")
  # 
  # groupings_unique_elec <- unique(groupings_elec)
  # 
  
  ## set up the cluster 
  print("Setting up electives cluster")
  electives_cluster <- snow::makeCluster(spec = length(year_admi_halves))
  electives_dat <- snow::clusterSplit(electives_cluster, year_admi_halves)
  
  tic("Requiring packages")
  snow::clusterEvalQ(electives_cluster, require(data.table))
  snow::clusterEvalQ(electives_cluster, require(stringr))
  toc()
  print("Copying data")
  tic("Copying over electives data")
  snow::clusterExport(electives_cluster, "hes_data", envir = environment())
  snow::clusterExport(electives_cluster, "costs_data", envir = environment())
  toc()
  tic("COpying over function")
  snow::clusterExport(electives_cluster, "hacky_merging")
  toc()
  print("Running electives on the cluster now")
  tic("Running electives")
  electives_merging <-snow::clusterApply(electives_cluster, electives_dat,
                                                fun = hacky_merging,
                                                hes_data = hes_data,
                                                costs_data = costs_data)
  snow::clusterEvalQ(electives_cluster, rm(electives_data))
  snow::clusterEvalQ(electives_cluster, rm(costs_data))
  snow::clusterEvalQ(electives_cluster, gc())
  snow::stopCluster(electives_cluster)
  toc()
  electives_out <- dplyr::bind_rows(electives_merging)
  
  ## Now for the emergencies
  
    
  return(electives_out)
  
}

transitions_data_costs <- hacky_merging_cluster(transitions_data_1213, only_int_cols)

hes_data <- transitions_data_0912[,c("procode3","admimeth_C","SUSHRGep","EpiEnd_FY")]



merging_left_func <- function(hes_data, costs_data, loop_iter){

  nrow_splits <- seq(1,nrow(transitions_data_whole), length.out = loop_iter + 1 )
  hes_data <- hes_data[,c("procode3","admimeth_C","EpiEnd_FY","SUSHRGep")]
  
  merged_dat <- NULL

  for(k in 1:loop_iter){
    
    print(k)
    current_trans_data <- hes_data[nrow_splits[k]:nrow_splits[k + 1],]
    
    start_time <- Sys.time()
    
    transitions_data_whole_costs <- dplyr::left_join(current_trans_data, y = costs_data,
                                                   by = c("procode3" = "ORG.CODE",
                                                          "admimeth_C" = "admimeth_c",
                                                          "EpiEnd_FY" = "EpiENd_FY",
                                                          "SUSHRGep" = "CURRENCY.CODE"))
    merged_dat <- rbind.data.frame(merged_dat, transitions_data_whole_costs)
    end_time <- Sys.time()
    print(end_time - start_time)
    
  }

  return(merged_dat)
}




transitions_data_merged <- merging_left_func(transitions_data_whole, only_int_cols, 12)




transitions_data_whole_procode <- dplyr::left_join(transitions_data_whole, only_int_cols,
                                                   by = c("procode3 = admimeth_c"))


## Lets try this with the dtplyr package

lazy_trans_whole <- transitions_data_whole %>% dtplyr::lazy_dt()
lazy_costs_whole <- only_int_cols %>% dtplyr::lazy_dt()

joined_whole <- left_join(x = lazy_trans_whole, y = lazy_costs_whole,
                          by = c("procode3" = "ORG.CODE",
                                 "admimeth_C" = "admimeth_c",
                                 "EpiEnd_FY" = "EpiENd_FY"))




transitions_data_whole_costs <- dplyr::left_join(trans_whole, y = only_int_cols,
                                                 by = c("procode3 = ORG.CODE",
                                                        "admimeth_C = admimeth_c",
                                                        "EpiEnd_FY = EpiENd_FY"))
transitions_data_whole_procode <- dplyr::left_join(transitions_data_whole, only_int_cols,
                                                   by = c("procode3" = "ORG.CODE"))

transitions_data_whole$procode3

## try out dplyr left join 

test_df <- data.frame(matrix(ncol = 4,nrow = 5))
test_df_2 <- data.frame(matrix(ncol = 3, nrow = 10))

colnames(test_df) <- c("id","cost","sum","mean")
colnames(test_df_2) <- c("id","cost","mode")


test_df[,1] <- seq(1,5)
test_df_2[,1] <- seq(1,length.out = 5,by = 2)

test_df[,2] <- seq(1.1,5.1, by = 1)
test_df_2[,2] <- seq(0.1,length.out = 10,by = 1)

test_df[,3] <- seq(2.2,6.2, by = 1)
test_df_2[,3] <- seq(100,length.out = 10,by = 2)

test_df[,4] <- seq(3.3,7.3, by = 1)

test_left_join <- dplyr::left_join(test_df, test_df_2, by = c("id", "cost"))
test_left_join

## That looks reasonable lets try it out with the narrowed costs data YAY!!!!!

####################################################################################################
##### START OF COST MERGING DONE BY KL AND DR...JOSH WHAT HAVE YOU DONE ABOVE THIS???
mergedcosts_orglevel <- left_join(transitions_data_1213, only_int_cols, 
                         by = c("procode3" = "ORG.CODE",
                               "admimeth_C" = "admimeth_c",
                               "EpiEnd_FY" = "EpiEnd_FY",
                               "SUSHRGep" = "CURRENCY.CODE"))

# subset episodes who have UNIT.COST in the mergedcosts dataset (the unit cost are org-hrg level)
notNA_merged_orglevel <- mergedcosts_orglevel[!is.na(mergedcosts_orglevel$UNIT.COST),]

# subset episodes with UNIT.COST as NA for mergedcosts_orglevel (procode3-SUSHRGep pair does not exist in costs data but does exist in hes data)
NA_merged_orglevel <- mergedcosts_orglevel[is.na(mergedcosts_orglevel$UNIT.COST),]

# remove UNIT.COST column from unmerged_orglevel
NA_merged_orglevel <- select(NA_merged_orglevel, -UNIT.COST)

# match national reference cost schedule data to the unmerged episodes
merged_rcs1213level <- left_join(NA_merged_orglevel, national_rcs1213,
                          by = c("admimeth_C" = "admimeth_C",
                                 "SUSHRGep" = "CURRENCY.CODE"))

# subset episdoes with and without UNIT.COST as NA for merged_rcs1213level
notNA_merged_rcs1213level <- merged_rcs1213level[!is.na(merged_rcs1213level$UNIT.COST),]
NA_merged_rcs1213level <- merged_rcs1213level[is.na(merged_rcs1213level$UNIT.COST),]

# remove UNIT.COST column from unmerged_rcs1213level
NA_merged_rcs1213level <- select(NA_merged_rcs1213level, -UNIT.COST)

# match unmerged_rcs1213level with national_rcs1112
merged_rcs1112level <- left_join(NA_merged_rcs1213level, national_rcs1112,
                                 by = c("admimeth_C" = "admimeth_C",
                                 "SUSHRGep" = "CURRENCY.CODE"))

# subset episodes with and without UNIT.COST as NA for merged_rcs1112level
notNA_merged_rcs1112level <- merged_rcs1112level[!is.na(merged_rcs1112level$UNIT.COST),]
NA_merged_rcs1112level <- merged_rcs1112level[is.na(merged_rcs1112level$UNIT.COST),]

# remove UNIT.COST column
NA_merged_rcs1112level <- select(NA_merged_rcs1112level, -UNIT.COST)

# match unmerged_rcs1112level with national_rcs1011
merged_rcs1011level <- left_join(NA_merged_rcs1112level, national_rcs1011,
                                 by = c("admimeth_C" = "admimeth_C",
                                 "SUSHRGep" = "CURRENCY.CODE"))

notNA_merged_rcs1011level <- merged_rcs1011level[!is.na(merged_rcs1011level$UNIT.COST),]
NA_merged_rcs1011level <- merged_rcs1011level[is.na(merged_rcs1011level$UNIT.COST),]
NA_merged_rcs1011level <- select(NA_merged_rcs1011level, - UNIT.COST)

merged_rcs0910level <- left_join(NA_merged_rcs1011level, national_rcs0910,
                                 by = c("admimeth_C" = "admimeth_C",
                                 "SUSHRGep" = "CURRENCY.CODE"))

# append merged_rcs0910level to merged_orglevel to form complete hes-cost merge (still missing 16,635 hrg-cost matches, 0.27% of HES data - this is as good as it gets)
full_mergedcosts <- rbind(merged_orglevel,
                          notNA_merged_rcs1213level,
                          notNA_merged_rcs1112level,
                          notNA_merged_rcs1011level,
                          merged_rcs0910level)

# drop NAs from full_mergedcosts$UNIT.COST
full_mergedcosts_noNA <- full_mergedcosts[!is.na(full_mergedcosts$UNIT.COST),]
  
# calculate average unit cost per patient group (mean, sd, min, p25, p50, p75, max)
quantilecost <- full_mergedcosts_noNA %>% group_by(ICD, agegrp_v3, admimeth_C) %>%
  do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

avgcost <- do.call(data.frame, aggregate(full_mergedcosts_noNA$UNIT.COST, 
                                         list(full_mergedcosts_noNA$ICD, full_mergedcosts_noNA$agegrp_v3, full_mergedcosts_noNA$admimeth_C),
                                         function(x) c(mean = mean(x), sd = sd(x))))

colnames(avgcost) <- c("ICD","agegrp","admimeth","avg_cost","sd_cost")

costs <- left_join(avgcost, quantilecost,
                   by = c("ICD" = "ICD",
                          "agegrp" = "agegrp_v3",
                          "admimeth" = "admimeth_C"))
costs <- rename(costs, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)

costs$disease <- "ICD"

costs$age <- paste("_AGE",costs$agegrp)

costs$p <- paste0(costs$disease, costs$ICD,costs$age)
costs$S <- "G"
costs$a <- ifelse(costs$admimeth == 1, "N", "E")

costs_export <- select(costs, -c(ICD, agegrp, admimeth, disease, age))

write.csv(costs_export,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/costs.csv",
          row.names = FALSE,
          quote = FALSE)

# calculate average unit cost per patient group across WT for electives (eg. for p25 WT, what is the cost p25, p50, p75, mean, sd)
electives_merged <- subset(merged, admimeth_C == 1)

quantilewt <- electives_merged %>% group_by(ICD, agegrp_v3) %>%
  do(data.frame(t(quantile(.$WaitingTime, props = c(0.25, 0.50, 0.75)))))

quantilewt <- rename(quantilewt, min = X0., p25 = X25., median = X50., p75 = X75., max = X100.)

# at each quantile, what is the avg, sd, min, 25, 50, 75, max cost?

# identify which electives fall within each quantile
electives_wtquant <- left_join(electives_merged, quantilewt,
                               by = c("ICD","agegrp_v3"))

electives_wtquant$p025 <- ifelse(electives_wtquant$WaitingTime < electives_wtquant$p25, 1, 0)
electives_wtquant$p2550 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$p25 & 
                                     electives_wtquant$WaitingTime < electives_wtquant$median), 1, 0)
electives_wtquant$p5075 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$median & 
                                     electives_wtquant$WaitingTime < electives_wtquant$p75), 1, 0)
electives_wtquant$p75100 <- ifelse((electives_wtquant$WaitingTime >= electives_wtquant$p75), 1, 0)

electives_p025 <- subset(electives_wtquant, p025 == 1)
electives_p2550 <- subset(electives_wtquant, p2550 == 1)
electives_p5075 <- subset(electives_wtquant, p5075 == 1)
electives_p75100 <- subset(electives_wtquant, p75100 == 1)

# quantiles of costs within each WT quantile
quant_costs_wtp025 <- electives_p025 %>% group_by(ICD, agegrp_v3, p025) %>%
    do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

quant_costs_wtp2550 <- electives_p2550 %>% group_by(ICD, agegrp_v3, p2550) %>%
    do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

quant_costs_wtp5075 <- electives_p5075 %>% group_by(ICD, agegrp_v3, p5075) %>%
    do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

quant_costs_wtp75100 <- electives_p75100 %>% group_by(ICD, agegrp_v3, p75100) %>%
    do(data.frame(t(quantile(.$UNIT.COST, props = c(0.25, 0.50, 0.75)))))

# mean and sd cost within each WT quantile
avg_costs_wtp025 <- do.call(data.frame, aggregate(electives_p025$UNIT.COST, 
                                                  list(electives_p025$ICD, electives_p025$agegrp_v3),
                                        function(x) c(mean = mean(x), sd = sd(x))))
colnames(avg_costs_wtp025) <- c("ICD","agegrp_v3","avg_cost","sd_cost")

avg_costs_wtp2550 <- do.call(data.frame, aggregate(electives_p2550$UNIT.COST, 
                                                  list(electives_p2550$ICD, electives_p2550$agegrp_v3),
                                        function(x) c(mean = mean(x), sd = sd(x))))
colnames(avg_costs_wtp2550) <- c("ICD","agegrp_v3","avg_cost","sd_cost")

avg_costs_wtp5075 <- do.call(data.frame, aggregate(electives_p5075$UNIT.COST, 
                                                  list(electives_p5075$ICD, electives_p5075$agegrp_v3),
                                        function(x) c(mean = mean(x), sd = sd(x))))
colnames(avg_costs_wtp5075) <- c("ICD","agegrp_v3","avg_cost","sd_cost")

avg_costs_wtp75100 <- do.call(data.frame, aggregate(electives_p75100$UNIT.COST, 
                                                  list(electives_p75100$ICD, electives_p75100$agegrp_v3),
                                        function(x) c(mean = mean(x), sd = sd(x))))
colnames(avg_costs_wtp75100) <- c("ICD","agegrp_v3","avg_cost","sd_cost")

# combine mean, sd, and quantile costs at each WT quantile
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
write.csv(costs_wtp025,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/costs_p025.csv", 
          row.names = FALSE,
          quote = FALSE)

write.csv(costs_wtp2550,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/costs_p2550.csv", 
          row.names = FALSE,
          quote = FALSE)

write.csv(costs_wtp5075,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/costs_p5075.csv", 
          row.names = FALSE,
          quote = FALSE)

write.csv(costs_wtp75100,
          file = "D:/Dropbox/COVID19/Overflow/Wolfram/costs_p75100.csv", 
          row.names = FALSE,
          quote = FALSE)














