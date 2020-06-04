###############################################################################
## Master Overflows analysis script ###########################################
###############################################################################

require(vroom)
require(dplyr)
require(plyr)
require(ggplot2)
require(ggpubr)
require(snow)
require(tidyverse)
require(lubridate)
require(KFAS)
require(forecast)



## First step to load up HES dataset ##

input_args <- commandArgs(trailingOnly = TRUE)

hes_data <- vroom(input_args[1], delim = ",")

cohort_allocation <- cohort_set_up(num_cores = 16)
                                   


###############################################################################
## Function list here: 

source(file = "C:/Users/dr314/Dropbox/COVID19/Overflow/Rscript/cohort_identification.R")
source(file = "C:/Users/dr314/Dropbox/COVID19/Overflow/Rscript/time_series_creator.R")  
