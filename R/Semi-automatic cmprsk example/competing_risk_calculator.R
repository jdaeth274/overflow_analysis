options(warn = -1)

#Reading data for competing risk analysis in R

mydata<-read.csv(file = "C:/Users/kryst/Desktop/cr_age.csv")

  str(mydata)
  
### create factor
  
  mydata$category <-factor(mydata$cat, levels=c(0,1,2), labels=c("<25", "25-64","65+"))
  mydata$status_f <-factor(mydata$status, levels=c(1,2), 
                                labels=c("Survived","Death"))
  mydata$wt <- sample(100, size = nrow(mydata), replace = TRUE)
  
## trial crr run
  library(cmprsk)
  attach(mydata)
  cr1 <- crr(ftime = ftime, fstatus = status, cov1 = wt)
  predict.crr(cr1,cov1=0)
  test <- do.call(rbind.data.frame, p)
  plot.predict.crr(p,lty=1:(ncol(p)-1), color=1,  
                   ylim=c(0, max(p[, -1])), xmin=0, xmax=max(p[, 1]))
  
  library("timereg")
  cr2 <- comp.risk(Event(ftime, status)~wt,data=mydata, cause = 1)

## basic descriptive summaries
  attach(mydata)
  
  table(category, status_f)

  round(tapply(ftime, list(category, status_f), mean), digits = 2)
  
  # Interactive table
  
  library(rpivotTable)
  rpivotTable(data=mydata, rows="category",cols = "status_f",
            vals = "ftime", aggregatorName = "Average",rendererName = "Table" )
  
# Fitting the cumulative incidence function
# The function 'CumIncidence' will output: plots, statistic, time-point estimates, 
#                                          standard errors and 95% CI (if specified)
  
  source ("/Users/Monkeyface/Dropbox/00_COVID/Rcode/Overflows/Tutorial/CumIncidence.R")

  col<-c("darkblue","darkred","darkgreen")
  
  fit<-CumIncidence (mydata$ftime, mydata$status_f, mydata$category, 
                     cencode = 0, 
                     xlab = "Days",
                     col=col)  
  
  # If you want to plot specific time steps only
  fit<-CumIncidence(mydata$ftime, mydata$status_f, mydata$category, cencode=0,
                    xlab = "Days",
                    t=c(0, 1, 2, 3, 4, 5, 7, 10, 14, 21, 40, 60),
                    col=col)

  # Add 95% confidence interval and sepparate into pannels
  
  par(mfrow=c(2,3))
  
  fit<-CumIncidence(mydata$ftime, mydata$status_f, mydata$category, cencode=0,
                    xlab = "Months",
                    t=c(0, 1, 2, 3, 4, 5, 7, 10, 14, 21, 40, 60),
                    col=col,level=0.95)
  
  
  par(mfrow=c(1,1))
  
  
  