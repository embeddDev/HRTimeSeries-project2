#TASK 1 - project 2
############################################################################################
#course:  Time Series Analysis, T-862-TIMA
#Date:    /09/2014
#Students: Daniel Bergmann Sigtryggsson, Lilja Bjorg Gudmundsdottir, Jon Vilberg Georgsson 
#
############################################################################################
require("forecast")
data = read.csv("veks.csv")
time = as.POSIXct("1960-1-1") + (data$jdate*24 + data$hh) *3600
#Consider the time series of heat consumption.
par(mfrow=c(2,1))
plot(data$HC.f ~time,
     type='l',
     main='Original Heat consumption(GJ/h) data ',
     xlab="running time",
     ylab="Heat consumption [GJ/h]",
     col="red")
grid()
#Comment: we can see that the data set is nonstationary. Thus, we difference the series.
HC_diff = diff(data$HC.f, difference=1)
plot(HC_diff ~time[2:length(time)],
     type='l',
     main='Differenced Heat consumption(GJ/h) data ',
     xlab="running time",
     ylab="Heat consumption [GJ/h]",
     col="red")
grid()
#Estimate the autocovariance, autocorrelation and partial autocorrelation
#functions for the heat consumption.

WD = NULL
for(i in 1:length(data$ds.dow)){
  if(data$ds.tod[i] == 1){
    WD[i] = 0
  }
  else{
    WD[i] = 1
  }  
}
data$WDa <- WD    #Working days: 0 or 1
training_set = subset(data, data$row.names <= 6000)
time_training = time[1:6000]
test_set = subset(data,data$row.names > 6000)
time_test = time[6001:length(time)]
attach(training_set)
xreg1 = WDa
xreg4 = cbind(WDa,Ta.f)  #passar við fit4
xreg5 = cbind(WDa,Ta.f,GR.f)
xreg6 = cbind(WDa,Ta.f,GR.f,W.f)
detach(training_set)
dev.off()
par(mfrow = c(3,1))
acf(HC_diff)
pacf(HC_diff)
acf(HC_diff, type="covariance")
#Estimate the spectrum for the variations of the heat consumption.
dev.off()
spectrum(HC_diff,
         method = "ar")

#You should comment on the estimated correlation functions
#and spectrum, e.g. what can you tell about the process itself,
#about a good choice of model structure for modelling the process, etc..

#Comment: From the ACF and PACF, we assume we can model the seris as ARIMA(8,2) from table 6.1
#     This is probably too many parameters. Thus we suspect there is seaonality, we need to look at seasonal variations
#     The spectrum reinforces this theory that there are still some variations in higher frequencies.
par(mfrow = c(2,1))
acf(HC_diff, lag.max=300)
pacf(HC_diff, lag.max=300)
#Comment: We can see 24hour periodic correlations in the data, and some weekvariations
# We want our model to include 24hour and weekly seasonality 
# We will accomodate weekend seasonality with a vactor of external regressors.
HC_reg = NULL
for(i in 1:length(training_set$ds.dow)){
  if(training_set$ds.tod[i] == 1){
    HC_reg[i] = 0
  }
  else{
    HC_reg[i] = 1
  }  
}
# test_reg = NULL
# for(i in 1:length(test_set$ds.dow)){
#   if(test_set$ds.tod[i] ==1){  #if working day
#     test_reg[i] = 0
#   }
#   else{
#     test_reg[i] = 1            #if not working day
#   }  
#   
#   
# }
#HC_reg = t(HC_reg) # transpose to a coloumn vector
dev.off()
fit1 = Arima(training_set$HC.f, 
             order=c(1,1,1),
             seasonal = list(order = c(1,0,1),period = 24),
             xreg = HC_reg)
fit1
acf(fit1$residuals,lag.max=50)
#We discard this model, because the residuals are not within the confidence interval
fit2 = Arima(training_set$HC.f, 
             order=c(1,1,2),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = xreg1)
acf(fit2$residuals, lag.max=50)
fit2
#We discard this model, because the residuals are not within the confidence interval

x =[HC_reg,training_set$Ta.f]

fit3 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = xreg1)
acf(fit3$residuals)
fit3

fit4 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = xreg4)
acf(fit4$residuals)
fit4

fit5 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = xreg5)
acf(fit5$residuals)
fit5

fit6 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = xreg6)
acf(fit6$residuals)
fit6

#We like this this model, because the residuals are more or less within the confidence interval and AIC has the smallest value
# and are thus deemed to be white noise.
dev.off()
par(mfrow = c(2,1))
Y1 = forecast.Arima(fit3,fan=TRUE, h=1, xreg= test_reg[1])

plot.forecast(Y1,
              include=24, 
              type='l',
              fcol="red",
              col='blue',
              main="Forecast, 1 hour ahead")
lines(c(training_set$HC.f, test_set$HC.f[1]))
Y2= forecast.Arima(fit3,fan=TRUE, h=6, xreg= test_reg[1:6])
plot.forecast(Y2,
              include=24,
              type='l',
              fcol="red",
              col='blue',
              main="Forecast, 6 hour ahead")
lines(c(training_set$HC.f, test_set$HC.f[1:6]))
#------------------------Task 3 (15%)----------------------------------------
plot(data$Ta.f ~time,
     main="Ambient air temperature data",
     type='l',
     ylab='Centigrade',
     col='red')
grid()
plot(diff(data$Ta.f,difference=1)
     ~time[2:length(data$Ta.f)],
     main="Ambient air temperature data (first difference)",
     type='l',
     ylab='Centigrade',
     col='red')
grid()
plot(diff(data$Ta.f,difference=2)
     ~time[3:length(data$Ta.f)],
     main="Ambient air temperature data (second difference)",
     type='l',
     ylab='Centigrade',
     xlab='time',
     col='red')
grid()
dev.off()
ccf(diff(data$HC.f,difference=1),diff(data$Ta.f,difference=1),
    main="cross correlation, HC vs amb temp",
    col='red')
#Clearly, ambient temp shares common trend/seasonality with heatconsumtion data. Thus, we wish to prewhiten the amb temp series.

Ta.filtered = filter(training_set$Ta.f,
                filter=c(-0.8967,1,-0.8109,-0.244123),
                sides=2)
ccf(fit3$residuals, Ta.filtered,na.action=na.omit)
grid()
#Lag -2,lag -3, lag-4,lag-5
a = cbind(training_set$Ta.f,
          lag2x= lag(fit3,-2),
          lag3x= lag(fit3,-3),
          lag4x= lag(fit3,-4),
          lag5x= lag(fit3,-5))
lm(training_set$Ta.f ~lag2x+lag3x+lag4x+lag5x,
   data=a, na.action=na.omit)

dev.off()
plot(data$W.f ~time,
     main="Wind speed data",
     type='l',
     ylab='m/s',
     col='red')
grid()
par(mfrow=c(2,1))
ccf(data$W.f,diff(data$HC.f,difference=1),  #cross corr, wind speed vs. HC
    col='red',
    main="cross correlation, HC vs wind speed")  
ccf(data$W.f,diff(data$HC.f,difference=1),  #cross corr, wind speed vs. HC
    lag.max=300,
    col='red',
    main="cross correlation, HC vs wind speed")  
dev.off()
plot(data$GR.f,
     main="Global radiation data",
     type='l',
     ylab='W/m2',
     col='red')
grid()
plot(diff(data$GR.f,difference=1),
     main="Global radiation data",
     type='l',
     ylab='W/m2',
     col='red')
grid()
par(mfrow=c(2,1))
ccf(diff(data$GR.f,difference=1),diff(data$HC.f,difference=1),  #cross corr, global radiation vs. HC
    col='red',
    main="cross correlation, HC vs Global radiation")
ccf(diff(data$GR.f,difference=1),diff(data$HC.f,difference=1),  #cross corr, global radiation vs. HC
    lag.max=300,
    col='red',
    main="cross correlation, HC vs Global radiation")
