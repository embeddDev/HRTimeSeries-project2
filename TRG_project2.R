#TASK 1 - project 2
############################################################################################
#course:  Time Series Analysis, T-862-TIMA
#Date:    /09/2014
#Students: Daniel Bergmann Sigtryggsson, Lilja Bjorg Gudmundsdottir, Jon Vilberg Georgsson 
#
############################################################################################
require("forecast")
data = read.csv("veks.csv")
data$HC.f = ts(data$HC.f,frequency = 24, start = c(1995,((data$ds.diy[1]*24)+data$ds.hh[1])))
data$Ta.f = ts(data$Ta.f, frequency = 24, start =c( 1995,((data$ds.diy[1]*24)+data$ds.hh[1])))
data$W.f = ts(data$W.f, frequency = 24, start = c(1995,((data$ds.diy[1]*24)+data$ds.hh[1])))
data$GR.f = ts(data$GR.f, frequency =24, start = c(1995,((data$ds.diy[1]*24)+data$ds.hh[1])))
time = as.POSIXct("1960-1-1") + (data$jdate*24 + data$hh) *3600
# ------------ Task 1 ----------------------
#Consider the time series of heat consumption.
par(mfrow=c(2,1))
plot.ts(data$HC.f ,
     type='l',
     main='Original Heat consumption(GJ/h) data ',
     xlab="running time",
     ylab="Heat consumption [GJ/h]",
     col="red")
grid()
#Comment: we can see that the data set is nonstationary. Thus, we difference the series.

plot.ts(diff(data$HC.f, difference=1) ,
     type='l',
     main='Differenced Heat consumption(GJ/h) data ',
     xlab="running time",
     ylab="Heat consumption [GJ/h]",
     col="red")
grid()

#Estimate the autocovariance, autocorrelation and partial autocorrelation
#functions for the heat consumption.
dev.off()
par(mfrow = c(2,2))
acf(diff(data$HC.f, difference=1))
pacf(diff(data$HC.f, difference=1))
acf(diff(data$HC.f, difference=1), type="covariance")
#Estimate the spectrum for the variations of the heat consumption.
spectrum(diff(data$HC.f, difference=1),
         method = "ar")

#You should comment on the estimated correlation functions
#and spectrum, e.g. what can you tell about the process itself,
#about a good choice of model structure for modelling the process, etc..

#Comment: From the ACF and PACF, we assume we can model the seris as ARIMA(8,2) from table 6.1
#     This is probably too many parameters. Thus we suspect there is seaonality, we need to look at seasonal variations
#     The spectrum reinforces this theory that there are still some variations in higher frequencies.

# -------------------------- Task 2 ---------------------------
# A model for predicting the heat consumption several hours ahead should be formulated,
# and the prediction performance should be analyzed. In particular the prediction 
# performance for a one hour and a six hour prediction shall be described. In this 
# task only the measured heat consumption (and possibly the time) should be used in the model.
# Make a vector for week/weekend days
WD = NULL                        
for(i in 1:length(data$ds.dow)){
  if(data$ds.tod[i] == 1){
    WD[i] = 0
  }
  else{                    
    WD[i] = 1
  }  
}
# making a matrix with the regressors we want to use, do it here, before we cut to training
# and test set.
data$ext_regressors <- cbind(WD=data$WD,Ta.f=data$Ta.f,GR.f=data$GR.f,W.f=data$W.f)
# split the data into training set and test set
training_set = subset(data, data$row.names <= 6000)
test_set = subset(data,data$row.names > 6000)
time_training = time[1:6000]
time_test = time[6001:length(time)]

# Elaborate on your choice of model.  Your predictions should 
# include prediction intervals and an assessment of their quality.

par(mfrow = c(2,1))
acf(diff(data$HC.f, difference=1), lag.max=300)
pacf(diff(data$HC.f, difference=1), lag.max=300)

#Comment: We can see 24hour periodic correlations in the data, and some weekly variations
# We want our model to include 24hour and weekly seasonality 
# We will accomodate weekend seasonality with a vector of external regressors.
# We have made this regressor earlier in data$ext_regressors[,1]
dev.off()
fit1 = Arima(training_set$HC.f, 
             order=c(1,1,1),
             seasonal = list(order = c(1,0,1),period = 24),
             xreg = training_set$ext_regressors[,1] #working days regressors
)
fit1
acf(fit1$residuals,lag.max=50)
#We discard this model, because the residuals are not within the confidence interval
fit2 = Arima(training_set$HC.f, 
             order=c(1,1,2),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = training_set$ext_regressors[,1]  # working days regressor
)
acf(fit2$residuals, lag.max=50)
fit2
#We discard this model, because the residuals are not within the confidence interval

fit3 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = training_set$ext_regressors[,1] #working days regressor
)
acf(fit3$residuals)
fit3
#We like this this model, because the residuals are more or less within the confidence interval and AIC has the smallest value
# and are thus deemed to be white noise.
dev.off()
par(mfrow = c(2,1))
forecast_1ahead = forecast.Arima(fit3,fan=TRUE, h=1, xreg= test_set$ext_regressors[1,1])

plot.forecast(forecast_1ahead,
              include=24, 
              type='l',
              fcol="red",
              col='blue',
              main="Forecast, 1 hour ahead")
lines(c(training_set$HC.f, test_set$HC.f[1]))

forecast_6ahead= forecast.Arima(fit3,fan=TRUE, h=6, xreg= test_set$ext_regressors[1:6,1])
plot.forecast(forecast_6ahead,
              include=24,
              type='l',
              fcol="red",
              col='blue',
              main="Forecast, 6 hour ahead")
lines(c(training_set$HC.f, test_set$HC.f[1:6]))

# EIGUM EFTIR AÐ GERA PREDICTION PERFOMRMANS Á SPÁNNI OG
# META QUALITY Á PREDICTION INTERVAL.

#------------------------Task 3 (15%)----------------------------------------
dev.off()
plot.ts(data$Ta.f,
     main="Ambient air temperature data",
     type='l',
     ylab='Centigrade',
     col='red')
grid()
plot.ts(diff(data$Ta.f,difference=1),
     main="Ambient air temperature data (first difference)",
     type='l',
     ylab='Centigrade',
     col='red')
grid()
dev.off()
ccf(diff(data$HC.f,difference=1),diff(data$Ta.f,difference=1),
    main="cross correlation, HC vs amb temp",
    col='red')
#Clearly, ambient temp shares common trend/seasonality with heatconsumtion data. Thus, we wish to prewhiten the amb temp series.


ambTemp_filtered = Arima(training_set$Ta.f,
                          model=fit3, 
                          fixed = c( 0.9766,  -0.2084, -0.8984,  0.9961,  -0.9224 ,-14.2153),
                         )

#here we have the differences between observed fit3 valuse and estimated, Ta.f values based on the fit3 model
ccf(fit3$residuals, residuals(ambTemp_filtered), na.action=na.omit)
grid()
#We obsere non-zero values at lag 4 and perhaps other lags beyond 4, with no particular pattern
# correlations near 9 at other lags.
#>>> possible regression terms in model: lag(HC.f,4), perhaps additional lags of HC.f
# no lags of Ta.f
HC_lag4 = lag(data$HC.f,4)
data$ext_regressors = cbind(data$ext_regressors,HC_lag4 )


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


#TASK4
dev.off()
#fit4 is the arma model with 2 external regressors, work days and amb temp
fit4 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = training_set$ext_regressors[,1:2])
acf(fit4$residuals)
fit4
#fit 5 is the arma model with 4 external regressors, work-days, amb temperature and wind-speed
fit5 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg= training_set$ext_regressors[,1:3])
acf(fit5$residuals)
fit5
#fit 6 is the arma model with 5 external regressors, work-days, amb temperature, wind-speed and solar radiation.
fit6 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg= training_set$ext_regressors[,1:3])
acf(fit6$residuals)
fit6

#Sources: STAT 510, Penn state,lesson 9.1
#         https://onlinecourses.science.psu.edu/stat510/node/75



## ccf á residuals þá geta se´ð laggið betur
