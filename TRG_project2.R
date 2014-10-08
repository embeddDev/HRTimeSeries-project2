#TASK 1 - project 2
############################################################################################
#course:  Time Series Analysis, T-862-TIMA
#Date:    /09/2014
#Students: Daniel Bergmann Sigtryggsson, Lilja Bjorg Gudmundsdottir, Jon Vilberg Georgsson 
#
############################################################################################
require("forecast")
install.packages("ggplot2")
require("ggplot2")
data = read.csv("veks.csv")
#data$HC.f = ts(data$HC.f,frequency = 24, start = c(1995,((data$ds.diy[1]*24)+data$ds.hh[1])))
#data$Ta.f = ts(data$Ta.f, frequency = 24, start =c( 1995,((data$ds.diy[1]*24)+data$ds.hh[1])))
#data$W.f = ts(data$W.f, frequency = 24, start = c(1995,((data$ds.diy[1]*24)+data$ds.hh[1])))
#data$GR.f = ts(data$GR.f, frequency =24, start = c(1995,((data$ds.diy[1]*24)+data$ds.hh[1])))
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
data$ext_regressors =NULL
data$ext_regressors <- cbind(WD=WD)
# split the data into training set and test set
training_set = subset(data, data$row.names <= 6000)
test_set = subset(data,data$row.names > 6000)
time_training = time[1:6000]
time_test = time[6001:length(time)]

par(mfrow = c(2,1))
acf(diff(data$HC.f, difference=1), lag.max=300)
pacf(diff(data$HC.f, difference=1), lag.max=300)
# Elaborate on your choice of model.  Your predictions should 
# include prediction intervals and an assessment of their quality.

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
fit3 = NULL
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


# META QUALITY Á PREDICTION INTERVAL.
last <- function(x) { tail(x, n = 1) }
compare.q = seq(50,99,by=1)
sigma.eps_1ahead = sd(last(forecast_1ahead$residuals))
sigma.eps_6ahead = sd(tail(forecast_6ahead$residuals,n=6))
obs.p_1ahead = c()
obs.p_6ahead = c()
for(ii in compare.q){
  obs.p_1ahead <- c(obs.p_1ahead,mean(pnorm(last(forecast_1ahead$residuals),sd=sigma.eps_1ahead) <= ii))
  obs.p_6ahead <- c(obs.p_1ahead,mean(pnorm(tail(forecast_6ahead$residuals,n=6),sd=sigma.eps_6ahead) <= ii))
}

dat.obs_1ahead <- data.frame(Observed=obs.p_1ahead,Theoretical=compare.q)
dat.theo_1ahead <- data.frame(Observed=compare.q,Theoretical=compare.q)
dat.obs_6ahead <- data.frame(Observed=obs.p_6ahead,Theoretical=compare.q)
dat.theo_6ahead <- data.frame(Observed=compare.q,Theoretical=compare.q)

ggplot_object = ggplot(data=dat.obs_1ahead )
ggplot_object+
  geom_line(data=dat.theo_1ahead,mapping=aes(x=Theoretical,y=Observed),size=2)+
  geom_point(data=dat.obs_1ahead,mapping=aes(x=Theoretical,y=Observed),size=5,colour="blue")

ggplot_object = ggplot(data=dat.obs_6ahead )
ggplot_object+
  geom_line(data=dat.theo_1ahead,mapping=aes(x=Theoretical,y=Observed),size=2)+
  geom_point(data=dat.obs_6ahead,mapping=aes(x=Theoretical,y=Observed),size=5,colour="blue")

#---------------------------------------------------------------



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
acf(diff(training_set$Ta.f))
#We make new model for external regressor Ta.f
temp.mdl <- arima(training_set$Ta.f,order=c(2,1,0),seasonal=list(order=c(1,0,0),period=24))
#Now we use this model to prewhiten our data (heat consumtion)
hct.res <- residuals(arima(training_set$HC.f,order=c(2,1,0),seasonal=list(order=c(1,0,0),period=24),fixed=temp.mdl$coef))


ccf(hct.res,temp.mdl$residuals)
grid()
Ta_lag1 = lag(training_set$Ta.f,1)
training_set$ext_regressors = cbind(WD=WD[1:6000],Ta_lag1=Ta_lag1) 


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

#Clearly, Wind data shares common trend/seasonality with heatconsumtion data. Thus, we wish to prewhiten the wind series.
par(mfrow=c(2,1))
acf(diff(training_set$W.f))
pacf(diff(training_set$W.f))
#We make new model for external regressorb W.f
temp.mdl <- arima(training_set$W.f,order=c(2,0,0))
#Now we use this model to prewhiten our data (heat consumtion)
HC_Wind.res <- residuals(arima(training_set$HC.f,order=c(2,1,0),seasonal=list(order=c(1,0,0),period=24),fixed=c(temp.mdl$coef)))
HC_Wind.res = na.action(HC_Wind.res,na.rm=TRUE)
dev.off()
ccf(HC_Wind.res,temp.mdl$residuals)  #FUCKED LAG
grid()
W_lag = lag(training_set$W.f,2)
training_set$ext_regressors = cbind(WD=WD[1:6000],Ta_lag1=Ta_lag1, W_lag=W_lag)


plot(data$GR.f ~time,
     main="Global radiation data",
     type='l',
     ylab='W/m2',
     col='red')
grid()
plot(diff(data$GR.f,difference=1) ~time[1:(length(data$GR.f)-1)],
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
dev.off()
#Clearly, solar radiation data shares common trend/seasonality with heatconsumtion data. Thus, we wish to prewhiten the Radiation series.
par(mfrow=c(2,1))
acf(diff(training_set$GR.f))
pacf(diff(training_set$GR.f))
#We make new model for external regressor GR.f
temp.mdl <- arima(training_set$GR.f,order=c(1,1,0),seasonal=list(order=c(1,0,0),period=24))
#Now we use this model to prewhiten our data (heat consumtion)
HC_GR.res <- residuals(arima(training_set$HC.f,order=c(1,1,0),seasonal=list(order=c(1,0,0),period=24),fixed=temp.mdl$coef))

dev.off()
ccf(HC_GR.res,temp.mdl$residuals)  #FUCKED LAG
grid()
GR_lag = lag(training_set$W.f,1)
training_set$ext_regressors = cbind(WD=WD[1:6000],Ta_lag1=Ta_lag1, W_lag=W_lag,GR_lag=GR_lag)



#TASK4
dev.off()
#fit4 is the arma model with 2 external regressors, work days and amb temp
fit4 = Arima(training_set$HC.f,
             order=c(2,1,1),
             seasonal = list(order = c(1,0,1), period = 24),
             xreg = training_set$ext_regressors[,1:4])
acf(fit4$residuals)
fit4



#TASK 5