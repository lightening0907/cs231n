
setwd("/Users/ChiYuan/Documents/learning/analysis/orbital/")
car_data = read.csv("data.csv",header=TRUE)
plot(car_data$date,car_data$car.count,ylab = 'car count',xy.lines=TRUE,type = "l")
lines(car_data$date,car_data$car.count)

car_data$dayofYear = strftime(car_data$date,format='%j')
agg_by_day = aggregate(car_data$car.count,by=list(Category=car_data$dayofYear),FUN=mean)
plot(agg_by_day$Category,agg_by_day$x,xlab='Day of Year',ylab='car count',xy.lines=TRUE,type="l")
lines(agg_by_day$Category,agg_by_day$x)
str(car_data)
car_data[is.na(car_data$car.count),]
# install.packages("zoo")
library(zoo)
# duplicated(car_data$date)
car_ts <- zoo(car_data$car.count,order.by=as.Date(as.character(car_data$date),format='%Y-%m-%d'))
car_ts <- ts(car_ts,frequency=365,start=2010)

# install.packages("imputeTS")
library(imputeTS)
car_ts_ip = na.interpolation(car_ts,option='spline')

acf(car_data$car.count,na.action = na.pass,main="car count")
pacf(car_data$car.count,na.action = na.pass,main="car count")
decom_car_ts <- decompose(car_ts_ip,type='additive')
plot(decom_car_ts)
decom_car_ts$random
acf(decom_car_ts$random)
#complete.cases(decom_car_ts$seasonal)
#nrow(decom_car_ts$seasonal[complete.cases(decom_car_ts$seasonal)])


pacf(decom_car_ts$random,na.action=na.pass)
stl(car_ts,"periodic",na.action=na.pass)

n_ele = length(decom_car_ts$seasonal)
FF = abs(fft(decom_car_ts$seasonal)/sqrt(n_ele))^2
P = (4/n_ele)*FF[1:(n_ele%/%2)]
f = (0:(n_ele%/%2-1))/n_ele
plot(f,P,type="l")
period = 1/f[P==max(P)]

trend = time(car_ts) - mean(time(car_ts))
trend2 = trend^2
regtrend = lm(car_ts_ip~trend+trend2)
summary(regtrend)

plot(as.numeric(trend),as.numeric(car_ts_ip),ylim=c(0,270),xlab="Time", ylab="Car Count")
lines(as.numeric(trend),as.numeric(predict(regtrend)),type="l",col='red',pch=10)

n_ele = length(residuals(regtrend))
FF = abs(fft(residuals(regtrend))/sqrt(n_ele))^2
P = (4/n_ele)*FF[1:(n_ele%/%2)]
f = (0:(n_ele%/%2-1))/n_ele
plot(f,P,type="l")
period = 1/f[P==max(P)]
period

library("astsa")
acf2(residuals(regtrend),main='residual after de-trend')
acf2(diff(residuals(regtrend),91,1))
de_period = diff(residuals(regtrend),91,1)
acf2(de_period)
pacf(de_period)
pacf(diff(de_period,91,1),lag.max=100)
adj_reg = sarima(car_ts_ip,p=0,d=0,q=0,D=1,S=91,xreg=cbind(trend,trend2))

adj_reg
forc_sarima = sarima.for(residuals(regtrend),30,p=0,d=0,q=0,D=1,S=91)
plot(adj_reg)
#?sarima

summary(adj_reg)
predict(adj_reg)
Box.test(de_period,60)
acf(de_period)
pacf(de_period,lag.max=100)

diff(residuals(regtrend),91,1)
ts_resi_detrend = ts(diff(residuals(regtrend),91,1),frequency = 365)
de_ts_resi_detrend = decompose(ts_resi_detrend)
plot(de_ts_resi_detrend)


library(forecast)
nativeR_sarima = arima(car_ts_ip[1:365], order = c(0,0,0), seasonal = list(order = c(0,1,0), period = 91),xreg=cbind(trend[1:365],trend2[1:365]))
# nativeR_sarima

nativeR_predict = predict(nativeR_sarima, newxreg = cbind(trend[366:380],trend2[366:380]),n.ahead=15)
accuracy(nativeR_predict$pred,as.numeric(car_ts_ip[366:380]))
accuracy(nativeR_sarima)


y_wq <- msts(residuals(regtrend), seasonal.periods=c(7,91.3))
fit_wq <- tbats(y_wq)
fc_wq<-forecast(fit_wq)
plot(fc_wq)
summary(fit_wq)
Box.test(fit_wq$errors,30,type = "Ljung")

y_q <- msts(residuals(regtrend), seasonal.periods=c(91.3))
fit_q <- tbats(y_q)
fc_q<-forecast(fit_q)
plot(fc_q)
summary(fit_q)
Box.test(fit_q$errors,30,type = "Ljung")
accuracy(fc_wq)
fit_q$AIC

seasonal.ma1.sim <- arima.sim(list(order = c(0,12,0)), n = 200)
plot(seasonal.ma1.sim)
?arima.sim
?arima
