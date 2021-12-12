library(tseries)
library(forecast)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(GGally)

### import and format data
csv2020<- read.csv("Data/1090006_36.01_-78.94_tmy-2020.csv",stringsAsFactors = F)#,header= FALSE)
basicinfo <- csv2020[1,]
data<-csv2020[c(-1,-2),seq(1,14,1)]
data<-mutate_all(data, function(x) as.numeric(as.character(x)))
colnames(data)<- csv2020[2,seq(1,14,1)]

data$date<-paste(paste("2020",data$Month,data$Day,sep="-"), paste(data$Hour,data$Minute,sep=":"))
data$time <-strptime(data$date,format="%Y-%m-%d  %H:%M")
data[,seq(1,14,1)]<- sapply(data[,seq(1,14,1)],as.numeric)
tsGHI<-ts(data$GHI)#,frequency=24*365)
tsTemp<-ts(data$Temperature,frequency=24*365)

tsdata <- lapply(data[,seq(6,14,1)],function(t) ts(t,frequency=24*365))

GHIhead <- window(tsGHI,start=c(1,1),end=c(1,8040))
GHItail <- window(tsGHI,start=c(1,8040),end=c(1,8760))

##############################################################################################################
# EDA


### plot ts
GHIdf <- data.frame(GHI = tsGHI, hour = seq(1,length(tsGHI),1))
ggplot(data=GHIdf[seq(12,8760,24),],aes(x=hour,y=GHI))+geom_point()+labs(x="GHI at noon over a year")

ts.plot(tsGHI,col="red3",gpars=list(xlab="Hour since Jan.1st 00:30", ylab="Hourly GHI", lty=c(1:3)))
ts.plot(tsGHI[1:2000],col="red3",gpars=list(xlab="hour", ylab="GHI", lty=c(1:3)))
ts.plot(tsGHI[1:120],col="red3",gpars=list(xlab="hour", ylab="GHI", lty=c(1:3)))
ts.plot(tsGHI[1200:1300],col="red3",gpars=list(xlab="Hour since Jan.1st 00:30", ylab="Hourly GHI", lty=c(1:3)))


### EDA for multiple linear regression
qplot(data$DHI,data$GHI)+ylab("GHI (W/m^2)")+xlab("DHI")
qplot(data$DNI,data$GHI)+ylab("GHI (W/m^2)")+xlab("DNI")
qplot(data$Temperature,data$GHI)+ylab("GHI (W/m^2)")+xlab("Temperature (Celsius)")
qplot(data$`Dew Point`,data$GHI)+ylab("GHI (W/m^2)")+xlab("Dew Point")
qplot(data$`Wind Direction`,data$GHI)+ylab("GHI (W/m^2)")+xlab("Wind Direction")
qplot(data$`Wind Speed`,data$GHI)+ylab("GHI (W/m^2)")+xlab("Wind Speed")
qplot(data$Pressure,data$GHI)+ylab("GHI (W/m^2)")+xlab("Pressure")
qplot(data$`Surface Albedo`, data$GHI) +ylab("GHI (W/m^2)")+xlab("Surface Albedo")#+ geom_point(alpha = .5,colour="blue4") #+  geom_smooth(method = (y~x))


### make solar albdeo categorical
hist(data$`Surface Albedo`)
hist(data$`Surface Albedo`[data$`Surface Albedo`<0.8])
length(data$`Surface Albedo`[data$`Surface Albedo`>0.8])


### fit linear regression model and select predictor
data$SurfaceAlbedoLevel[data$`Surface Albedo`==0.11]<-"0.11"
data$SurfaceAlbedoLevel[data$`Surface Albedo`==0.12]<-"0.12"
data$SurfaceAlbedoLevel[data$`Surface Albedo`==0.13]<-"0.13"
data$SurfaceAlbedoLevel[data$`Surface Albedo`==0.14]<-"0.14"
data$SurfaceAlbedoLevel[data$`Surface Albedo`==0.15]<-"0.15"
data$SurfaceAlbedoLevel[data$`Surface Albedo`==0.87]<-"0.87"
data$SurfaceAlbedoLevel <- as.factor(data$SurfaceAlbedoLevel)
data$hour_since <- seq(1,8760,1)

linearmodel1 <- lm(GHI~Temperature+SurfaceAlbedoLevel+Pressure+`Dew Point`+`Wind Speed`+`Wind Direction`,data=data)
summary(linearmodel1)
plot(linearmodel1)
ggplot(data, aes(x=as.numeric(hour_since), y=linearmodel1$residual)) +
  geom_point(alpha = .5,colour="blue4") +
  geom_hline(yintercept=0,col="red3") + labs(title="Residuals vs Time")

linearmodel2 <- lm(GHI~Temperature+SurfaceAlbedoLevel+Pressure+`Dew Point`+`Wind Speed`+`Wind Direction`+Hour+Month,data=data)
summary(linearmodel2)
plot(linearmodel2)

ggplot(data, aes(x=as.numeric(time), y=linearmodel2$residual)) +
  geom_point(alpha = .5,colour="blue4") +
  geom_hline(yintercept=0,col="red3") + labs(title="Residuals vs Time")


### ts stationary test
tsGHI %>% ggtsdisplay(lag.max=80,main="Hourly GHI")
adf_test <- adf.test(tsGHI,alternative = 'stationary')
print(adf_test)
kpss_test <- kpss.test(tsGHI)
print(kpss_test)
acf(tsGHI[1:8760],lag.max=80)
pacf(tsGHI[1:8760],lag.max=80)


### first-order differencing
tsDiffOne <- diff(tsGHI)
autoplot(tsDiffOne)
tsDiffOne%>%ggtsdisplay()
adf_test <- adf.test(tsDiffOne,alternative = 'stationary')
print(adf_test)
kpss_test <- kpss.test(tsDiffOne)
print(kpss_test)
acf(tsDiffOne,lag.max=200)
pacf(tsDiffOne,lag.max=100)


### Seasonal differencing using lag= 24 on top of 1-st order differencing
diff(tsGHI,24)%>%diff() ->   tsDiff24
tsDiff24%>%ggtsdisplay()
adf_test <- adf.test(tsDiff24,alternative = 'stationary')
print(adf_test)
kpss_test <- kpss.test(tsDiff24)
print(kpss_test)
acf(tsDiff24,lag.max=80,main="Seasonal differenceing with lag=24")
pacf(tsDiff24,lag.max=80,main="Seasonal differenceing with lag=24")


################################################################################################
# Model fitting

### Fit AR (26,0), takes a long time to run
# armodel <- arima(tsGHI,order=c(26,1,0))#method"CSS"cross validation
# armodel
# checkresiduals(armodel)
# acf(armodel$residuals,lag.max=80)
# pacf(armodel$residuals,lag.max=80)

### arima model on the training data
autoarimaModel1<-auto.arima(GHIhead)
autoarimaModel1
arimamodel_refit <- arima(GHIhead,order=c(4,1,2))
arimamodel_refit
arimafcast <- forecast(arimamodel_refit,h=2190)
autoplot(GHIhead) +
  autolayer(arimafcast, series="forecast", PI=FALSE)
arimafcast %>% autoplot()
accuracy(arimafcast,GHItail)

### arima (26,1,2) takes too long to run
# arima27_refit <- arima(GHIhead,order=c(26,1,2))
# arima27_refit
# checkresiduals(arima27_refit)
# acf(arima27_refit $residuals,lag.max=80)
# pacf(arima27_refit  $residuals,lag.max=80)
# autoplot(GHIhead) +
#   autolayer(ar_arimafcast, series="forecast", PI=FALSE)
# ar_arimafcast %>% autoplot()
# accuracy(ar_arimafcast,GHItail)

### arima with predictors
xreg <- cbind(Temp=data$Temperature,DewPoint=data$`Dew Point`,SurfaceAlbedoLevel=data$SurfaceAlbedoLevel,WindDir =data$`Wind Direction`,WindSpeed=data$`Wind Speed`)
model1_predictor <- auto.arima(data$GHI,xreg = xreg)
model1_predictor
checkresiduals(model1_predictor)
acf(residuals(model1_predictor))
pacf(residuals(model1_predictor))

datahead <- data[1:8040,]
datatail <- data[8040:8760,]
xreg <- cbind(Temp=datahead$Temperature,DewPoint=datahead$`Dew Point`,SurfaceAlbedoLevel=datahead$SurfaceAlbedoLevel,WindDir =datahead$`Wind Direction`,WindSpeed=datahead$`Wind Speed`)

refit_predictor <- Arima(GHIhead, order=c(4,1,2),xreg = xreg)
checkresiduals(refit_predictor,lag.max=100)
arimapredictorfcast <- forecast(refit_predictor,
                                xreg = cbind(Temp=datatail$Temperature,DewPoint=datatail$`Dew Point`,SurfaceAlbedoLevel=datatail$SurfaceAlbedoLevel,WindDir=datatail$`Wind Direction`,WindSpeed=datatail$`Wind Speed`
                                ))
autoplot(GHIhead) +
  autolayer(arimapredictorfcast, series="forecast", PI=FALSE)

arimapredictorfcast %>% autoplot(ylab="Hourly GHI")
accuracy(arimapredictorfcast,GHItail)

### Seasonal Arima based on EDA, can take a long time to run
seasonaldiffARIMA1 <- Arima(tsGHI,order=c(4,1,3),seasonal=list(order=c(1,1,1),period=24))
seasonaldiffARIMA1
checkresiduals(seasonaldiffARIMA1,lag.max=100)
acf(seasonaldiffARIMA1$residuals,lag.max=80)
pacf(seasonaldiffARIMA1$residuals,lag.max=80)
# compare using mor lags in AR
seasonaldiffARIMA2 <- Arima(tsGHI,order=c(7,1,3),seasonal=list(order=c(1,1,1),period=24))
seasonaldiffARIMA2
checkresiduals(seasonaldiffARIMA2,lag=80)
acf(seasonaldiffARIMA2$residuals,lag.max=80)
pacf(seasonaldiffARIMA2$residuals,lag.max=80)
# forecasting
seasonalarimamodel_refit <- Arima(GHIhead,order=c(4,1,3),seasonal=list(order=c(1,1,1),period=24))
seasonalarimafcast <- forecast(seasonalarimamodel_refit,h=30*24)
autoplot(GHIhead) +
  autolayer(seasonalarimafcast, series="forecast", PI=FALSE)
seasonalarimafcast %>% autoplot()
accuracy(seasonalarimafcast,GHItail)


### STL decomposition
mstsGHI<-msts(data$GHI,seasonal.periods = c(24,8760/12))#first is the daily seasonality, then the monthly: period1 =24 hours, period2 = 24*30 =8760/12
mstlGHI <- mstl(mstsGHI)
mstlGHI %>% autoplot()
mstsGHI_train <- window(mstsGHI,start=c(1,1),end=c(1,8040))
mstsGHI_test <- window(mstsGHI,start=c(1,8040),end=c(1,8760))

adjusted<- mstlGHI_train%>%seasadj()
adjustedDiff1 <- diff(adjusted)
acf(adjustedDiff1,lag.max=80)
pacf(adjustedDiff1,lag.max=80)
# arima model on the adjusted component
adjusted_model <- auto.arima(adjusted)
checkresiduals(adjustedModel,lag.max=150)

### final model
seasonalSTL0<-stlm(mstsGHI_train,modelfunction=Arima,order=c(5,1,2))#,seasonal=list(order=c(2,0,0),period=24))
seasonalSTLmodel0<-seasonalSTL0$model
seasonalSTLmodel0
checkresiduals(seasonalSTLmodel0,lag.max=150)
acf(residuals(seasonalSTLmodel0),lag.max=100)
AIC(seasonalSTL0$model)
ms_fcast <-forecast(seasonalSTL0,h=720)
ms_fcast%>%autoplot()
accuracy(ms_fcast,mstsGHI_test)

### weekahead forecast
weekahead_train <- window(mstsGHI,start=c(1,1),end=c(1,8592))
weekahead_test <- window(mstsGHI,start=c(1,8592),end=c(1,8760))
weekaheadGHI_refit <- stlm(weekahead_train,modelfunction=Arima,order=c(5,1,2))
# mstsGHI_refit<-mstsGHI_refit["model"]$model
weekahead_fcast <-forecast(weekaheadGHI_refit,h=24*7)
autoplot(weekahead_train) +
  autolayer(weekahead_fcast, series="forecast", PI=FALSE)
weekahead_fcast %>% autoplot()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
#+scale_x_date(limits =  as.Date(c("01-01", "12-31"),"%m-%d"))#,minor_breaks = as.Date(c("01-01", "04-01",\"08-01", "12-01"),"%m-%d"))

accuracy(weekahead_fcast,weekahead_test)

weekahead <- data.frame(hour=seq(1,24*7,1),ghi=weekahead_fcast$mean,ghi_actual = tsGHI[8593:8760])
weekahead$ghi[weekahead$ghi<0] <- 0
ggplot(weekahead,aes(x=hour,y=ghi,color="predicted"))+geom_line(lwd=1.5)+geom_line(aes(x=hour,y=ghi_actual,color="actual"),lwd=1)+labs(title="One-week-ahead forecast for Dec.24th-31st",y="Hourly GHI")#+geom_smooth()
