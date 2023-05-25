library(tseries)
library(forecast)
library(lmtest)
library(TSA)
library(car)
library(xts)
library(zoo)
library(imputeTS)
library(DataCombine)
library(Metrics)

setwd("C:/Users/devan/Documents/School/Time Series/Final")

# Import data
cost_data <- read.csv('costs.csv', header=TRUE)
colnames(cost_data)[1] <- 'Date'
fuel_data <- read.csv('fuel.csv', header=TRUE)
colnames(fuel_data)[1] <- 'Date'
ltratio_data <- read.csv('ltratio.csv', header=TRUE)
colnames(ltratio_data)[1] <- 'Date'

# Convert dates to date format
cost_data$Date <- as.Date(cost_data$Date, '%Y-%m-%d')
fuel_data$Date <- as.Date(fuel_data$Date, '%Y-%m-%d')
ltratio_data$Date <- as.Date(ltratio_data$Date, '%Y-%m-%d')

# Check for missing values
sum(is.na(cost_data$LH.CPM))
sum(is.na(fuel_data$DOE.Avg))
sum(is.na(ltratio_data$LT.Ratio))

# Convert each series into an xts
cost_xts <- xts(cost_data$LH.CPM, order.by = cost_data$Date)
attr(cost_xts, 'frequency') <- 52
colnames(cost_xts)[1] <- 'LH.CPM'
fuel_xts <- xts(fuel_data$DOE.Avg, order.by = fuel_data$Date)
attr(fuel_xts, 'frequency') <- 52
colnames(fuel_xts)[1] <- 'DOE.Avg'
lt_xts <- xts(ltratio_data$LT.Ratio, order.by = ltratio_data$Date)
attr(lt_xts, 'frequency') <- 52
colnames(lt_xts)[1] <- 'LT.Ratio'

# View each series
plot(cost_xts, main = "LH CPMs 2019 to April 2023")
plot(fuel_xts, main = "Fuel 2019 to April 2023")
plot(lt_xts, main = "LT Ratio 2019 to April 2023")

hist(cost_xts$LH.CPM, main = "LH CPM Distribution", xlab = "LH CPM")

# Interpolate gaps in cost and load to truck ratios, review results
cost_interp <- na_seadec(cost_xts)
lt_interp <- na_seadec(lt_xts)

plot(cost_interp, ylab='LH CPM', main = 'LH CPM with Interpolation', col='red')
lines(cost_xts, lwd=2, col='black')
legend(x="top", legend = c("Actual", "Imputed"), col = c("black", "red"), lty = 1, xpd = TRUE)

plot(lt_interp, ylab='LT Ratio', main = 'Load to Truck Ratio with Interpolation', col='red')
lines(lt_xts, lwd=2, col='black')
legend(x="top", legend = c("Actual", "Imputed"), col = c("black", "red"), lty = 1, xpd = TRUE)

# Review stationarity of variables

adf.test(cost_interp$LH.CPM)
kpss.test(cost_interp$LH.CPM, null="Trend")
acf(cost_interp$LH.CPM, main='LH CPM ACF')

adf.test(lt_interp$LT.Ratio)
kpss.test(lt_interp$LT.Ratio, null="Trend")

adf.test(fuel_data$DOE.Avg)
kpss.test(fuel_data$DOE.Avg, null="Trend")

# Find lag of LT with most explanatory power and view correlation between fuel and cost 

lt_ccf_result <- ccf(as.ts(cost_interp$LH.CPM), as.ts(lt_interp$LT.Ratio))
plot(lt_ccf_result, main = "L:T and LH CPM Cross-Correlation Function")
lt_significant_lags <- lt_ccf_result$lag[which.max(abs(lt_ccf_result$acf))]
lt_significant_lags*52

cor(cost_interp$LH.CPM,fuel_xts$DOE.Avg)

# Construct dummy variable for COVID impact on supply chains

min_index <- which.min(cost_interp$LH.CPM)
cost_interp[min_index,]

cost_interp_cut <- window(cost_interp, start = as.Date("2022-07-02"), end = as.Date("2023-01-01"))
min_index_cut <- which.min(cost_interp_cut$LH.CPM)
cost_interp_cut[min_index_cut,]

# will use 2020-05-07 through 2022-11-19 as bounds for COVID impact based on cost time points
# and shift of series values near these boundaries

cost_data$CovDum <- ifelse(cost_data$Date >= as.Date("2020-05-07") & cost_data$Date <= as.Date("2022-11-19"), 1, 0)
covdum_xts <- xts(cost_data$CovDum, order.by = cost_data$Date)
colnames(covdum_xts)[1] <- "cov19_dummy"

# Create slid L:T variable, fill in N/As at front of series, re-impute missing values

ltslide <- slide(ltratio_data, Var = "LT.Ratio", TimeVar="Date", NewVar='LT.Ratiol2', slideBy = 2)
ltslide$LT.Ratiol2[224] <- 1.91
ltslide$LT.Ratiol2[225] <- 1.99
ltslide$LT.Ratiol2[226] <- 3.00

ltl2_xts <- xts(ltslide$LT.Ratiol2, order.by = ltslide$Date)
attr(ltl2_xts, 'frequency') <- 52
colnames(ltl2_xts)[1] <- 'LT.Ratio'
ltl2_interp <- na_seadec(ltl2_xts)

# Create fully combined data set then test train split

data_comb <- merge(cost_interp, ltl2_interp, fuel_xts, covdum_xts)
cost_train <- window(data_comb, end = as.Date("2023-01-31"))
cost_test <- window(data_comb, start = as.Date("2023-01-31"))

# Fit SARIMA Model & Review Residuals

sarima <- auto.arima(cost_train$LH.CPM, seasonal = TRUE, max.order = 6)
summary(sarima)
plot(forecast(sarima, h=13))

acf(sarima$residuals, main = 'SARIMA Residuals ACF')
pacf(sarima$residuals)
mean(sarima$residuals)
Box.test(sarima$residuals, lag=10, type = "Ljung-Box")
shapiro.test(sarima$residuals)
bptest(sarima$residuals ~ fitted(sarima))


# Fit Dynamic Regression Model & Review Residuals

ex_regs_train = as.matrix(cost_train[,c('LT.Ratio','DOE.Avg','cov19_dummy')])
dyn_reg <- auto.arima(cost_train$LH.CPM, ic = 'aicc', xreg = ex_regs_train, max.order = 6)
summary(dyn_reg)
ex_regs_test = as.matrix(cost_test[,c('LT.Ratio','DOE.Avg','cov19_dummy')])
plot(forecast(dyn_reg, h = 13, xreg = ex_regs_test))

acf(dyn_reg$residuals, main = 'Dynamic Regression Residuals ACF')
pacf(dyn_reg$residuals)
mean(dyn_reg$residuals)
Box.test(dyn_reg$residuals, lag=10, type = "Ljung-Box")
shapiro.test(dyn_reg$residuals)
bptest(dyn_reg$residuals ~ ex_regs_train)

# Compare sMAPE accuracy of two models & Plot forecasts over true data

sarima_preds <- forecast(sarima, h=13)$mean
smape(as.numeric(cost_test$LH.CPM), sarima_preds)

dynreg_preds <- forecast(dyn_reg, h=13, xreg = ex_regs_test)$mean
smape(as.numeric(cost_test$LH.CPM), dynreg_preds)

cost_finals <- cost_test$LH.CPM
cost_finals$sarima.cpm <- sarima_preds
cost_finals$dynreg.cpm <- dynreg_preds

plot(cost_interp, main = "LH CPM with Forecasts")
lines(cost_finals$sarima.cpm, lwd=2, col='blue')
lines(cost_finals$dynreg.cpm, lwd=2, col='red')
legend(x="top", legend = c("Actual", "SARIMA", "Dyn Reg"), col = c("black", "blue", "red"), lty = 1, xpd = TRUE)
