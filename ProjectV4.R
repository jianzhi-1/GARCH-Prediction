library(TTR)
library(xts)
library(astsa)
library(quantmod)
library(ggplot2) # plotting
library(naniar) # fill na
library(visdat)
library(dplyr)
library(ggmosaic) # Mosaic plot
library(gridExtra) # multipanel plots
library(MASS) # mvrnorm
library(mvtnorm)
library(tseries)
library(fGarch)
library(zoo)

### 2. Task II: Data extraction and processing and visualisation
getSymbols('SPY')
df.spy = as.data.frame(SPY)
logReturn = diff(df.spy$SPY.Close, 1)/df.spy$SPY.Close[-length(df.spy$SPY.Close)]
df = data.frame(logReturn)
df$logReturn2 = df$logReturn^2
df$date = rownames(df.spy)[2:length(rownames(df.spy))]
# write.csv(df, file="data.csv", row.names=F)

# Calculates rolling volatility
calculate.rolling.volatility <- function(returns, window_size) {
  volatility <- rollapply(returns^2, window_size, mean, align = "right", fill = NA, partial = TRUE)
  volatility <- sqrt(volatility)
  return(volatility)
}

window.size = 31
df$volatility = calculate.rolling.volatility(df$logReturn, window.size) # we actually don't observe this
df$volatility = c(rep(NA, (window.size - 1)/2), df$volatility[window.size:nrow(df)], rep(NA, (window.size - 1)/2)) # we actually don't observe this
df = na.omit(df)

# Density plot of log returns
ggplot(df, aes(x = logReturn)) +
  geom_density(fill = "skyblue", color = "black") +
  labs(title = "Density Plot of Numbers", x = "Numbers", y = "Density") +
  theme_minimal()

# ACF and PACF Plots
acf.r = acf(df$logReturn)
acf.r2 = acf(df$logReturn2)
pacf.r = pacf(df$logReturn)
pacf.r2 = pacf(df$logReturn2)

# Try the ARIMA models suggested from ACF and PACF plots
sarima(df$logReturn, 1, 0, 0) # AIC = -5.909567  AICc = -5.909567  BIC = -5.905137 
sarima(df$logReturn, 0, 0, 1) # AIC -5.910135  AICc = -5.910135  BIC = -5.905705 
sarima(df$logReturn, 1, 0, 1) # AIC = -5.909991  AICc = -5.90999  BIC = -5.904084 
# Conclusion is to use AR(1) for now. Doesn't seem to matter

### 3. Task III Fitting a ARMA(1, 0) + GARCH(1, 1) Model

# Warning: do not include the intercept
garch_model <- garchFit(~arma(1, 0) + garch(1, 1), data = df$logReturn, include.mean = FALSE)
summary(garch_model)
# AIC       BIC       SIC      HQIC 
# -6.454745 -6.448872 -6.454746 -6.452672 

# unstandardised residual plot
unstandardised.residuals <- data.frame(x = 1:length(garch_model@residuals), y = garch_model@residuals)

# Plot of unstandardised residuals
ggplot(unstandardised.residuals, aes(x = x, y = y)) +
  geom_line() +  # Add a line plot
  labs(x = "Index", y = "Value") +  # Add axis labels
  ggtitle("Plot of Unstandardised Residuals")  # Add title

# fitted returns plot
# WHY DO I NEED A 20* to make the plots look similar size? Did I miss something?
fitted.return.df <- data.frame(x = 1:length(garch_model@fitted), y = 20*garch_model@fitted, group = "Line 1")
actual.return.df <- data.frame(x = 1:length(garch_model@fitted), y = df$logReturn, group = "Line 2")
combined_data <- rbind(fitted.return.df, actual.return.df)

ggplot(combined_data, aes(x = x, y = y, color =group)) +
  geom_line(linewidth = 0.1) +
  scale_color_manual(values = c("blue", "red")) +
  labs(x = "X Values", y = "Y Values", title = "Two Line Graphs") +
  theme_minimal()

# sigma.t plot (estimated standard deviation)
estimated.std <- data.frame(x = 1:length(garch_model@sigma.t), y = garch_model@sigma.t)
ggplot(estimated.std, aes(x = x, y = y)) +
  geom_line() +  # Add a line plot
  labs(x = "Index", y = "Value") +  # Add axis labels
  ggtitle("Plot of Estimated Standard Deviation")  # Add title

# h.t plot (estimated variance) (Note: this is just sigma.t^2)
estimated.var <- data.frame(x = 1:length(garch_model@h.t), y = garch_model@h.t)
ggplot(estimated.var, aes(x = x, y = y)) +
  geom_line() +  # Add a line plot
  labs(x = "Index", y = "Value") +  # Add axis labels
  ggtitle("Plot of Estimated Variance")  # Add title

estimated.volatility.df <- data.frame(x = 1:length(garch_model@sigma.t), y = garch_model@sigma.t, group="Line 1")
actual.volatility.df <- data.frame(x = 1:length(garch_model@sigma.t), y = df$volatility, group = "Line 2")
combined_data <- rbind(estimated.volatility.df, actual.volatility.df)
ggplot(combined_data, aes(x = x, y = y, color =group)) +
  geom_line(linewidth = 0.1) +
  scale_color_manual(values = c("blue", "red")) +
  labs(x = "X Values", y = "Y Values", title = "Volatility") +
  theme_minimal()

# Notes on attributes and methods
# garch_model@fit (don't care)
# garch_model@model (don't care)
# garch_model@residuals (unstandardised residuals)
# garch_model@fitted (fitted returns)
# garch_model@sigma.t (estimated standard deviation)
# garch_model@h.t (estimated variance)


### 4. Task IV: Bootstrapping the Residuals
residuals.iv = garch_model@residuals
sigma.iv = garch_model@sigma.t
standardised.residuals.iv = residuals.iv/sigma.iv

# Test for normality
shapiro.test(standardised.residuals.iv) # p < 0.0001, not normal
qqline(standardised.residuals.iv) # distribution is leptokurtic i.e. kurtosis > normal

current.coef.matrix = garch_model@fit$matcoef
current.parameters.estimate = current.coef.matrix[,1]
current.parameters.stderror = current.coef.matrix[,2]

ntrials = 100 # number of bootstrap runs
param.arr = matrix(0, nrow = ntrials, ncol = length(current.parameters.estimate))
for (i in 1:ntrials){
  # 4.1 Build out the bootstrap
  bootstrap.return = sigma.iv*sample(standardised.residuals.iv, length(sigma.iv), replace=T)
 
  # 4.2 Fit GARCH
  garch_model.bootstrap <- garchFit(~arma(1, 0) + garch(1, 1), data = bootstrap.return, include.mean = FALSE, trace=F)
  
  # 4.3 Collect Garch Parameters
  param.arr[i,] = garch_model.bootstrap@fit$matcoef[,1]
}

# plot a distribution of the parameters
ggplot() + geom_density(aes(as.numeric(param.arr[,1]))) +
  geom_vline(xintercept = current.parameters.estimate[1], color = "blue", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = mean(param.arr[,1]), color = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("col")

ggplot() + geom_density(aes(as.numeric(param.arr[,2]))) +
  geom_vline(xintercept = current.parameters.estimate[2], color = "blue", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = mean(param.arr[,2]), color = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("col")

ggplot() + geom_density(aes(as.numeric(param.arr[,3]))) +
  geom_vline(xintercept = current.parameters.estimate[3], color = "blue", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = mean(param.arr[,3]), color = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("col")

ggplot() + geom_density(aes(as.numeric(param.arr[,4]))) +
  geom_vline(xintercept = current.parameters.estimate[4], color = "blue", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = mean(param.arr[,4]), color = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("col")

# 5. Task IV: Predictions
garch_predictions <- predict(garch_model, n.ahead = 10)

