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
n = nrow(df) # number of data points

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
  xlab("ar1")

ggplot() + geom_density(aes(as.numeric(param.arr[,2]))) +
  geom_vline(xintercept = current.parameters.estimate[2], color = "blue", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = mean(param.arr[,2]), color = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("omega")

ggplot() + geom_density(aes(as.numeric(param.arr[,3]))) +
  geom_vline(xintercept = current.parameters.estimate[3], color = "blue", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = mean(param.arr[,3]), color = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("alpha1")

ggplot() + geom_density(aes(as.numeric(param.arr[,4]))) +
  geom_vline(xintercept = current.parameters.estimate[4], color = "blue", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = mean(param.arr[,4]), color = "red", linetype = "dashed", linewidth = 0.5) +
  xlab("beta1")

### 5. Task V: Predictions
N.simul = 100
n.test = 500
n.train = n - n.test
apparent.error.arr = integer(N.simul)
prediction.error.arr = integer(N.simul)
apparent.error.bootstrap.arr = integer(N.simul)
prediction.error.bootstrap.arr = integer(N.simul)

for (i in 1:N.simul){
  
  ### 5.1: Simulate according to ARMA(1, 0) + GARCH(1, 1) model
  garch.spec.temp = garchSpec(model = list(ar=c(0, current.parameters.estimate[1]), alpha=current.parameters.estimate[3], beta=current.parameters.estimate[4], omega=current.parameters.estimate[2]))
  simulated_garch = garchSim(spec = garch.spec.temp, n = n, extended=T)
  # plot(simulated_garch, main = "Simulated GARCH(1,1) Process", ylab = "Volatility")
  # ggplot() + geom_line(aes(x=1:length(simulated_garch$garch),y=simulated_garch$garch)) + geom_line(aes(x=1:length(df$logReturn),y=df$logReturn), color="red", linewidth=0.1)
  # ggplot() + geom_line(aes(x=1:length(df$volatility),y=df$volatility), color = "red")+ geom_line(aes(x=1:length(simulated_garch$sigma),y=simulated_garch$sigma))
  
  ### 5.2: Simulate according to epsilon_t ~ N(0, sigma_t^2)
  simulated_garch.prime = integer(n)
  simulated_garch.h.prime = integer(n)
  for (j in 1:n){
    if (j == 1){
      simulated_garch.h.prime[j] = current.parameters.estimate[2]
      simulated_garch.prime[j] = sqrt(current.parameters.estimate[2])*sample(standardised.residuals.iv, 1, replace=T)
    } else {
      # rt = ar*r_{t - 1} + Et (Et ~ sigma_t epsilon_t)
      simulated_garch.prime[j] = current.parameters.estimate[1]*simulated_garch.prime[j - 1] + sqrt(simulated_garch.h.prime[j - 1])*sample(standardised.residuals.iv, 1, replace=T)
      # ht = omega + beta*h_{t-1} + alpha*r_{t-1}
      simulated_garch.h.prime[j] = current.parameters.estimate[2] + current.parameters.estimate[3]*simulated_garch.prime[j - 1]^2 + current.parameters.estimate[4]*simulated_garch.h.prime[j - 1]
    }
  }
  
  ### 6. Task VI: Measure empirical prediction error
  
  ### 6.1 Method I
  train.set.r.temp = head(simulated_garch$garch, n.train)
  train.set.vol.temp = head(simulated_garch$sigma, n.train)
  test.set.r.temp = tail(simulated_garch$garch, n.test)
  test.set.vol.temp = tail(simulated_garch$sigma, n.test)
  
  garch_model.simul <- garchFit(~arma(1, 0) + garch(1, 1), data = train.set.r.temp, include.mean = FALSE)
  summary(garch_model.simul)
  coef.simul = garch_model.simul@fit$matcoef[,1]
  
  # pred.temp = predict(garch_model.simul, n.ahead = n.test) # this doesn't work - this only predicts returns
  estimated.volatility.simul <- garch_model.simul@h.t
  pred.h.simul = integer(n.test)
  pred.sigma.simul = integer(n.test)
  for (ii in 1:n.test){
    if (ii == 1){
      # warning: square the vol
      pred.h.simul[ii] = coef.simul["omega"] + train.set.r.temp[n.train]^2*coef.simul["alpha1"] + train.set.vol.temp[n.train]^2*coef.simul["beta1"]
    } else {
      # we observed test.r[ii - 1]
      pred.h.simul[ii] = coef.simul["omega"] + test.set.r.temp[ii - 1]^2*coef.simul["alpha1"] + pred.h.simul[ii - 1]*coef.simul["beta1"]
    }
  }
  
  pred.sigma.simul = sqrt(pred.h.simul)
  
  ggplot() + geom_line(aes(x=1:length(pred.sigma.simul),y=pred.sigma.simul), color = "red")+ geom_line(aes(x=1:length(test.set.vol.temp),y=test.set.vol.temp))
  ggplot() + geom_line(aes(x=1:length(train.set.vol.temp),y=train.set.vol.temp), color = "red")+ geom_line(aes(x=1:length(garch_model@sigma.t),y=garch_model@sigma.t))
  prediction.error.arr[i] = sum(pred.sigma.simul - test.set.vol.temp)^2/n.test
  apparent.error.arr[i] = sum(train.set.vol.temp - garch_model.simul@sigma.t)^2/n.train
  
  # Density plot of prediction error (looks normal)
  # ggplot() + geom_density(aes(pred.sigma.simul - test.set.vol.temp)) + xlab("col")
  
  ### 6.2 Method II
  train.set.r.temp = head(simulated_garch.prime, n.train)
  train.set.vol.temp = head(sqrt(simulated_garch.h.prime), n.train)
  test.set.r.temp = tail(simulated_garch.prime, n.test)
  test.set.vol.temp = tail(sqrt(simulated_garch.h.prime), n.test)
  
  garch_model.simul <- garchFit(~arma(1, 0) + garch(1, 1), data = train.set.r.temp, include.mean = FALSE)
  summary(garch_model.simul)
  coef.simul = garch_model.simul@fit$matcoef[,1]
  
  # pred.temp = predict(garch_model.simul, n.ahead = n.test) # this doesn't work - this only predicts returns
  estimated.volatility.simul <- garch_model.simul@h.t
  pred.h.simul = integer(n.test)
  pred.sigma.simul = integer(n.test)
  for (ii in 1:n.test){
    if (ii == 1){
      # warning: square the vol
      pred.h.simul[ii] = coef.simul["omega"] + train.set.r.temp[n.train]^2*coef.simul["alpha1"] + train.set.vol.temp[n.train]^2*coef.simul["beta1"]
    } else {
      # we observed test.r[ii - 1]
      pred.h.simul[ii] = coef.simul["omega"] + test.set.r.temp[ii - 1]^2*coef.simul["alpha1"] + pred.h.simul[ii - 1]*coef.simul["beta1"]
    }
  }
  
  pred.sigma.simul = sqrt(pred.h.simul)
  
  # ggplot() + geom_line(aes(x=1:length(pred.sigma.simul),y=pred.sigma.simul), color = "red")+ geom_line(aes(x=1:length(test.set.vol.temp),y=test.set.vol.temp))
  # ggplot() + geom_line(aes(x=1:length(train.set.vol.temp),y=train.set.vol.temp), color = "red")+ geom_line(aes(x=1:length(garch_model@sigma.t),y=garch_model@sigma.t))
  prediction.error.bootstrap.arr[i] = sum(pred.sigma.simul - test.set.vol.temp)^2/n.test
  apparent.error.bootstrap.arr[i] = sum(train.set.vol.temp - garch_model.simul@sigma.t)^2/n.train
  
}

ggplot() + geom_density(aes(apparent.error.arr), color="blue") + geom_density(aes(apparent.error.bootstrap.arr), color="red") + xlab("Training Error") + labs(title = "Distribution of Training Error (Bootstrap in Red)")
ggplot() + geom_density(aes(prediction.error.arr), color="blue") + geom_density(aes(prediction.error.bootstrap.arr), color="red") + xlab("Test Error") + labs(title = "Distribution of Test Error (Bootstrap in Red)")
# Results show that the bootstrap procedure is pretty good

### 6.3: Apply to actual data

train.set.r.actual = head(df$logReturn, n.train)
train.set.vol.actual = head(df$volatility, n.train)
test.set.r.actual = tail(df$logReturn, n.test)
test.set.vol.actual = tail(df$volatility, n.test)

garch_model.actual <- garchFit(~arma(1, 0) + garch(1, 1), data = train.set.r.actual, include.mean = FALSE)
summary(garch_model.actual)
coef.actual = garch_model.actual@fit$matcoef[,1]

# pred.temp = predict(garch_model.simul, n.ahead = n.test) # this doesn't work - this only predicts returns
estimated.volatility.actual <- garch_model.actual@h.t
pred.h.actual = integer(n.test)
pred.sigma.actual = integer(n.test)
for (ii in 1:n.test){
  if (ii == 1){
    # warning: square the vol
    pred.h.actual[ii] = coef.actual["omega"] + train.set.r.actual[n.train]^2*coef.actual["alpha1"] + train.set.vol.actual[n.train]^2*coef.actual["beta1"]
  } else {
    # we observed test.r[ii - 1]
    pred.h.actual[ii] = coef.actual["omega"] + test.set.r.actual[ii - 1]^2*coef.actual["alpha1"] + pred.h.actual[ii - 1]*coef.actual["beta1"]
  }
}

pred.sigma.actual = sqrt(pred.h.actual)

# ggplot() + geom_line(aes(x=1:length(pred.sigma.simul),y=pred.sigma.simul), color = "red")+ geom_line(aes(x=1:length(test.set.vol.temp),y=test.set.vol.temp))
# ggplot() + geom_line(aes(x=1:length(train.set.vol.temp),y=train.set.vol.temp), color = "red")+ geom_line(aes(x=1:length(garch_model@sigma.t),y=garch_model@sigma.t))
prediction.error.actual = sum(pred.sigma.actual - test.set.vol.actual)^2/n.test
apparent.error.actual = sum(train.set.vol.actual - garch_model.actual@sigma.t)^2/n.train
# prediction.error.actual = 9.855875e-06
# apparent.error.actual = 0.0001105498

combined_df_plot = data.frame(c1 = seq(1, length(pred.sigma.actual)), c2 = pred.sigma.actual, c3 = test.set.vol.actual)

# Plotting both data sets together
ggplot(combined_df_plot) +
  geom_line(aes(x = c1, y = c2, color = "Predicted Sigma Actual")) +
  geom_line(aes(x = c1, y = c3, color = "Test Set Vol Actual")) +
  labs(title = "Predicted Sigma Actual vs Test Set Vol Actual",
       x = "Index",
       y = "Value") +
  scale_color_manual(name = "Data", values = c("Predicted Sigma Actual" = "blue", "Test Set Vol Actual" = "red"))

### 7: Task VII Apply Shrinkage
ccf(pred.sigma.actual - test.set.vol.actual, test.set.vol.actual)

N.simul.lambda = 10
lambda.sweep = seq(0, 1.3, 0.05)
n.sweep = length(lambda.sweep)

shrinkage.prediction.error.arr = integer(n.sweep)
shrinkage.prediction.error.bootstrap.arr = integer(n.sweep)

counter = 0
for (lambda in lambda.sweep){
  counter = counter + 1
  for (i in 1:N.simul.lambda){
    
    ### 7.1: Simulate according to ARMA(1, 0) + GARCH(1, 1) model
    garch.spec.temp = garchSpec(model = list(ar=c(0, current.parameters.estimate[1]), alpha=current.parameters.estimate[3], beta=current.parameters.estimate[4], omega=current.parameters.estimate[2]))
    simulated_garch = garchSim(spec = garch.spec.temp, n = n, extended=T)
    
    ### 7.2: Simulate according to epsilon_t ~ N(0, sigma_t^2)
    simulated_garch.prime = integer(n)
    simulated_garch.h.prime = integer(n)
    for (j in 1:n){
      if (j == 1){
        simulated_garch.h.prime[j] = current.parameters.estimate[2]
        simulated_garch.prime[j] = sqrt(current.parameters.estimate[2])*sample(standardised.residuals.iv, 1, replace=T)
      } else {
        # rt = ar*r_{t - 1} + Et (Et ~ sigma_t epsilon_t)
        simulated_garch.prime[j] = current.parameters.estimate[1]*simulated_garch.prime[j - 1] + sqrt(simulated_garch.h.prime[j - 1])*sample(standardised.residuals.iv, 1, replace=T)
        # ht = omega + beta*h_{t-1} + alpha*r_{t-1}
        simulated_garch.h.prime[j] = current.parameters.estimate[2] + current.parameters.estimate[3]*simulated_garch.prime[j - 1]^2 + current.parameters.estimate[4]*simulated_garch.h.prime[j - 1]
      }
    }
    
    ### 7.3 Method I
    train.set.r.temp = head(simulated_garch$garch, n.train)
    train.set.vol.temp = head(simulated_garch$sigma, n.train)
    test.set.r.temp = tail(simulated_garch$garch, n.test)
    test.set.vol.temp = tail(simulated_garch$sigma, n.test)
    
    garch_model.simul <- garchFit(~arma(1, 0) + garch(1, 1), data = train.set.r.temp, include.mean = FALSE)
    summary(garch_model.simul)
    coef.simul = garch_model.simul@fit$matcoef[,1]
    
    # pred.temp = predict(garch_model.simul, n.ahead = n.test) # this doesn't work - this only predicts returns
    estimated.volatility.simul <- garch_model.simul@h.t
    pred.h.simul = integer(n.test)
    pred.sigma.simul = integer(n.test)
    for (ii in 1:n.test){
      if (ii == 1){
        # warning: square the vol
        # apply the shrinkage
        pred.h.simul[ii] = coef.simul["omega"] + train.set.r.temp[n.train]^2*coef.simul["alpha1"] + train.set.vol.temp[n.train]^2*coef.simul["beta1"]
      } else {
        # we observed test.r[ii - 1]
        # apply the shrinkage
        pred.h.simul[ii] = coef.simul["omega"] + test.set.r.temp[ii - 1]^2*coef.simul["alpha1"] + pred.h.simul[ii - 1]*coef.simul["beta1"]
      }
    }
    
    pred.sigma.simul = sqrt(pred.h.simul)
    
    ggplot() + geom_line(aes(x=1:length(pred.sigma.simul),y=pred.sigma.simul), color = "red")+ geom_line(aes(x=1:length(test.set.vol.temp),y=test.set.vol.temp))
    shrinkage.prediction.error.arr[counter] = shrinkage.prediction.error.arr[counter] + sum(pred.sigma.simul - test.set.vol.temp)^2/n.test
 
    # Density plot of prediction error (looks normal)
    # ggplot() + geom_density(aes(pred.sigma.simul - test.set.vol.temp)) + xlab("col")
    
    ### 7.4 Method II
    train.set.r.temp = head(simulated_garch.prime, n.train)
    train.set.vol.temp = head(sqrt(simulated_garch.h.prime), n.train)
    test.set.r.temp = tail(simulated_garch.prime, n.test)
    test.set.vol.temp = tail(sqrt(simulated_garch.h.prime), n.test)
    
    garch_model.simul <- garchFit(~arma(1, 0) + garch(1, 1), data = train.set.r.temp, include.mean = FALSE)
    summary(garch_model.simul)
    coef.simul = garch_model.simul@fit$matcoef[,1]
    
    # pred.temp = predict(garch_model.simul, n.ahead = n.test) # this doesn't work - this only predicts returns
    estimated.volatility.simul <- garch_model.simul@h.t
    pred.h.simul = integer(n.test)
    pred.sigma.simul = integer(n.test)
    for (ii in 1:n.test){
      if (ii == 1){
        # warning: square the vol
        pred.h.simul[ii] = coef.simul["omega"] + train.set.r.temp[n.train]^2*coef.simul["alpha1"] + train.set.vol.temp[n.train]^2*coef.simul["beta1"]
      } else {
        # we observed test.r[ii - 1]
        pred.h.simul[ii] = coef.simul["omega"] + test.set.r.temp[ii - 1]^2*coef.simul["alpha1"] + pred.h.simul[ii - 1]*coef.simul["beta1"]
      }
    }
    
    pred.sigma.simul = sqrt(pred.h.simul)
    
    shrinkage.prediction.error.bootstrap.arr[counter] = sum(pred.sigma.simul - test.set.vol.temp)^2/n.test
  }
  shrinkage.prediction.error.arr[counter] = shrinkage.prediction.error.arr[counter] / N.simul.lambda
  shrinkage.prediction.error.bootstrap.arr[counter] = shrinkage.prediction.error.bootstrap.arr[counter] / N.simul.lambda
}

combined_df_plot = data.frame(c1 = lambda.sweep, c2 = shrinkage.prediction.error.arr, c3 = shrinkage.prediction.error.bootstrap.arr)

# Plotting both data sets together
ggplot(combined_df_plot[17:23,]) +
  geom_line(aes(x = c1, y = c2, color = "Shrinkage Method I")) +
  geom_line(aes(x = c1, y = c3, color = "Shrinkage Method II")) +
  labs(title = "MSPE under Shrinkage",
       x = "lambda",
       y = "MSPE") +
  scale_color_manual(name = "Data", values = c("Shrinkage Method I" = "blue", "Shrinkage Method II" = "red"))

### 7.5 Predict on real data set
# -- OMIT -- Naive shrinkage does not improve prediction

### 7.6 Thresholding lambda
# -- OMIT -- Naive shrinkage does not improve prediction

### 8: Forecasting with Different Time Horizons (Investigate the Stability of Parameters)
rolling.window.size = 200
param.arr.window = matrix(0, nrow = (n.train - rolling.window.size), ncol = length(current.parameters.estimate))
for (i in 1:(n.train - rolling.window.size)){
  # 8.1 Fit GARCH
  garch_model.bootstrap.window <- garchFit(~arma(1, 0) + garch(1, 1), data = df$logReturn[i:(i+rolling.window.size)], include.mean = FALSE, trace=F)
  
  # 8.2 Collect Garch Parameters
  param.arr.window[i,] = garch_model.bootstrap.window@fit$matcoef[,1]
}

# 8.3 Plot of parameters
# ar1
ggplot() + geom_line(aes(x=1:length(param.arr.window[,1]),y=param.arr.window[,1])) + xlab("ar1")+
  geom_hline(yintercept = mean(param.arr[,1]), color = "red")

# omega
ggplot() + geom_line(aes(x=1:length(param.arr.window[,2]),y=param.arr.window[,2])) + xlab("omega")+
  geom_hline(yintercept = mean(param.arr[,2]), color = "red")

# alpha1
ggplot() + geom_line(aes(x=1:length(param.arr.window[,3]),y=param.arr.window[,3])) + xlab("alpha1")+
  geom_hline(yintercept = mean(param.arr[,3]), color = "red")

# beta1
ggplot() + geom_line(aes(x=1:length(param.arr.window[,4]),y=param.arr.window[,4])) + xlab("beta1")+
  geom_hline(yintercept = mean(param.arr[,4]), color = "red")

# See a decreasing trend in beta1
# See an increasing trend in alpha1
# Magnitude of omega stays more or less 0
# AR1 parameter having some fluctuations

combined.alpha1 = data.frame(matrix(nrow = length(param.arr.window[,3]), ncol = 0))
combined.alpha1$y = param.arr.window[,3]
combined.alpha1$x = 1:length(param.arr.window[,3])
lm.alpha1 = lm(y~., data=combined.alpha1)
summary(lm.alpha1)

ggplot() + geom_line(aes(x=1:length(param.arr.window[,3]),y=param.arr.window[,3])) + xlab("alpha1")+
  geom_line(aes(x=1:length(param.arr.window[,3]),y=lm.alpha1$fitted.values), color="blue")+
  geom_hline(yintercept = mean(param.arr[,3]), color = "red")

combined.beta1 = data.frame(matrix(nrow = length(param.arr.window[,4]), ncol = 0))
combined.beta1$y = param.arr.window[,4]
combined.beta1$x = 1:length(param.arr.window[,4])
lm.beta1 = lm(y~., data=combined.beta1)
summary(lm.beta1)

ggplot() + geom_line(aes(x=1:length(param.arr.window[,4]),y=param.arr.window[,4])) + xlab("beta1")+
  geom_line(aes(x=1:length(param.arr.window[,4]),y=lm.beta1$fitted.values), color="blue")+
  geom_hline(yintercept = mean(param.arr[,4]), color = "red")

# 8.4 Idea is to shrink towards the regression line










