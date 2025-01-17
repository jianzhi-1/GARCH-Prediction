# Applications of Bootstrap and Shrinkage in the GARCH Model
### UC Berkeley STAT 215B Final Project

## Methodology
- Data: SPY index daily log return (2008 - 2024)
- Prediction: 1 day and 30 day volatility

1. First, an AR(1) model is fitted and partialled out from the log return. Then, a window approach is used to ensure that prediction (1 day and 30-day) at each step is made by a model trained on the same quantity of samples. The GARCH model is fitted to the window and used for prediction. Point estimates of the model parameters are obtained.
2. Repeat the above procedure via bootstrapping the innovation terms. Also repeat the procedure by sampling the standardised innovation terms from a normal distribution. Construct prediction interval using these sets of model parameter estimates.
3. Compare coverage ratio.

## Conclusion
The bootstrap method provides prediction intervals that gave a higher coverage ratio compared to resampling the innovation terms from a normal distribution. The technique of bootstrap can also test for the misspecification of the GARCH model.

## Report
See [report](https://github.com/jianzhi-1/GARCH-Prediction/blob/main/STAT215BReport.pdf) for details and citations.