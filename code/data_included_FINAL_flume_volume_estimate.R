# flows in ml/s: 
flow <- c(131, 132, 132, 131, 134, 129, 130, 129, 130, 129, 130, 127, 127, 127, 126, 123, 126, 125, 124, 125, 127, 125, 123, 122, 125, 124, 118, 123, 125, 125, 120, 122, 124, 126, 122, 117, 123, 118, 120, 120, 121, 118, 118, 118, 119, 118, 120, 117, 117, 114, 115, 115, 116, 112, 112, 110, 111, 108, 105, 105, 102, 100, 100, 102)

# times in fractional hrs of day:  
time <- c(12.13333, 12.21667, 12.28333, 12.41667, 12.5, 12.56667, 12.6, 13.11667, 13.15, 13.18333, 13.58333, 13.61667, 13.65, 14.25, 14.28333, 14.31667, 14.35, 14.38333, 14.68333, 14.71667, 14.75, 14.78333, 14.81667, 14.85, 14.88333, 14.93333, 15.35, 15.38333, 15.41667, 15.45, 15.48333, 15.76667, 15.8, 15.86667, 15.9, 15.93333, 15.96667, 16.18333, 16.23333, 16.28333, 16.31667, 16.33333, 16.36667, 16.55, 16.58333, 16.61667, 16.81667, 16.85, 16.9, 16.93333, 16.98333, 17.03333, 17.06667, 17.2, 17.23333, 17.35, 17.36667, 17.53333, 17.58333, 17.61667, 17.73333, 17.75, 17.78333, 17.8)

#flow start time:
t_0 <- 11.59166667

#time at which flume was drained to top of test-array holding basin
t_b <- 16.96666667

#time at which basin was drained
t_f <- 17.825

#polynomial model fits:
fits <- lapply(1:4, function(x) {lm(flow ~ poly(time, x, raw = TRUE))})

#finite integral calculator
poly_fin_integr <- function(model,lim1,lim2) {
  sum(coefficients(model)/(1:length(coefficients(model)))*
        lim2^(1:length(coefficients(model))))-
  sum(coefficients(model)/(1:length(coefficients(model)))*
        lim1^(1:length(coefficients(model))))
}

estimates <- lapply(fits, function(x) poly_fin_integr(x, lim1 = t_0, lim2 = t_b))

#convert from mL/s*h to cubic m: 1 mL/s*h = 3600 mL = 3.6 L = 0.0036 cubic m

estimates <- unlist(estimates)*.0036

names(estimates) <- c("linear", "quadratic", "cubic", "4th order")

estimates
