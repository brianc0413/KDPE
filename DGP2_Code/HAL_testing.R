# HAL testing for longitudinal setting
source("long_helper.R")

nsims <- 350

hal_data <- data.frame(ate = rep(0, nsims), rr = rep(0, nsims), or = rep(0, nsims))


for (j in 1:nsims){
  set.seed(j)
  dat <- generate_long_data(n)
  # Base estimators that need debiasing
  Q_Y <- fit_hal(Y = dat$Y, X= dat[,1:4], family = "binomial",
                 fit_control = list(cv_select = TRUE,
                                    use_min = TRUE))
  
  Q_L1 <- fit_hal(Y = dat$L1, X = dat[, 1:2], family = "binomial",
                  fit_control = list(cv_select = TRUE,
                                     use_min = TRUE))
  
  
  # Nuisance Estimators that we do not care about
  g1 <- SuperLearner(dat$A1, 
                     dat[,1:3],
                     SL.library = SL.lib,
                     family = binomial())
  
  g0 <- SuperLearner(dat$A0,
                     data.frame(W=dat[,1]),
                     SL.library = SL.lib,
                     family = binomial())
  
  
  
  #######
  
  # Plug-in estimate of ATE
  all_zero <- dat
  all_zero$A0 <- 0
  all_zero$A1 <- 0
  
  all_zero_0 <- all_zero
  all_zero_1 <- all_zero
  all_zero_0$L1 <- 0
  all_zero_1$L1 <- 1
  
  all_zero$Y <- predict(object = Q_L1, new_data = all_zero[,1:2])  * predict(object = Q_Y, new_data = all_zero_1[,1:4]) + 
    (1-predict(object = Q_L1, new_data = all_zero[,1:2]))  * predict(object = Q_Y, new_data = all_zero_0[,1:4]) 
  
  #mean(all_zero$Y)
  all_ones <- dat
  all_ones$A0 <- 1
  all_ones$A1 <- 1
  all_ones_0 <- all_ones
  all_ones_1 <- all_ones
  all_ones_0$L1 <- 0
  all_ones_1$L1 <- 1
  
  
  
  all_ones$L1 <- predict(object = Q_L1, new_data = all_ones[,1:2]) 
  all_ones$Y <- predict(object = Q_L1, new_data = all_ones[,1:2]) * predict(object = Q_Y, new_data = all_ones_1[,1:4]) + 
    (1-predict(object = Q_L1, new_data = all_ones[,1:2]) ) * predict(object = Q_Y, new_data = all_ones_0[,1:4])
  mean(all_ones$Y)
  mu1_hal <- mean(all_ones$Y)
  mu0_hal <- mean(all_zero$Y)
  hal_data[j,] <- c(mu1_hal-mu0_hal, mu1_hal/mu0_hal, (mu1_hal/(1-mu1_hal))/(mu0_hal/(1-mu0_hal)))
  
  
}

rmse(hal_data$ate, true_mu_11-true_mu_00)
rmse(hal_data$rr, true_mu_11/true_mu_00)
rmse(hal_data$or, (true_mu_11/(1-true_mu_11))/(true_mu_00/(1-true_mu_00)))


true_mu_11 <- 0.420
true_mu_00 <- 0.527
true_ate <- true_mu_11-true_mu_00
true_var_ate <- 0.00085

x_fit <- seq(-0.4, 
             0.4, length = 100)
y_fit <- dnorm(x_fit, mean = true_ate, sd = sqrt(true_var_ate))
x_breaks <- seq(from = -1*0.4,to= 0.4,length.out= 30)

hist(hal_data$ate, prob = TRUE,
     xlab = "HAL Plug-in ATE Estimates", 
     breaks = x_breaks,
     xlim = c(-0.4, 0.4), main = "HAL Plug-in ATE Distribution",
     col = "orange", ylim = c(0, 15))
lines(x_fit, y_fit, col = "purple", lwd = 2)+ 
  legend("topright", c("HAL", "Limit"), fill=c("orange", "purple"))

hist(hal_data$ate)

