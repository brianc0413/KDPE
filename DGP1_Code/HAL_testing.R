# kdpe main simulation (Section 4.1)
library(SuperLearner)
library(tmle)
library(CVXR)
library(doParallel)
library(foreach)
library(hal9001)
source("helper.R")

# Initialize Settings 
n <- 300
n_sims <- 500

hal_ate <- rep(0,n_sims)
hal_rr <- rep(0, n_sims)
hal_or <- rep(0, n_sims)


SL.lib <- c("SL.randomForest", "SL.glm", "SL.mean")

for (i in 1:n_sims){
  set.seed(i)
  library(SuperLearner)
  library(tmle)
  library(CVXR)
  
  dat <- generate_data(n)
  Y <- dat$Y
  A <- dat$A
  X <- dat$X
  
  #Q <- SuperLearner(Y,
  #                  data.frame(A,X), 
  #                  SL.library = SL.lib, 
  #                  family = binomial())
  
  Q <- fit_hal(X = data.frame(A,X), Y = Y, family = "binomial",
               fit_control = list(cv_select = TRUE,
                                  use_min = FALSE))
  
  g <- SuperLearner(A,
                    data.frame(X),
                    SL.library = SL.lib,
                    family = binomial())
  
  
  Q_mat <- cbind(Q0W =predict(object = Q, new_data = data.frame(A=0,X)),
                 Q1W =  predict(object = Q, new_data = data.frame(A=1,X)))
  g_est <- predict(g, newdata = data.frame(X))$pred
  
  # sl fix
  hal_ate[i] <- mean(Q_mat[,2]-Q_mat[,1])
  hal_rr[i] <- mean(Q_mat[,2])/mean(Q_mat[,1])
  mu1sl <- mean(Q_mat[,2])
  mu0sl <- mean(Q_mat[,1])
  hal_or[i] <- (mu1sl/(1-mu1sl))/(mu0sl/(1-mu0sl))
  
  print(hal_ate)
  
}



# true values of mu1, mu0 based on simulations
X <- runif(1E7)
A_0 <- 0
A_1 <- 1
Y_0 <- rbinom(1E7,
              size = 1,
              prob = 0.4 + A_0*(X-0.3)^2 + 1/4*sin(40*X/pi))
Y_1 <- rbinom(1E7,
              size = 1,
              prob = 0.4 + A_1*(X-0.3)^2 + 1/4*sin(40*X/pi))
mu1 <- mean(Y_1)
mu0 <- mean(Y_0)

# what should the actual ate be? 0.123333333
ate_true <- mu1 - mu0
rr_true <- mu1/mu0
or_true <- (mu1 / (1-mu1))/(mu0 / (1-mu0))

rmse <- function(values, truth){
  return( sqrt( mean( (values - truth)^2 ) ) )
}



print(rmse(hal_ate, ate_true))
print(rmse(hal_rr, rr_true))
print(rmse(hal_or, or_true))


x_fit <- seq(-0.1, 
             0.4, length = 100)
y_fit <- dnorm(x_fit, mean = ate_true, sd = sqrt(0.002515497))
x_breaks <- seq(from = -1*0.1,to= 0.4,length.out= 25)
hist(hal_ate, prob = TRUE,
     xlab = "HAL Plug-in ATE Estimates", 
     breaks = x_breaks,
     xlim = c(-0.1, 0.4), main = "HAL Plug-in ATE Distribution",
     col = "orange", ylim = c(0, 8))
lines(x_fit, y_fit, col = "purple", lwd = 2)+ 
  legend("topright", c("HAL", "Limit"), fill=c("orange", "purple"))
