# KDPE for longitudinal data
library(ltmle)
library(rje)
library(dplyr)
library(CVXR)
library(doParallel)
# Roadmap: iteratively regress, using 
#' Y = f(W,A0, L1, A1)
#' L1 = f(W,A0)

# generate data
ncores <- detectCores()-1
source("long_helper.R")
set.seed(100)
kdpe_data <- data.frame(ATE=NA, RR=NA, OR=NA, iters = NA)
#kdpe_reg_data <- data.frame(ATE=NA, RR=NA, OR=NA)

my.cluster <- parallel::makeCluster(ncores, type = "PSOCK")
print(my.cluster)
registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
# Base Learners
library(SuperLearner)
SL.lib <- c("SL.randomForest", "SL.glm", "SL.mean")
n <- 300
lambda <- c(rep(15, 100),
            rep(20, 100),
            rep(25, 100),
            rep(30, 100))

kdpe_data_lambda_tuning <- foreach(
  i = 1:350,
  .combine = rbind) %dopar% {
    source("long_helper.R")
    library(SuperLearner)
    #library(tmle)
    library(CVXR)
    library(rje)
    library(dplyr)
    set.seed(i)
    
    dat <- generate_long_data(n)
    # Base estimators that need debiasing
    Q_Y <- SuperLearner(dat$Y,
                        dat[,1:4], 
                        SL.library = SL.lib, 
                        family = binomial())
    
    Q_L1 <- SuperLearner(dat$L1, 
                         dat[, 1:2], 
                         SL.library = SL.lib,
                         family = binomial())
    
    
    # Nuisance Estimators that we do not care about
    g1 <- SuperLearner(dat$A1, 
                       dat[,1:3],
                       SL.library = SL.lib,
                       family = binomial())
    
    g0 <- SuperLearner(dat$A0,
                       data.frame(W=dat[,1]),
                       SL.library = SL.lib,
                       family = binomial())
    
    kdpe(dat,Q_Y, Q_L1,g0,g1,
         density_bound = 0.001,
         converge_tol = 0.0001,
         lambda = lambda[i])

  }
stopCluster(my.cluster)
save(kdpe_data_lambda_tuning, file = "kdpe_data_lambda_tuning.RData")

for (j in 1:350){
  set.seed(j)
  dat <- generate_long_data(n)
  # Base estimators that need debiasing
  Q_Y <- SuperLearner(dat$Y,
                      dat[,1:4], 
                      SL.library = SL.lib, 
                      family = binomial())
  
  Q_L1 <- SuperLearner(dat$L1, 
                       dat[, 1:2], 
                       SL.library = SL.lib,
                       family = binomial())
  
  
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
  
  all_zero$L1 <- as.numeric(predict(Q_L1, newdata = all_zero[,1:2])$pred)
  all_zero$Y <- predict(Q_Y, newdata = all_zero[,1:4])$pred
  
  #mean(all_zero$Y)
  all_ones <- dat
  all_ones$A0 <- 1
  all_ones$A1 <- 1
  all_ones$L1 <- as.numeric(predict(Q_L1, newdata = all_ones[,1:2])$pred)
  all_ones$Y <- predict(Q_Y, newdata = all_ones[,1:4])$pred 
  mean(all_ones$Y)
  mu1_sl <- mean(all_ones$Y)
  mu0_sl <- mean(all_zero$Y)
  sl_data <- rbind(sl_data, c(mu1_sl-mu0_sl, 
                              mu1_sl/mu0_sl, 
                              (mu1_sl/(1-mu1_sl))/(mu0_sl/(1-mu0_sl))))
  
  
  
  
  tmle_results <- ltmle(data = dat, Anodes = c("A0", "A1"),
                        Lnodes = c("L1"),
                        Ynodes = c("Y"),
                        abar = list(c(1,1), c(0,0)),
                        survivalOutcome = F,
                        SL.library = SL.lib)
  
  tmle_data <- rbind(tmle_data, c(summary(tmle_results)[[2]]$ATE$estimate,
                                  summary(tmle_results)[[2]]$RR$estimate,
                                  summary(tmle_results)[[2]]$OR$estimate))
  
  
  kdpe_results <- kdpe(dat,Q_Y,
                       Q_L1,
                       g0,
                       g1,
                       density_bound = 0.001,
                       converge_tol = 0.0001,
                       lambda = 15)

  kdpe_data <- rbind(kdpe_data, kdpe_results)
  
  print((kdpe_data))
  
  
}
hist(kdpe_data$ATE)
hist(tmle_data$ATE)
hist(sl_data$ATE)

# values to create limit distributions
# code from Section 4 reused to generate the same plots

true_mu_11 <- 0.420
true_mu_00 <- 0.527

true_var_ate <- 0.00085









