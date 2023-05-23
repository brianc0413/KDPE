# bootstrap variance as m changes
source("helper.R")
library(SuperLearner)
library(tmle)
library(CVXR)

n <- 300
n_bootstrap <- 300
dat <- generate_data(n)
SL.lib <- c("SL.randomForest", "SL.glm", "SL.mean")
plot_dat <- data.frame(method = NA, 
                       Value = NA,
                       estimate = NA,
                       lower = NA,
                       upper = NA)

set.seed(2023)

# now run tmle, kdpe
Y <- dat$Y
A <- dat$A
X <- dat$X

Q <- SuperLearner(Y,
                  data.frame(A,X), 
                  SL.library = SL.lib, 
                  family = binomial())

g <- SuperLearner(A,
                  data.frame(X),
                  SL.library = SL.lib,
                  family = binomial())

kdpe_results <- kdpe(Y,A,X,g,Q, density_bound = 0.002, converge_tol = 0.001)

Q_mat <- cbind(Q0W = predict(Q, newdata = data.frame(A=0,X))$pred,
               Q1W = predict(Q, newdata = data.frame(A=1,X))$pred)
g_est <- predict(g, newdata = data.frame(X))$pred

tmle_results <- tmle(Y=dat$Y, A=dat$A, W=as.matrix(X), 
                     Q=Q_mat, 
                     g1W = g_est,
                     family = "binomial",
                     gbound = 0)



bs_kdpe_ate <- rep(0,n_bootstrap)
bs_kdpe_rr <- rep(0, n_bootstrap)
bs_kdpe_or <- rep(0, n_bootstrap)
  
for (i in 1:n_bootstrap){
  bstp_data <- dat[sample(nrow(dat),n,replace=TRUE),]
  Y <- bstp_data$Y
  A <- bstp_data$A
  X <- bstp_data$X
    
  Q <- SuperLearner(Y,
                    data.frame(A,X), 
                    SL.library = SL.lib, 
                    family = binomial())
    
  g <- SuperLearner(A,
                    data.frame(X),
                    SL.library = SL.lib,
                    family = binomial())
    
  bs_kdpe_results <- kdpe(Y,A,X,g,Q,
                          density_bound = 0.002,
                          converge_tol = 0.001)
  bs_kdpe_ate[i] <- bs_kdpe_results[1]
  bs_kdpe_rr[i] <- bs_kdpe_results[2]
  bs_kdpe_or[i] <- bs_kdpe_results[3]
  print(c(i, bs_kdpe_ate[i]))
}
  
bs_results <- data.frame(ATE = bs_kdpe_ate,
                           RR = bs_kdpe_rr,
                           OR = bs_kdpe_or)

#var_ests <- var(bs_results)


 ### Delete this part later
bs_results <- bs_results[1:158, ]
bs_results <- na.omit(bs_results)

bs_results <- bs_results - cbind(rep(mean(bs_results$ATE), nrow(bs_results)), rep(mean(bs_results$RR), nrow(bs_results)), rep(mean(bs_results$OR), nrow(bs_results)))
bs_results_2 <- bs_results_2 - cbind(rep(mean(bs_results_2$ATE), nrow(bs_results_2)), rep(mean(bs_results_2$RR), nrow(bs_results_2)), rep(mean(bs_results_2$OR), nrow(bs_results_2)))
bs_results <- rbind(bs_results, bs_results_2)
### Fill it it

nrow(bs_results)



bs_50 <- bs_results[1:50, ]
bs_100 <- bs_results[1:100, ]
bs_150 <- bs_results[1:150, ]
bs_200 <- bs_results[1:200, ]
bs_250 <- bs_results[1:250, ]
  

var_50 <- var(bs_50)
var_100 <- var(bs_100)
var_150 <- var(bs_150)
var_200 <- var(bs_200)
var_250 <- var(bs_250)

# get variances for tmle
tmle_ate_var <-  tmle_results$estimates$ATE$var.psi[1,1]
tmle_rr_var <- ((tmle_results$estimates$RR$CI[2]-tmle_results$estimates$RR$CI[1])/(1.96*2))^2
tmle_or_var <- ((tmle_results$estimates$OR$CI[2]-tmle_results$estimates$OR$CI[1])/(1.96*2))^2

ate_vars <- c(var_50[1,1], var_100[1,1], var_150[1,1], var_200[1,1], var_250[1,1], tmle_ate_var)
rr_vars <- c(var_50[2,2], var_100[2,2], var_150[2,2], var_200[2,2], var_250[2,2], tmle_rr_var)
or_vars <- c(var_50[3,3], var_100[3,3], var_150[3,3], var_200[3,3], var_250[3,3], tmle_or_var)

# confidence intervals for results
plot_df <- data.frame(Method = c(rep("KDPE", 5), "TMLE"), 
                      m = c(50, 100, 150, 200, 250, NA), 
                      ATE = ate_vars,
                      RR = rr_vars,
                      OR = or_vars)

plot_df
library("xtable")
print(xtable(plot_df, digits=5))
