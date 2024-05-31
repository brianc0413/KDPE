# bootstrap confidence interval
source("helper.R")
library(SuperLearner)
library(tmle)
library(CVXR)
library(doParallel)

set.seed(2023)
ncores <- detectCores()-1
n <- 300
n_bootstrap <- 100
n_sims <- 240
SL.lib <- c("SL.randomForest", "SL.glm", "SL.mean")
plot_dat <- data.frame(method = NA, 
                       Value = NA,
                       estimate = NA,
                       lower = NA,
                       upper = NA)

for (j in 1:n_sims){
  dat <- generate_data(n)
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
  
  kdpe_results <- kdpe(Y,A,X,g,Q)
  
  # numerical issues
  if (is.na(kdpe_results[1])){
    next
  }
  
  Q_mat <- cbind(Q0W = predict(Q, newdata = data.frame(A=0,X))$pred,
                 Q1W = predict(Q, newdata = data.frame(A=1,X))$pred)
  g_est <- predict(g, newdata = data.frame(X))$pred
  
  tmle_results <- tmle(Y=dat$Y, A=dat$A, W=as.matrix(X), 
                       Q=Q_mat, 
                       g1W = g_est,
                       family = "binomial",
                       gbound = 0)
  
  
  # get bootstrap intervals
  # initialize parallelization
  my.cluster <- parallel::makeCluster(ncores, type = "PSOCK")
  print(my.cluster)
  registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  bs_results <- foreach(
    i = 1:n_bootstrap,
    .combine = rbind
  ) %dopar% {
    source("helper.R")
    library(SuperLearner)
    library(tmle)
    library(CVXR)
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
    kdpe(Y,A,X,g,Q)
  }
  # end parallelization
  stopCluster(my.cluster)
  
  bs_results <- bs_results[,1:3]
  rownames(bs_results) <- NULL
  colnames(bs_results) <- c("ATE", "RR", "OR")
  bs_results <- na.omit(bs_results)
  
  var_ests <- var(bs_results)
  # confidence intervals for results
  # ATE
  ate_kdpe <- kdpe_results[1]
  ate_kdpe_ci <- var_ests[1,1]
  ate_tmle <- tmle_results$estimates$ATE$psi
  ate_tmle_ci <- tmle_results$estimates$ATE$var.psi[1,1]
  
  plot_dat <- rbind(plot_dat, c("KDPE", 
                                "ATE", 
                                ate_kdpe, 
                                ate_kdpe - 1.96 * sqrt(ate_kdpe_ci), 
                                ate_kdpe + 1.96 * sqrt(ate_kdpe_ci)))
  plot_dat <- rbind(plot_dat, c("TMLE", 
                                "ATE", 
                                ate_tmle, 
                                tmle_results$estimates$ATE$CI[1], 
                                tmle_results$estimates$ATE$CI[2]))
  
  
  # relative risk
  rr_kdpe <- kdpe_results[2]
  rr_kdpe_ci <- var_ests[2,2]
  rr_tmle <- tmle_results$estimates$RR$psi
  
  plot_dat <- rbind(plot_dat, c("KDPE", 
                                "RR", 
                                rr_kdpe, 
                                rr_kdpe - 1.96 * sqrt(rr_kdpe_ci),
                                rr_kdpe + 1.96 * sqrt(rr_kdpe_ci)))
  
  plot_dat <- rbind(plot_dat, c("TMLE", 
                                "RR", 
                                rr_tmle, 
                                tmle_results$estimates$RR$CI[1],
                                tmle_results$estimates$RR$CI[2]))
  
  
  
  # OR
  or_kdpe <- kdpe_results[3]
  or_kdpe_ci <- var_ests[3,3]
  or_tmle <- tmle_results$estimates$OR$psi
  
  plot_dat <- rbind(plot_dat, c("KDPE", 
                                "OR", 
                                or_kdpe, 
                                or_kdpe - 1.96 * sqrt(or_kdpe_ci),
                                or_kdpe + 1.96 * sqrt(or_kdpe_ci)))
  plot_dat <- rbind(plot_dat, c("TMLE", 
                                "OR", 
                                or_tmle, 
                                tmle_results$estimates$OR$CI[1],
                                tmle_results$estimates$OR$CI[2]))
  print(plot_dat)
}

library(ggplot2)

plot_dat_full <- plot_dat[-1,]
plot_dat_full[,3:5] <- sapply(plot_dat_full[,3:5], as.numeric)

plot_dat_ate <- plot_dat_full[plot_dat_full$Value == "ATE", ]
plot_dat_rr <- plot_dat_full[plot_dat_full$Value == "RR", ]
plot_dat_or <- plot_dat_full[plot_dat_full$Value == "OR", ]

true_ate <- 0.1233
true_rr <- 0.52359/0.40026
true_or <- 0.52359/(1-0.52359)/(0.40026/(1-0.40026))
plot_dat_ate$cover <- (plot_dat_ate$lower <= true_ate) & 
  (plot_dat_ate$upper >= true_ate)
plot_dat_rr$cover <- (plot_dat_rr$lower <= true_rr) & 
  (plot_dat_rr$upper >= true_rr)
plot_dat_or$cover <- (plot_dat_or$lower <= true_or) & 
  (plot_dat_or$upper >= true_or)

mean(plot_dat_ate$cover[plot_dat_ate$method == "KDPE"])
mean(plot_dat_rr$cover[plot_dat_rr$method == "KDPE"])
mean(plot_dat_or$cover[plot_dat_or$method == "KDPE"])

# get estimated variances, and estimated average length of CIs
# use full plot_dat
plot_dat_full$var <- ((plot_dat_full$estimate - plot_dat_full$lower)/1.96)^2
df <- expand.grid(unique(plot_dat_full$method), unique(plot_dat_full$Value))
colnames(df) <- c("Method", "Parameter")
est_var <- rep(0, 6)
avg_length <- rep(0,6)
for (i in 1:nrow(df)){
  method <- df[i,1]
  param <- df[i,2]
  temp_df <- plot_dat_full[(plot_dat_full$method == method) &
                             (plot_dat_full$Value == param), ]
  est_var[i] <- mean(temp_df$var)
  avg_length[i] <- mean(temp_df$upper - temp_df$lower)
}
df$est_var <- est_var
df$avg_length <- avg_length
df
