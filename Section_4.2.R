# bootstrap confidence interval
source("helper.R")
library(SuperLearner)
library(tmle)
library(CVXR)

n <- 300
n_bootstrap <- 50
n_sims <- 15
dat <- generate_data(n)
SL.lib <- c("SL.randomForest", "SL.glm", "SL.mean")
plot_dat <- data.frame(method = NA, 
                       Value = NA,
                       estimate = NA,
                       lower = NA,
                       upper = NA)

set.seed(2023)

for (j in 1:n_sims){
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
    
    kdpe_results <- kdpe(Y,A,X,g,Q)
    bs_kdpe_ate[i] <- kdpe_results[1]
    bs_kdpe_rr[i] <- kdpe_results[2]
    bs_kdpe_or[i] <- kdpe_results[3]
    print(c(i, bs_kdpe_ate[i]))
  }
  
  bs_results <- data.frame(ATE = bs_kdpe_ate,
                           RR = bs_kdpe_rr,
                           OR = bs_kdpe_or)
  bs_results <- na.omit(bs_results)
  var_ests <- var(bs_results)
  
  
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
  
  Q_mat <- cbind(Q0W = predict(Q, newdata = data.frame(A=0,X))$pred,
                 Q1W = predict(Q, newdata = data.frame(A=1,X))$pred)
  g_est <- predict(g, newdata = data.frame(X))$pred
  
  tmle_results <- tmle(Y=dat$Y, A=dat$A, W=as.matrix(X), 
                       Q=Q_mat, 
                       g1W = g_est,
                       family = "binomial",
                       gbound = 0)
  
  # confidence intervals for results
  
  
  
  # ATE
  ate_kdpe <- kdpe_results[1]
  ate_kdpe_ci <- var_ests[1,1]
  ate_tmle <- tmle_results$estimates$ATE$psi
  ate_tmle_ci <- tmle_results$estimates$ATE$var.psi[1,1]
  
  plot_dat <- rbind(plot_dat, c("UDMLE", 
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
  
  plot_dat <- rbind(plot_dat, c("UDMLE", 
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
  
  plot_dat <- rbind(plot_dat, c("UDMLE", 
                                "OR", 
                                or_kdpe, 
                                or_kdpe - 1.96 * sqrt(or_kdpe_ci),
                                or_kdpe + 1.96 * sqrt(or_kdpe_ci)))
  plot_dat <- rbind(plot_dat, c("TMLE", 
                                "OR", 
                                or_tmle, 
                                tmle_results$estimates$OR$CI[1],
                                tmle_results$estimates$OR$CI[2]))
}


library(ggplot2)



plot_dat_full <- plot_dat[-1,]
plot_dat_full
plot_dat_full


# get table
nrow(plot_dat_full)/6

plot_dat_full$method <- as.factor(plot_dat_full$method)
plot_dat_full$value <- as.factor(plot_dat_full$value)
plot_dat_full$Value

plot_dat_ate <- plot_dat_full[plot_dat_full$Value == "ATE", ]


summary_plot_dat <- data.frame(method = NA,
                               value = NA,
                               estimate = NA,
                               lower = NA,
                               upper = NA)
# plot_dat_ate
for (val in unique(plot_dat_full$Value)){
  for (method in unique(plot_dat_full$method)){
    df_temp <- plot_dat_full[plot_dat_full$method == method & plot_dat_full$Value ==val, ]
    temp <- c(method,
              val,
              estimate = mean(as.numeric(df_temp$estimate)),
              lower = mean(as.numeric(df_temp$lower)),
              upper = mean(as.numeric(df_temp$upper)))
    summary_plot_dat <- rbind(summary_plot_dat, temp)
  }
}

names(summary_plot_dat) <- c("Method", "Parameter", "Avg. Estimate", "90% LB", "90% UB")


summary_plot_dat <- summary_plot_dat[-1,]


rownames(summary_plot_dat) <- NULL

xtable(summary_plot_dat)




