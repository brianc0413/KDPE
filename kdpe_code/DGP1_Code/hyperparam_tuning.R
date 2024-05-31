# Hyperparameter Tuning
library(SuperLearner)
library(tmle)
library(CVXR)
source("helper.R")


n <- 300
n_sims <- 100

density_bound <- c(0.05, 0.01, 0.002)
gamma <- c(0.025, 0.005, 0.001)

hyperparams <- expand.grid(density_bound, gamma)
names(hyperparams) <- c("density_bound", "gamma")


tmle_data <- data.frame(db = NA,
                        gam = NA,
                        ate = NA,
                        rr = NA,
                        or = NA)
sl_data <- data.frame(db = NA,
                      gam = NA,
                        ate = NA,
                        rr = NA,
                        or = NA)

kdpe_data <- data.frame(db = NA,
                         gam = NA,
                        ate = NA,
                        rr = NA,
                        or = NA,
                        iter_count = NA)

SL.lib <- c("SL.randomForest", "SL.glm", "SL.mean")

for (i in 1:nrow(hyperparams)){
  
  db <- hyperparams$density_bound[i]
  gam <- hyperparams$gamma[i]
  
  for (i in 1:n_sims){
    dat <- generate_data(n)
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
    
    
    Q_mat <- cbind(Q0W = predict(Q, newdata = data.frame(A=0,X))$pred,
                   Q1W = predict(Q, newdata = data.frame(A=1,X))$pred)
    g_est <- predict(g, newdata = data.frame(X))$pred
    
    # no tmle fix
    
    sl_ate <- mean(Q_mat[,2]-Q_mat[,1])
    sl_rr <- mean(Q_mat[,2])/mean(Q_mat[,1])
    mu1sl <- mean(Q_mat[,2])
    mu0sl <- mean(Q_mat[,1])
    sl_or <- (mu1sl/(1-mu1sl))/(mu0sl/(1-mu0sl))
    sl_data <- rbind(sl_data, c(db, gam, sl_ate, sl_rr, sl_or))
    
    # tmle estimate
    tmle_results <- tmle(Y=dat$Y, A=dat$A, W=as.matrix(X), 
                         Q=Q_mat, 
                         g1W = g_est,
                         family = "binomial",
                         gbound = 0)
    tmle_ate <- tmle_results$estimates$ATE$psi
    tmle_rr <- tmle_results$estimates$RR$psi
    tmle_or <- tmle_results$estimates$OR$psi
    tmle_data <- rbind(tmle_data, c(db, gam, tmle_ate, tmle_rr, tmle_or))
    
    kdpe_results <- kdpe(Y,A,X,g,Q, density_bound = db, converge_tol = gam)
    kdpe_ate <- kdpe_results[1]
    kdpe_rr <- kdpe_results[2]
    kdpe_or <- kdpe_results[3]
    iter_counter <- kdpe_results[4]
    kdpe_data <- rbind(kdpe_data, c(db, gam, kdpe_ate, kdpe_rr, kdpe_or, iter_counter))
    
    
    print(kdpe_data)
  }
  
  
}


ate_true <- 1/3 - 0.3 + 0.09
rr_true <- 0.52359/0.40026
true_or <- 0.52359/(1-0.52359)/(0.40026/(1-0.40026))

hyperparams$failures <- rep(0, nrow(hyperparams))
hyperparams$bias_ate <- rep(0, nrow(hyperparams))
hyperparams$bias_rr <- rep(0, nrow(hyperparams))
hyperparams$bias_or <- rep(0, nrow(hyperparams))
hyperparams$avg_iters <- rep(0, nrow(hyperparams))

for (i in 1:nrow(hyperparams)){
  temp_data <- kdpe_data[(kdpe_data$db == hyperparams[i,1])
                          & (kdpe_data$gam == hyperparams[i,2]) , ]
  hyperparams$failures[i] <- mean(is.na(temp_data$ate))
  
  temp_data_clean <- na.omit(temp_data)
  hyperparams$bias_ate[i] <- mean(temp_data_clean$ate) - ate_true
  hyperparams$bias_rr[i] <- mean(temp_data_clean$rr) - rr_true
  hyperparams$bias_or[i] <- mean(temp_data_clean$or) - true_or
  hyperparams$avg_iters[i] <- mean(temp_data_clean$iter_count)
  
}
hyperparams

library(xtable)
print(xtable(hyperparams, digits = 5))


# plot ate as illustrative check for normality
kdpe_data_cleaned <- na.omit(kdpe_data)

library(ggplot2)
for (i in 1:nrow(hyperparams)){
  temp_data <- kdpe_data_cleaned[(kdpe_data_cleaned$db == hyperparams[i,1])
                          & (kdpe_data_cleaned$gam == hyperparams[i,2]) , ]
  temp_plot <- ggplot(temp_data, aes(x=ate))+
    geom_histogram(aes(y=..density..), 
                   lwd = 1, 
                   color = "black", 
                   fill = "white")+
    theme_bw() + 
    geom_vline(xintercept = ate_true, color = "blue", linetype = "dashed", lwd = 2) + 
    geom_density(alpha=.2, fill="#BB4444") + xlim(-0.05, 0.3) + ylim(0, 10)
  plot(temp_plot)

  
}


library(ggplot2)




