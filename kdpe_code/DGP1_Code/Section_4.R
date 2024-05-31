# kdpe main simulation (Section 4.1)
library(SuperLearner)
library(tmle)
library(CVXR)
library(doParallel)
library(foreach)
source("helper.R")

# Initialize Settings 
n <- 300
n_sims <- 500

tmle_ate <- rep(0,n_sims)
tmle_rr <- rep(0, n_sims)
tmle_or <- rep(0, n_sims)

sl_ate <- rep(0,n_sims)
sl_rr <- rep(0, n_sims)
sl_or <- rep(0, n_sims)

kdpe_ate <- rep(0,n_sims)
kdpe_rr <- rep(0, n_sims)
kdpe_or <- rep(0, n_sims)

SL.lib <- c("SL.randomForest", "SL.glm", "SL.mean")

for (i in 1:n_sims){
  library(SuperLearner)
  library(tmle)
  library(CVXR)
  
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
  
  # sl fix
  sl_ate[i] <- mean(Q_mat[,2]-Q_mat[,1])
  sl_rr[i] <- mean(Q_mat[,2])/mean(Q_mat[,1])
  mu1sl <- mean(Q_mat[,2])
  mu0sl <- mean(Q_mat[,1])
  sl_or[i] <- (mu1sl/(1-mu1sl))/(mu0sl/(1-mu0sl))
  
  # tmle estimate
  tmle_results <- tmle(Y=dat$Y, A=dat$A, W=as.matrix(X), 
                       Q=Q_mat, 
                       g1W = g_est,
                       family = "binomial",
                       gbound = 0)
  tmle_ate[i] <- tmle_results$estimates$ATE$psi
  tmle_rr[i] <- tmle_results$estimates$RR$psi
  tmle_or[i] <- tmle_results$estimates$OR$psi
  
  kdpe_results <- kdpe(Y,A,X,g,Q)
  kdpe_ate[i] <- kdpe_results[1]
  kdpe_rr[i] <- kdpe_results[2]
  kdpe_or[i] <- kdpe_results[3]
  
  
  print(c(i, kdpe_ate[i]))
  print(c(i, kdpe_rr[i]))
  
}


# Making Density Plots/ QQ Plots to Test Asymptotic Normality

# what should the actual ate be? 0.123333333
ate_true <- 1/3 - 0.3 + 0.09

ate_estimates <- data.frame(kdpe = kdpe_ate,
                        TMLE = tmle_ate,
                        SL = sl_ate)

rr_estimates <- data.frame(kdpe  = kdpe_rr,
                           TMLE = tmle_rr,
                           SL = sl_rr)

or_estimates <- data.frame(kdpe  = kdpe_or,
                           TMLE = tmle_or,
                           SL = sl_or)



ate_estimates <- na.omit(ate_estimates)
rr_estimates <- na.omit(rr_estmates)
or_estimates <- na.omit(or_estimates)


# Density Plots / QQ plots 

df <- data.frame(Method = c(rep("KDPE", nrow(ate_estimates)), rep("TMLE", nrow(ate_estimates)), rep("SL", nrow(ate_estimates))),
                 ATE = c(ate_estimates$kdpe, ate_estimates$TMLE, ate_estimates$SL))

library(ggplot2)
#ggplot(df, aes(x=ATE, fill=Method)) +
#  geom_histogram( color='#e9ecef', alpha=0.6, position='identity', bins = 10)

ggplot(df, aes(x=ATE, fill = Method, color = Method, ..count..))+
  geom_density(alpha=0.3, position = "identity", lwd = 1.2, n= 15)+
  theme_bw() + geom_vline(xintercept = ate_true, color = "red", linetype = "dashed")

mean(ate_estimates$kdpe)
mean(ate_estimates$TMLE)
mean(ate_estimates$SL)




df_2 <- data.frame(Method = c(rep("kdpe", nrow(rr_estimates)), rep("TMLE", nrow(rr_estimates)), rep("SL", nrow(rr_estimates))),
                 RR = c(rr_estimates$kdpe, rr_estimates$TMLE, rr_estimates$SL))

ggplot(df_2, aes(x=RR, fill = Method, color = Method, ..count..))+
  geom_density(alpha=0.3, position = "identity", lwd = 1.2, n= 15)+
  theme_bw()+ geom_vline(xintercept = 0.52359/0.40026, color="red", linetype = "dashed")


df_3 <- data.frame(Method = c(rep("kdpe", nrow(or_estimates)), rep("TMLE", nrow(or_estimates)), rep("SL", nrow(or_estimates))),
                   OR = c(or_estimates$kdpe, or_estimates$TMLE, or_estimates$SL))

true_or <- 0.52359/(1-0.52359)/(0.40026/(1-0.40026))

ggplot(df_3, aes(x=OR, fill = Method, color = Method, ..count..))+
  geom_density(alpha=0.3, position = "identity", lwd = 1.2, n= 15)+
  theme_bw()+ geom_vline(xintercept = true_or, color="red", linetype = "dashed")


# Histograms
ggplot(df, aes(x=ATE, fill = Method, color = Method, ..count..))+
  geom_density(alpha=0.3, position = "identity", lwd = 1.2, n= 15)+
  theme_bw() + geom_vline(xintercept = ate_true, color = "red", linetype = "dashed")
ggplot(df_2, aes(x=RR, fill = Method, color = Method, ..count..))+
  geom_density(alpha=0.3, position = "identity", lwd = 1.2, n= 15)+
  theme_bw()+ geom_vline(xintercept = 0.52359/0.40026, color="red", linetype = "dashed")
ggplot(df_3, aes(x=OR, fill = Method, color = Method, ..count..))+
  geom_density(alpha=0.3, position = "identity", lwd = 1.2, n= 15)+
  theme_bw()+ geom_vline(xintercept = true_or, color="red", linetype = "dashed")

# QQ Plots
ggplot(df, aes(sample=ATE,fill= Method, color = Method))+
  stat_qq() + 
  stat_qq_line() + theme_bw() +  labs(y= "Sample Quantiles for ATE", x = "Theoretical Quantiles")
ggplot(df_2, aes(sample=RR,fill= Method, color = Method))+
  stat_qq() + 
  stat_qq_line() + theme_bw() +  labs(y= "Sample Quantiles for RR", x = "Theoretical Quantiles")
ggplot(df_3, aes(sample=OR,fill= Method, color = Method))+
  stat_qq() + 
  stat_qq_line() + theme_bw() + labs(y= "Sample Quantiles for OR", x = "Theoretical Quantiles")
  
# ATE Histograms (compared to limit distribution)
var_ate_true <- 0.002515497
# true variance and ate
true_var <- var(if_vec)/450
true_ate <- ate_true

dat_KDPE <- df[df$Method=="KDPE", 2]
x_fit <- seq(-0.1, 
             0.4, length = 40)
y_fit <- dnorm(x_fit, mean = true_ate, sd = sqrt(true_var))
x_breaks <- seq(from = -1*0.1,to= 0.4,length.out= 20)

hist(dat_KDPE, prob = TRUE,
     xlab = "KDPE ATE Estimates", 
     breaks = x_breaks,
     xlim = c(-0.1, 0.4), main = "KDPE ATE Distribution",
     col = "red", ylim = c(0, 7.5))
lines(x_fit, y_fit, col = "purple", lwd = 2)+ 
  legend("topright", c("KDPE", "Limit"), fill=c("red", "purple"))


dat_tmle <- df[df$Method=="TMLE",2]
hist(dat_tmle, prob = TRUE,
     xlab = "TMLE ATE Estimates", 
     breaks = x_breaks,
     xlim = c(-0.1, 0.4), main = "TMLE ATE Distribution",
     col = "blue", ylim = c(0, 7.5))
lines(x_fit, y_fit, col = "purple", lwd = 2)+ 
  legend("topright", c("TMLE", "Limit"), fill=c("blue", "purple"))


dat_sl <- df[df$Method=="SL",2]
hist(dat_sl, prob = TRUE,
     xlab = "SL ATE Estimates", 
     breaks = x_breaks,
     xlim = c(-0.1, 0.4), main = "SL ATE Distribution",
     col = "green", ylim = c(0, 7.5))
lines(x_fit, y_fit, col = "purple", lwd = 2)+ 
  legend("topright", c("SL", "Limit"), fill=c("green", "purple"))


# Dataframe for printing table
table_df <- data.frame(ATE = c(mean(ate_estimates$kdpe),
                               mean(ate_estimates$TMLE),
                               mean(ate_estimates$SL)),
                       RR = c(mean(rr_estimates$kdpe),
                              mean(rr_estimates$TMLE),
                              mean(rr_estimates$SL)),
                       OR = c(mean(or_estimates$kdpe),
                              mean(or_estimates$TMLE),
                              mean(or_estimates$SL)))
print(table_df)


# QQ plots to show normality
qqnorm(ate_estimates$kdpe)
qqnorm(ate_estimates$TMLE)
qqnorm(ate_estimates$SL)

qqnorm(rr_estimates$kdpe)
qqnorm(rr_estimates$TMLE)
qqnorm(rr_estimates$SL)

qqnorm(or_estimates$kdpe)
qqnorm(or_estimates$TMLE)
qqnorm(or_estimates$SL)

table_df$ATE <- table_df$ATE - ate_true
table_df$RR <- table_df$RR - 0.52359/0.40026
table_df$OR <- table_df$OR- true_or

print(table_df)

