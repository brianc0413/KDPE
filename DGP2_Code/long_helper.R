library(SuperLearner)
library(ltmle)
library(CVXR)
library(rje)

generate_long_data <- function(n){
  W <- runif(n, min = 0, max = 8)
  A0 <- rbinom(n, size=1, prob = 0.5) 
  L1 <- as.numeric(3 + A0 - 0.75*W + rnorm(n) > 0)
  A1 <- as.numeric(expit(-3 + 0.5*W + 0.4*L1) >= 0.3)
  Y <- as.numeric(expit((W-4) + 0.5-0.3* A0 - 
                          0.5*L1 - 0.5*A1 + rnorm(n)) >= 0.5)
  
  return(data.frame(W, A0, L1, A1, Y))
}


mu_00_sims <- function(n){
  W <- runif(n, min = 0, max = 8)
  A0 <- rep(0, times = n) 
  L1 <- as.numeric(3 + A0 - 0.75*W + rnorm(n) > 0)
  A1 <- rep(0, times=n)
  Y <- as.numeric(expit((W-4) + 0.5-0.3* A0 - 
                          0.5*L1 - 0.5*A1 + rnorm(n)) >= 0.5)
  
  return(mean(Y))
}

mu_11_sims <- function(n){
  W <- runif(n, min = 0, max = 8)
  A0 <- rep(1, times = n) 
  L1 <- as.numeric(3 + A0 - 0.75*W + rnorm(n) > 0)
  A1 <- rep(1, times=n)
  Y <- as.numeric(expit((W-4) + 0.5-0.3* A0 - 
                          0.5*L1 - 0.5*A1 + rnorm(n)) >= 0.5)
  
  return(mean(Y))
}

#' P is a table of conditional probabilities, which should represent 
#' P(Y, A_1, L1, A0| W), P(Y|A_1,L_1, A_0, W), P(L_1|A_0, W)
generate_P <- function(data,Q_Y,Q_L1,g0,g1){
  # table of conditional probabilities
  
  # first 2n rows are observed data, and its direct opposite
  n <- nrow(data)
  opp_data <- data
  opp_data$Y <- 1-data$Y
  P <- rbind(data, opp_data)
  
  # bottom rows never need updating, only needed for integral
  w <- data$W
  a0 <- c(0,1)
  l1 <- c(0,1)
  a1 <- c(0,1)
  y <- c(0,1)
  all_P <- expand.grid(w,a0,l1,a1,y)
  colnames(all_P) <- c("W", "A0", "L1", "A1", "Y")
  
  P <- rbind(P, anti_join(all_P, P))
  
  # calculate probabilities for all data
  g0_vec <- as.vector(predict(g0, data.frame(W= P[,1]))$pred)
  g0_vec[P$A0 == 0] <- 1-g0_vec[P$A0 == 0]
  l1_vec <- as.vector(predict(Q_L1, P[,1:2])$pred)
  l1_vec[P$L1 == 0] <- 1-l1_vec[P$L1 == 0]
  a1_vec <- as.vector(predict(g1, P[,1:3])$pred)
  a1_vec[P$A1 == 0] <- 1-a1_vec[P$A1==0]
  y_vec <- as.vector(predict(Q_Y, P[,1:4])$pred)
  y_vec[P$Y == 0] <- 1-y_vec[P$Y == 0]
  
  p <- g0_vec * l1_vec * a1_vec * y_vec
  P$p <- p
  P$p_y <- y_vec
  P$p_a1 <- a1_vec
  P$p_l1 <- l1_vec
  return(P)
}


# storing all exponential differences between possible observations
store_exp <- function(P){
  hold <- list()
  for (i in 1:nrow(P)){
    w <- P$W[i]
    a0 <- P$A0[i]
    l1 <- P$L1[i]
    a1 <- P$A1[i]
    y <- P$Y[i]
    temp <- sweep(P[,1:5], 2, c(w,a0,l1,a1,y))^2
    hold[[i]] <- exp(-1*rowSums(temp))
  }
  return(hold)
}

#' X is list of covariates
#' P is a table of conditional probabilities, which should represent P(Y,A|X)
#' f is a matrix of 16n values, representing the f function
f <- function(P, store_exp){
  f_vec <- rep(0, times = nrow(P))
  n <- nrow(P)/16
  p <- P$p
  for (i in 1:nrow(P)){
    f_vec[i] <- 1/n*sum(p*store_exp[[i]])
  }
  return(f_vec)
}

#' normalization constant for mean-zero term
norm_term <- function(P, f_vec){
  return(16/nrow(P) * sum(P$p * f_vec))
}

#' kernel function (mean zero)
#' norm term = should only be calculated once
# Construct Kernel Matrix
k_mat <- function(P, stored_exp, f_vec, norm){
  n <- nrow(P)
  kmatrix <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in i:n){
      base <- stored_exp[[i]][j]
      f_i <- f_vec[i]
      f_j <- f_vec[j]
      kmatrix[i,j] <- base - (f_i * f_j / norm)
      kmatrix[j,i] <- kmatrix[i,j]
    }
  }
  return(kmatrix)
}

# works - gives uniform update for all values (unobserved and observed) for P
proj_y_update <- function(dat, kmat, P, alpha_opt){
  n <- nrow(dat)
  #t((t(alpha_opt) %*% k_matrix) + 1) 
  
  # get the indices of vector where Y=1, Y=0
  zero_mat <- matrix(data= NA, nrow = n, ncol = nrow(P))
  one_mat <- matrix(data= NA, nrow = n, ncol = nrow(P))
  P_0 <- c()
  for (i in 1:n){
    w <- dat$W[i]
    a0 <- dat$A0[i]
    l1 <- dat$L1[i]
    a1 <- dat$A1[i]
    zero_ind <- which(P$W ==w & P$A0 == a0 & P$L1 == l1 & P$A1 == a1 & P$Y == 0)
    one_ind <- which(P$W ==w & P$A0 == a0 & P$L1 == l1 & P$A1 == a1 & P$Y == 1)
    
    
    zero_mat[i,] <- kmat[zero_ind, ]
    one_mat[i,] <- kmat[one_ind, ]
    
  }
  
  P_0 <- P$p_y
  P_0[P$Y == 1] <- 1-P$p_y[P$Y==1]
  
  cond_exp <- as.matrix(P_0, nrow=n, ncol=1) * t(t(alpha_opt) %*% zero_mat) + 
    as.matrix((1-P_0), ncol=n, ncol=1) * t(t(alpha_opt) %*% one_mat)
  
  
  return( t(t(alpha_opt) %*% kmat[1:n, ]) - cond_exp)
}

### finish updating this (i.e. fix probability weights for each term)
proj_l1_update <- function(dat, kmat, P, alpha){
  n <- nrow(dat)
  m <- nrow(P)
  
  mat_000 <- matrix(data= NA, nrow = n, ncol = m)
  mat_001 <- matrix(data= NA, nrow = n, ncol = m)
  mat_010 <- matrix(data= NA, nrow = n, ncol = m)
  mat_100 <- matrix(data= NA, nrow = n, ncol = m)
  mat_111 <- matrix(data= NA, nrow = n, ncol = m)
  mat_110 <- matrix(data= NA, nrow = n, ncol = m)
  mat_101 <- matrix(data= NA, nrow = n, ncol = m)
  mat_011 <- matrix(data= NA, nrow = n, ncol = m)
  
  P_000 <- c()
  P_001 <- c()
  P_010 <- c()
  P_100 <- c()
  P_111 <- c()
  P_110 <- c()
  P_101 <- c()
  P_011 <- c()
  
  
  mat_l1_00 <- matrix(data= NA, nrow = n, ncol = m)
  mat_l1_11 <- matrix(data= NA, nrow = n, ncol = m)
  mat_l1_10 <- matrix(data= NA, nrow = n, ncol = m)
  mat_l1_01 <- matrix(data= NA, nrow = n, ncol = m)
  P_l1_00 <- c()
  P_l1_11 <- c()
  P_l1_10 <- c()
  P_l1_01 <- c()
  
  # get the correct probability values
  for (i in 1:nrow(P)){
    w <- P$W[i]
    a0 <- P$A0[i]
    l1 <- P$L1[i]
    
    ind_000 <- which(P$W ==w & P$A0 == a0 & P$L1 == 0 & P$A1 == 0 & P$Y == 0)
    ind_001 <- which(P$W ==w & P$A0 == a0 & P$L1 == 0 & P$A1 == 0 & P$Y == 1)
    ind_010 <- which(P$W ==w & P$A0 == a0 & P$L1 == 0 & P$A1 == 1 & P$Y == 0)
    ind_100 <- which(P$W ==w & P$A0 == a0 & P$L1 == 1 & P$A1 == 0 & P$Y == 0)
    ind_111 <- which(P$W ==w & P$A0 == a0 & P$L1 == 1 & P$A1 == 1 & P$Y == 1)
    ind_110 <- which(P$W ==w & P$A0 == a0 & P$L1 == 1 & P$A1 == 1 & P$Y == 0)
    ind_101 <- which(P$W ==w & P$A0 == a0 & P$L1 == 1 & P$A1 == 0 & P$Y == 1)
    ind_011 <- which(P$W ==w & P$A0 == a0 & P$L1 == 0 & P$A1 == 1 & P$Y == 1)
    
    P_000 <- c(P_000, P$p_y[ind_000] * P$p_a1[ind_000] * P$p_l1[ind_000])
    P_001 <- c(P_001, P$p_y[ind_001] * P$p_a1[ind_001] * P$p_l1[ind_001])
    P_010 <- c(P_010, P$p_y[ind_010] * P$p_a1[ind_010] * P$p_l1[ind_010])
    P_100 <- c(P_100, P$p_y[ind_100] * P$p_a1[ind_100] * P$p_l1[ind_100])
    P_111 <- c(P_111, P$p_y[ind_111] * P$p_a1[ind_111] * P$p_l1[ind_111])
    P_110 <- c(P_110, P$p_y[ind_110] * P$p_a1[ind_110] * P$p_l1[ind_110])
    P_101 <- c(P_101, P$p_y[ind_101] * P$p_a1[ind_101] * P$p_l1[ind_101])
    P_011 <- c(P_011, P$p_y[ind_011] * P$p_a1[ind_011] * P$p_l1[ind_011])
    
    if (l1 == 1){
      P_l1_00 <- c(P_l1_00, P$p_y[ind_100] * P$p_a1[ind_100])
      P_l1_11 <- c(P_l1_11, P$p_y[ind_111] * P$p_a1[ind_111])
      P_l1_10 <- c(P_l1_10, P$p_y[ind_110] * P$p_a1[ind_110])
      P_l1_01 <- c(P_l1_01, P$p_y[ind_101] * P$p_a1[ind_101])
    }
    else{
      P_l1_00 <- c(P_l1_00, P$p_y[ind_000] * P$p_a1[ind_000])
      P_l1_11 <- c(P_l1_11, P$p_y[ind_011] * P$p_a1[ind_011])
      P_l1_10 <- c(P_l1_10, P$p_y[ind_010] * P$p_a1[ind_010])
      P_l1_01 <- c(P_l1_01, P$p_y[ind_001] * P$p_a1[ind_001])
    }
    
  }
  
  # get the correct kernel matrix
  for (i in 1:n){
    w <- dat$W[i]
    a0 <- dat$A0[i]
    l1 <- dat$L1[i]
    
    ind_000 <- which(P$W ==w & P$A0 == a0 & P$L1 == 0 & P$A1 == 0 & P$Y == 0)
    ind_001 <- which(P$W ==w & P$A0 == a0 & P$L1 == 0 & P$A1 == 0 & P$Y == 1)
    ind_010 <- which(P$W ==w & P$A0 == a0 & P$L1 == 0 & P$A1 == 1 & P$Y == 0)
    ind_100 <- which(P$W ==w & P$A0 == a0 & P$L1 == 1 & P$A1 == 0 & P$Y == 0)
    ind_111 <- which(P$W ==w & P$A0 == a0 & P$L1 == 1 & P$A1 == 1 & P$Y == 1)
    ind_110 <- which(P$W ==w & P$A0 == a0 & P$L1 == 1 & P$A1 == 1 & P$Y == 0)
    ind_101 <- which(P$W ==w & P$A0 == a0 & P$L1 == 1 & P$A1 == 0 & P$Y == 1)
    ind_011 <- which(P$W ==w & P$A0 == a0 & P$L1 == 0 & P$A1 == 1 & P$Y == 1)
    

    
    mat_000[i,] <- kmat[ind_000, ]
    mat_001[i,] <- kmat[ind_001, ]
    mat_010[i,] <- kmat[ind_010, ]
    mat_100[i,] <- kmat[ind_100, ]
    mat_111[i,] <- kmat[ind_111, ]
    mat_110[i,] <- kmat[ind_110, ]
    mat_101[i,] <- kmat[ind_101, ]
    mat_011[i,] <- kmat[ind_011, ]
    
    
    if (l1 == 1){
      mat_l1_00[i,] <- kmat[ind_100, ]
      mat_l1_11[i,] <- kmat[ind_111, ]
      mat_l1_10[i,] <- kmat[ind_110, ]
      mat_l1_01[i,] <- kmat[ind_101, ]
      
      
    }
    else{
      mat_l1_00[i,] <- kmat[ind_000, ]
      mat_l1_11[i,] <- kmat[ind_011, ]
      mat_l1_10[i,] <- kmat[ind_010, ]
      mat_l1_01[i,] <- kmat[ind_001, ]
    }
    
  }
  
  term <- as.matrix(P_l1_00, nrow = m, ncol=1) * t( t(alpha) %*% mat_l1_00) + 
    as.matrix(P_l1_11, nrow=m, ncol=1) * t( t(alpha) %*% mat_l1_11) +
    as.matrix(P_l1_01,nrow=m, ncol=1) * t( t(alpha) %*% mat_l1_01) + 
    as.matrix(P_l1_10, nrow=m, ncol=1) * t( t(alpha) %*% mat_l1_10)
  
  cond_exp <- P_000 * t( t(alpha) %*% mat_000) + P_001 * t( t(alpha) %*% mat_001) +
    P_010 * t( t(alpha) %*% mat_010) + P_100 * t( t(alpha) %*% mat_100) + 
    P_100 * t( t(alpha) %*% mat_100) + P_111 * t( t(alpha) %*% mat_111) + 
    P_110 * t( t(alpha) %*% mat_110) + P_101 * t( t(alpha) %*% mat_101) + 
    P_011 * t( t(alpha) %*% mat_011)  
  return(term - cond_exp)
}

# function for calculating means given treatments
get_mean <- function(data, P, treat_vec){
  #index_y1 <- which(P$A0 == treat_vec[1] & 
  #                      P$A1 == treat_vec[2] & P$Y == 1)
  #index_all <- which(P$A0 == treat_vec[1] & 
  #                      P$A1 == treat_vec[2])
  #return(sum(P$p[index_y1])/sum(P$p[index_all]))

  
  w <- data$W
  value <- 0
  for (i in 1:nrow(data)){
    index_l1_0 <- which(P$W == w[i] & P$A0 == treat_vec[1] & 
                          P$A1 == treat_vec[2] & P$L1 ==0 & P$Y == 1)
    index_l1_1 <- which(P$W == w[i] & P$A0 == treat_vec[1] & 
                          P$A1 == treat_vec[2] & P$L1 ==1 & P$Y == 1)
    value <- value + 
      ( (P$p_l1[index_l1_0] + (1-P$p_l1[index_l1_1]))/2 * P$p_y[index_l1_0] + 
         (1-P$p_l1[index_l1_0] + P$p_l1[index_l1_1])/2 * P$p_y[index_l1_1] )
  }
  return(value/nrow(data))
}
# kdpe code
kdpe <- function(data,Q_Y,Q_L1,g0,g1,
                 density_bound = 0.001,converge_tol = 0.001,lambda = 0.001){
  
  # initialize values
  n <- nrow(data)
  
  # construct kernel matrix
  P <- generate_P(data,Q_Y,Q_L1,g0,g1)
  stored_values <- store_exp(P)

  
  # set up optimization problem
  
  alpha_opt <- rep(100000, n)
  iter_count <- 0
  change <- 1E6
  
  # to keep track of multiplicative factor
  
  while (change > converge_tol){
    iter_count <- iter_count+1
    
    # iteration too much then terminate
    if (iter_count > 100){
      return(NA)
    }
    # construct necessary kernel functions
    f_vec <- f(P, stored_values)
    norm <- norm_term(P, f_vec)
    kmat <- k_mat(P, stored_values, f_vec, norm)
    p_term <- P$p[1:n]
    alpha <- Variable(n)
    K <- kmat[1:n, 1:n]
    obj <- (sum(log((t(alpha) %*% K + 1) * t(p_term) ))-lambda*quad_form(alpha, K)) 
    
    # constraint: must be a distribution
    py_constr <- t((t(alpha) %*% kmat[1:n,1:n]) + 1) * P$p[1:n]
    constr <- list(py_constr >= density_bound, 1-py_constr >= density_bound)
    prob <- Problem(Maximize(obj), constr)
    result <- psolve(prob)
    
    if (result$status == "solver_error"){
      return(NA)
    }
    
    
    alpha_opt <- as.vector(result$getValue(alpha)) # kinda crazy values 
    # get final values (update values)
    # P(Y=y|A=a,X=x) update
    
    h_l1 <- proj_l1_update(data, kmat, P, alpha_opt)
    h_y <- proj_y_update(data, kmat, P, alpha_opt)
    
    # update distribution P$p
    old_p <- P$p
    old_p_y <- P$p_y
    old_p_l1 <- P$p_l1
    P$p_y <- P$p_y * (h_y + 1)
    P$p_l1 <- P$p_l1 * (h_l1 + 1)
    P$p <- P$p/(old_p_y * old_p_l1) * P$p_y * P$p_l1
    
    change <- mean((P$p - old_p)^2)
    print(c("Alpha_Sqaured:", sum(alpha_opt^2)))
    print(c("Change:", change))
    
  }
  
  # get final estimates
  mu_11 <- get_mean(data, P, c(1,1))
  mu_00 <- get_mean(data, P, c(0,0))
  
  ate_estimates <- mu_11-mu_00
  rr_estimates <- mu_11/mu_00
  or_estimates <- (mu_11/(1-mu_11))/(mu_00/(1-mu_00))
  return(c(ate_estimates, rr_estimates, or_estimates, iter_count, lambda))
}














