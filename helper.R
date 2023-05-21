library(SuperLearner)
library(tmle)
library(CVXR)
library(foreach)

generate_data <- function(n){
  X <- runif(n)
  A <- rbinom(n,
              size = 1,
              prob = 0.5+1/3*sin(50*X/pi))
  Y <- rbinom(n,
              size = 1,
              prob = 0.4 + A*(X-0.3)^2 + 1/4*sin(40*X/pi))
  return(data.frame(Y,A,X))
}



#' P is a table of conditional probabilities, which should represent P(Y,A|X)
generate_P <- function(Y,A,X,g,Q){
  # table of conditional probabilities
  n <- length(X)
  P <- data.frame(Y=rep(c(rep(0,n), rep(1,n)) , 2), 
                  A = c(rep(0, 2*n), rep(1, 2*n)),
                  X = rep(X,4), 
                  p = 0)
  p <- rep(NA, nrow(P))
  
  # estimate g, Q values for observed data
  g_vec <- predict(g, newdata = data.frame(X=P[,c(3)]))$pred
  Q_vec <- predict(Q, newdata = P[,c(2,3)])$pred
  
  p_a <- g_vec
  p_a[P$A==0] <- 1-g_vec[P$A==0]  
  p_y <- Q_vec
  p_y[P$Y==0] <- 1-Q_vec[P$Y==0] 
  
  # get conditional probabilities
  p <- p_a*p_y
  P$p <- p
  return(P)
}

#' X is list of covariates
#' P is a table of conditional probabilities, which should represent P(Y,A|X)
f <- function(y,a,x,P){
  n <- length(X)
  p <- P$p
  hold <- sweep(P[,1:3], 2, c(y,a,x))^2
  dist <- exp(-1*rowSums(hold))
  value<- 1/n*sum(p*dist)
  return(value)
}

#' normalization constant for mean-zero term
norm_term <- function(P){
  f_vals <- rep(0,nrow(P))
  for (i in 1:nrow(P)){
    y <- P[i,1]
    a <- P[i,2]
    x <- P[i,3]
    f_vals[i] <- f(y,a,x,P)
  }
  
  return(4/nrow(P) * sum(P$p * f_vals))
}

#' kernel function (mean zero)
#' v1 = tuple of y,a,x
#' v2 = tuple of y,a,x
#' norm term = should only be calculated once
k <- function(v1, v2, norm, P){
  
  k_base <- exp(-1 * sum((v1-v2)^2))
  
  y_1 <- v1[1]
  a_1 <- v1[2]
  x_1 <- v1[3]
  
  y_2 <- v2[1]
  a_2 <- v2[2]
  x_2 <- v2[3]
  
  return(k_base - f(y_1,a_1,x_1, P)*f(y_2,a_2,x_2, P)/norm)
}


k_mat <- function(Y,A,X, norm, P){
  n <- length(X)
  kmatrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in i:n){
      v1 <- c(Y[i], A[i], X[i])
      v2 <- c(Y[j], A[j], X[j])
      kmatrix[i,j] <- k(v1, v2, norm, P)
    }
  }
  
  for (i in 2:n){
    for (j in 1:(i-1)){
      kmatrix[i,j] <- kmatrix[j,i]
    }
  }
  
  return(kmatrix)
}

# # check to see if this is mean zero
# v1 <- unlist(P[1, 1:3])
# k_hold <- rep(0, nrow(P))
# for (i in 1:nrow(P)){
#   v2 <- unlist(P[i, 1:3])
#   k_hold[i] <- k(v1,v2, norm, P)
# }
# sum(k_hold * P$p) ## indeed mean 0!!!! so we're good to go!!!


# atmle code: regularized approach with estimator
kdpe <- function(Y,A,X,g,Q, density_bound = 0.01, converge_tol = 0.001){
  
  # initialize values
  n <- length(X)
  
  # construct kernel matrix
  P <- generate_P(Y,A,X,g,Q)
  
  # get Q function for values
  Q0 <- predict(Q, newdata = data.frame(A=0,X))$pred
  Q1 <- predict(Q, newdata = data.frame(A=1,X))$pred
  
  # set up optimization problem
  Q_base <- predict(Q, newdata = data.frame(A,X))$pred
  p_y <- Q_base
  p_y[Y==0] <- 1-Q_base[Y==0]
  
  alpha_opt <- rep(100000, n)
  iter_count <- 0
  change <- 1E6
  #while (mean(alpha_opt^2) > 1/n){
  while (change > converge_tol){
    iter_count <- iter_count+1
    
    # iteration too much then terminate
    if (iter_count > 100){
      return(NA)
    }
    # construct necessary kernel functions
    norm <- norm_term(P)
    k_matrix <- k_mat(Y,A,X, norm, P)
    k_matrix <- k_matrix[which(!duplicated(k_matrix)), ]
    n_unique <- nrow(k_matrix)
    
    alpha <- Variable(n_unique)
    obj <- sum(log((t(alpha) %*% k_matrix + 1) * t(p_y) ) ) 
    
    # constraint: must be a distribution
    py_constr <- t((t(alpha) %*% k_matrix) + 1) * p_y
    constr <- list(py_constr >= density_bound, 1-py_constr >= density_bound)
    prob <- Problem(Maximize(obj), constr)
    result <- psolve(prob)
    
    if (result$status == "solver_error"){
      return(NA)
    }
    
    
    alpha_opt <- as.vector(result$getValue(alpha)) # kinda crazy values 
    # get final values (update values)
    # P(Y=y|A=a,X=x) update
    old_p_y <- p_y
    p_y <- p_y * t((t(alpha_opt) %*% k_matrix) + 1) 
    # P(Y=1 |A=1,X=x) update ... what about the unobserved quantities?
    # as of right now, not updated
    Q1_update <- Q1
    Q1_update[Y==1 & A==1] <- p_y[Y==1 & A==1]
    Q1_update[Y==0 & A==1] <- (1-p_y[Y==0 & A==1])
    
    Q0_update <- Q0  
    Q0_update[Y==1 & A==0] <- p_y[Y==1 & A==0]
    Q0_update[Y==0 & A==0] <- 1-p_y[Y==0 & A==0]
    
    # Update P accordingly
    for (i in 1:n){
      index_matched <- (P$A==A[i]) & (P$X == X[i]) & (P$Y == Y[i])
      index_opposite <- (P$A==A[i]) & (P$X == X[i]) & (P$Y == 1-Y[i])
      P[index_matched,4] <- p_y[i]
      P[index_opposite,4] <- 1-p_y[i]
    }
    
    change <- mean((p_y - old_p_y)^2)
    print(c("Norm:", sum(alpha_opt^2)))
    print(c("Change:", change))
    
  }
  ate_estimates <- mean(Q1_update-Q0_update)
  rr_estimates <- mean(Q1_update)/mean(Q0_update)
  mu1 <- mean(Q1_update)
  mu0 <- mean(Q0_update)
  or_estimates <- (mu1/(1-mu1))/(mu0/(1-mu0))
  return(c(ate_estimates, rr_estimates, or_estimates))
}

