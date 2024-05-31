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
# Base estimators that need debiasing

tmle_data <- data.frame(ATE=NA, RR=NA, OR=NA)
set.seed(100)
SL.lib <- c("SL.randomForest", "SL.glm", "SL.mean")

for (i in 1:200){
  set.seed(i)
  dat <- generate_long_data(800)
  tmle_results <- ltmle(data = dat, Anodes = c("A0", "A1"),
                        Lnodes = c("L1"),
                        Ynodes = c("Y"),
                        abar = list(c(1,1), c(0,0)),
                        survivalOutcome = F,
                        SL.library = SL.lib)
  
  tmle_data <- rbind(tmle_data, c(summary(tmle_results)[[2]]$ATE$estimate,
                                  summary(tmle_results)[[2]]$RR$estimate,
                                  summary(tmle_results)[[2]]$OR$estimate))
}


tmle_data_800 <- tmle_data[-1,]
save(tmle_data_800, file="tmle_data_800.RData")

my.cluster <- parallel::makeCluster(ncores, type = "PSOCK")
print(my.cluster)
registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

library(SuperLearner)

temp_data <- foreach(
  i = 1:250,
  .combine = rbind) %dopar% {
    set.seed(i)
    source("long_helper.R")
    library(SuperLearner)
    library(ltmle)
    library(CVXR)
    dat <- generate_long_data(800)
    tmle_results <- ltmle(data = dat, Anodes = c("A0", "A1"),
                          Lnodes = c("L1"),
                          Ynodes = c("Y"),
                          abar = list(c(1,1), c(0,0)),
                          survivalOutcome = F,
                          SL.library = SL.lib)
    
    c(summary(tmle_results)[[2]]$ATE$estimate,
      summary(tmle_results)[[2]]$RR$estimate,
      summary(tmle_results)[[2]]$OR$estimate)
}

stopCluster(my.cluster)
rownames(temp_data) <- NULL
colnames(temp_data) <- c("ate", "rr", "or")
temp_data <- data.frame(temp_data)


#hist(tmle_data$ATE)
save(temp_data, file = "tmle_800_data.RData")
