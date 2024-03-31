library(quadprog)
library(Matrix)
library(mvtnorm)
library(ggplot2)
set.seed(654321)


lambda <- 3.0
skill <- 0.15
sim_runs <- 1000
period <- 252
runs <- 10000


Assets_df <- read.csv("Sample_AA.csv", header = FALSE)
Assets <- Assets_df$V1
nvar <- as.numeric(length(Assets))
sim_ER <- as.vector(Assets_df[ , 2])/period
sim_EV <- as.vector(Assets_df[ , 3])/sqrt(period)
sim_Cor <- as.matrix(Assets_df[ , 4:ncol(Assets_df)])
sim_Cov <- diag(sim_EV) %*% sim_Cor %*% diag(sim_EV)
sim_Cov <- as.matrix(nearPD(sim_Cov)$mat)
full_ret <- rmvnorm(n = runs * period * 2, mean = sim_ER, sigma = sim_Cov)


idx_mat <- matrix(0, nrow = runs, ncol = 4)
for (i in 1:nrow(idx_mat)) {
  idx_mat[i, 1] <- 1 + ((2 * period) * (i - 1))
  idx_mat[i, 2] <- idx_mat[i, 1] + period - 1
  idx_mat[i, 3] <- idx_mat[i, 2] + 1
  idx_mat[i, 4] <- idx_mat[i, 3] + period - 1
} 


Results <- as.data.frame(matrix(0, nrow = runs, ncol = nvar + 3))
names(Results) <- c(Assets, "Return", "Vol", "Ratio")
runs_mat <- matrix(0, nrow = sim_runs, ncol = nvar)


reg_vec <- as.vector(c(1, 0.9, 0.8, 0.7, 0.6, 0.5))
med_SR_vec <- numeric(length(reg_vec))


for (j in 1:length(reg_vec)) {
  resamp <- reg_vec[j]
  for (i in 1:nrow(Results)) {
    back_slice <- as.matrix(full_ret[idx_mat[i, 1]:idx_mat[i, 2], ])
    back_ER <- as.vector(colMeans(back_slice) * period)
    back_EV <- as.vector(apply(X = back_slice, MARGIN = 2, FUN = sd) * sqrt(period))
    back_Cor <- as.matrix(cor(back_slice))
    
    forw_slice <- as.matrix(full_ret[idx_mat[i, 3]:idx_mat[i, 4], ])
    forw_ER <- as.vector(colMeans(forw_slice) * period)
    forw_EV <- as.vector(apply(X = forw_slice, MARGIN = 2, FUN = sd) * sqrt(period))
    forw_Cor <- as.matrix(cor(forw_slice))
    
    final_ER <- as.vector(skill * forw_ER + (1 - skill) * back_ER)
    final_EV <- as.vector(skill * forw_EV + (1 - skill) * back_EV)
    final_Cor <- as.matrix(skill * forw_Cor + (1 - skill) * back_Cor)
    final_Cov <- diag(final_EV) %*% final_Cor %*% diag(final_EV)
    final_Cov <- as.matrix(nearPD(final_Cov)$mat)
    
    for (k in 1:nrow(runs_mat)) {
      picks <- sample(1:nvar, nvar * resamp)
      final_ER2 <- final_ER
      final_ER2[-picks] <- -10
      Dmat <- lambda * final_Cov
      dvec <- final_ER2
      Amat <- t(rbind(t(rep(1, nvar)), diag(nvar)))
      bvec <- c(1, rep(0, nvar))
      runs_mat[k, ] <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1, factorized = FALSE)$solution
    }
    
    Results[i, (1:nvar)] <- colMeans(runs_mat)
    
    ret_vec <- forw_slice %*% as.numeric(Results[i, (1:nvar)])
    Results[i, (nvar + 1)] <- mean(ret_vec - forw_slice[ , 5]) * period
    Results[i, (nvar + 2)] <- sd(ret_vec) * sqrt(period)
    Results[i, (nvar + 3)] <- Results[i, (nvar + 1)] / Results[i, (nvar + 2)]
  }
  med_SR_vec[j] <- median(Results$Ratio)
}



################################################################################



resamp <- 0.5

for (i in 1:nrow(Results)) {
  back_slice <- as.matrix(full_ret[idx_mat[i, 1]:idx_mat[i, 2], ])
  back_ER <- as.vector(colMeans(back_slice) * period)
  back_EV <- as.vector(apply(X = back_slice, MARGIN = 2, FUN = sd) * sqrt(period))
  back_Cor <- as.matrix(cor(back_slice))
  
  forw_slice <- as.matrix(full_ret[idx_mat[i, 3]:idx_mat[i, 4], ])
  forw_ER <- as.vector(colMeans(forw_slice) * period)
  forw_EV <- as.vector(apply(X = forw_slice, MARGIN = 2, FUN = sd) * sqrt(period))
  forw_Cor <- as.matrix(cor(forw_slice))
  
  final_ER <- as.vector(skill * forw_ER + (1 - skill) * back_ER)
  final_EV <- as.vector(skill * forw_EV + (1 - skill) * back_EV)
  final_Cor <- as.matrix(skill * forw_Cor + (1 - skill) * back_Cor)
  final_Cov <- diag(final_EV) %*% final_Cor %*% diag(final_EV)
  final_Cov <- as.matrix(nearPD(final_Cov)$mat)
  
  for (k in 1:nrow(runs_mat)) {
    picks <- sample(1:nvar, nvar * resamp)
    final_ER2 <- final_ER
    final_ER2[-picks] <- -10
    Dmat <- lambda * final_Cov
    dvec <- final_ER2
    Amat <- t(rbind(t(rep(1, nvar)), diag(nvar)))
    bvec <- c(1, rep(0, nvar))
    runs_mat[k, ] <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1, factorized = FALSE)$solution
  }
  
  Results[i, (1:nvar)] <- colMeans(runs_mat)
  
  ret_vec <- forw_slice %*% as.numeric(Results[i, (1:nvar)])
  Results[i, (nvar + 1)] <- mean(ret_vec - forw_slice[ , 5]) * period
  Results[i, (nvar + 2)] <- sd(ret_vec) * sqrt(period)
  Results[i, (nvar + 3)] <- Results[i, (nvar + 1)] / Results[i, (nvar + 2)]
}


write.csv(Results, "res_0.5.csv")
