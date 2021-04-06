# Required packages:

#library(glmnet)
#library(grplasso)
#library(MASS)

# Internal function
AR_cov <- function(p, acorr){
  cov_mat <- diag(rep(1, p))
  cor_vec <- acorr^c(1:p)
  for (i in 1:(p - 1)){
    cov_mat[i, (i + 1):p] <- cor_vec[1: (p - i)]
    cov_mat[(i + 1):p, i] <- cor_vec[1: (p - i)]
  }
  return(cov_mat)
}


# Generate_data: generate an individual level data example for trial.

### Input ###

# n: sample size at each site
# p: dimensionality
# M: number of sites
# magn_mu: magnitude of the homogeneous effect
# magn_alpha: magnitude of the heterogeneous effect


### Output ###

# X_lst: list of individual covariates matrices of the local sites.
# Y_lst: list of individual response vectors of the local sites.
# beta_true: list of true coefficients (beta) of the local sites.
# mu: vector of the true mu (homogeneous effect).
# alpha: list of true alpha (heterogeneous effect) of the local sites.

#' @importFrom stats rbinom
Generate_data <- function(n = 400, p = 400, M = 4, magn_mu = 0.5, magn_alpha = 0.35){

  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  r_lst = rep(0.3, M)

  # Generate the coefficients:

  s <- 6
  sign_lst <- c(1, -1, 1, -1, 1, -1)
  center_support <- c(1:6)
  beta_center <- rep(0, p)
  beta_center[center_support] <- sign_lst * rep(magn_mu, s)
  beta_center <- c(0, beta_center)

  center_support <- center_support + 1
  noise_support <- c(4:9)

  neg_indx <- sample(c(1:M), as.integer(M / 2))
  beta_delta_lst <- vector('list', M)
  beta_true_lst <- vector('list', M)

  sign_lst <-  c(1, 1, 1, -1, -1, -1)

  for (m in 1:M) {
    beta_delta_lst[[m]] <- rep(0, (p + 1))
    for (t in 1:length(noise_support)){
      sign_t <- sign_lst[t]
      if (m %in% neg_indx){
        beta_delta_lst[[m]][noise_support[t]] <- - sign_t * magn_alpha
      }else{
        beta_delta_lst[[m]][noise_support[t]] <- sign_t * magn_alpha
      }
    }

    beta_true_lst[[m]] <- beta_center + beta_delta_lst[[m]]
  }

  beta_center <- 0
  for (m in 1:M){
    beta_center <- beta_center + beta_true_lst[[m]]
  }
  beta_center <- beta_center / M
  for (m in 1:M){
    beta_delta_lst[[m]] <- beta_true_lst[[m]] - beta_center
  }


  # Generate X and Y.

  prop <- 15 / p
  for (m in 1:M){
    r <- r_lst[m]
    gamma_coef <- c()
    for (t in 1:s) {
      gamma_coef <- cbind(gamma_coef, r * stats::rbinom(p - s, 1, prop))
    }
    Sigma_remain <- AR_cov(p - s, r)
    Sigma_X <- rbind(cbind(diag(rep(1, s)) + t(gamma_coef) %*% Sigma_remain %*% gamma_coef,
                           t(Sigma_remain %*% gamma_coef)),
                     cbind(Sigma_remain %*% gamma_coef, Sigma_remain))
    X <- MASS::mvrnorm(n, rep(0, p), Sigma_X)
    odds <- exp(cbind(rep(1, n), X) %*% beta_true_lst[[m]])
    p_vec <- odds / (1 + odds)
    Y <- rbinom(n, 1, p_vec)
    X_lst[[m]] <- X
    Y_lst[[m]] <- Y
  }

  return(list(X = X_lst, Y = Y_lst, beta = beta_true_lst,
              mu = beta_center, alpha = beta_delta_lst))
}


#' @import Matrix
# Function to tune (internal)
Cal_GIC_pool <- function(U_lst, I_lst, X_enlarge, fit_coef, length_lst, lambda_g,
                         type = 'BIC'){
  M <- length(U_lst)
  p <- length(I_lst[[1]][1, ]) - 1
  mu <- fit_coef[1:(p + 1)]
  beta_fit <- c()
  alpha_fit <- c()

  for (i in 2:(M + 1)){
    beta_fit <- cbind(beta_fit, mu + fit_coef[((i - 1) * p + i):(i * p + i)])
    alpha_fit <- cbind(alpha_fit, fit_coef[((i - 1) * p + i):(i * p + i)])
  }

  norm_lst <- sqrt(rowSums(alpha_fit^2))
  S_alpha <- which(norm_lst != 0)
  S_mu <- which(mu != 0)
  S_full <- which(fit_coef != 0)
  # H_S <- t(X_enlarge[ ,S_full]) %*% X_enlarge[ ,S_full]
  return(list('M'=X_enlarge[ ,S_full], 'S'=S_full))
  H_S <- crossprod(X_enlarge[ ,S_full])
  partial <- diag(0, length(S_full), length(S_full))
  if (length(S_alpha) != 0){
    for (t in 1:length(S_alpha)) {
      j <- S_alpha[t]
      partial_j <- t + length(S_mu) + length(S_alpha) * c(0:(M - 1))
      partial[partial_j, partial_j] <- lambda_g / norm_lst[j] *
        (diag(1, M, M) - alpha_fit[j,] %*% t(alpha_fit[j,]) / norm_lst[j]^2)
    }
  }

  df <- sum(diag(solve(H_S + partial) %*% H_S))
  print(df)

  GIC <- 0
  for (i in 1:M){
    GIC <- GIC + t(beta_fit[,i]) %*% I_lst[[i]] %*% beta_fit[,i] - 2 * t(beta_fit[,i]) %*% U_lst[[i]]
  }

  if (type == 'BIC'){
    GIC <- GIC / sum(length_lst) + df * log(sum(length_lst)) / sum(length_lst)
  }
  if (type == 'mBIC'){
    GIC <- GIC / sum(length_lst) + df * log(max(exp(1), log(M * p))) * log(sum(length_lst)) / sum(length_lst)
  }
  if (type == 'AIC'){
    GIC <- GIC / sum(length_lst) + df * 2 / sum(length_lst)
  }
  if (type == 'RIC'){
    GIC <- GIC / sum(length_lst) + df * log(M * p) / sum(length_lst)
  }
  return(list(GIC = GIC, beta = beta_fit))
}

