#' Aggregates the locally derived summary data using the proposed method SHIR.
#'
#' @param H_lst List of locally derived Hessian matrices.
#' @param d_lst List of locally derived gradient vectors.
#' @param n_lst Vector containing the sample size at each local site.
#' @param lambda_lst Candidate set of the tuning parameters for mu.
#' It corresponds to sqrt(N)lambda in the SHIR paper.
#' If not specified or specified as NULL by the user, default value is 0.3*c(5:25)sqrt(nlog(p)).
#'
#' @param lambda_g_lst Candidate set of the tuning parameter for alpha.
#' It corresponds to lambda_g in the SHIR paper.
#' If not specified or specified as NULL by the user, default value is c(0.6, 0.9).
#'
#' @param tune Information criterion used for model selection.
#' Input options include 'AIC', 'BIC', 'mBIC' and 'RIC', which put different weights on the degree of freedom.
#' If not specified by the user, default value is 'BIC'.
#'
#' @importFrom Matrix bdiag
#' @importFrom grplasso grplasso LinReg
#'
#' @return min.lambda: Selected tuning parameters with the minimum information criterion.
#' @return min.beta: Matrices of the fitted SHIR estimator.
#' The m-th column contains the fitted coefficient beta of the m-th site, and
#' the first row contains the fitted intercepts.
#' @export
#'
SHIR_fit <- function(H_lst, d_lst, n_lst, lambda_lst = NULL,
                     lambda_g_lst = NULL, tune = 'BIC'){

  options(warn = -1)
  M <- length(d_lst)
  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  n <- sum(n_lst)
  p <- length(H_lst[[1]][1, ]) - 1


  if (is.null(lambda_lst)){
    lambda_lst = sqrt(n * log(p)) * 0.3 * c(5:25)
  }
  if (is.null(lambda_g_lst)){
    lambda_g_lst = c(0.6, 0.9)
  }

  initial <- rep(0, (M + 1) * (p + 1))

  for (i in 1:M){
    H <- H_lst[[i]]
    d <- d_lst[[i]]
    mat_all <- cbind(H, d)
    mat_all <- rbind(mat_all, t(c(d, max(mat_all) + n)))
    svd_result <- svd(mat_all)
    s_value <- svd_result$d
    s_mat <- diag(sqrt(s_value))[ ,1:(min(p + 1, n_lst[[i]]) + 1)]
    data_all <- svd_result$u %*% s_mat
    X <- t(data_all[-length(data_all[ ,1]),])
    Y <- data_all[length(data_all[ ,1]),]
    X_lst[[i]] <- X
    Y_lst[[i]] <- Y
  }

  X_all <- c()
  Y_all <- c()
  for (i in 1:M){
    X_all <- rbind(X_all, X_lst[[i]])
    Y_all <- c(Y_all, Y_lst[[i]])
  }

  X_alpha <- X_lst[[1]]
  for (i in 2:M){
    X_alpha <- bdiag(X_alpha, X_lst[[i]])
  }
  X_all <- cbind(X_all, X_alpha)

  enlarge_mat <- matrix(0, p + 1, p + 1)
  for (i in 1:M) {
    enlarge_mat <- cbind(enlarge_mat, diag(rep(1, p + 1)))
  }

  indx_lst <- c(c(NA, 1:p), rep(c(NA, (p + 1):(2 * p)), M))

  # tune with GIC
  min.lambda <- NULL
  min.GIC <- Inf
  min.coef <- NULL
  min.beta <- NULL

  lambda_all <- c()
  lambda_g_all <- c()
  fit_coef_mat <- c()

  total_num <- length(lambda_g_lst) * length(lambda_lst)
  z <- 0

  for (lambda in lambda_lst){
    for (lambda_g in lambda_g_lst){

      lambda_g <- lambda_g * lambda
      penscale <- function(x){
        ifelse(x == 1, 1, lambda_g / lambda)
      }

      Y_ <- 0
      rho_ <- n
      enlarge_Y <- rep(- Y_ / sqrt(rho_ / 2) / 2, p + 1)
      X_enlarge <- rbind(X_all, sqrt(rho_ / 2) * enlarge_mat)
      Y_enlarge <- c(Y_all, enlarge_Y)

      eta <- 2
      for (iter in 1:10) {
        invisible(utils::capture.output(
          fit_result <- grplasso::grplasso(x = as.matrix(X_enlarge), y = Y_enlarge,
                                           index = indx_lst, standardize = F,
                                           center = FALSE, lambda = lambda,
                                           penscale = penscale, model = LinReg(),
                                           coef.init = initial)))
        fit_coef <- fit_result$coefficients
        initial <- fit_coef

        alpha_fit <- c()
        for (i in 2:(M + 1)){
          alpha_fit <- cbind(alpha_fit, fit_coef[((i - 1) * p + i):(i * p + i)])
        }

        Y_ <- Y_ + rho_ * rowMeans(alpha_fit)
        rho_ <- rho_ * eta
        X_enlarge <- rbind(X_all, sqrt(rho_ / 2) * enlarge_mat)
        Y_enlarge <- c(Y_all, - Y_ / sqrt(rho_ / 2) / 2 * rep(1, p + 1))

      }

      lambda_all <- c(lambda_all, lambda)
      lambda_g_all <- c(lambda_g_all, lambda_g)
      fit_coef_mat <- cbind(fit_coef_mat, fit_coef)

      result_GIC <- Cal_GIC_pool(d_lst, H_lst, X_enlarge,
                                 fit_coef, n_lst, lambda_g / 2, type = tune)
      GIC <- result_GIC$GIC
      if (GIC < min.GIC){
        min.lambda <- c(lambda, lambda_g / lambda)
        min.GIC <- GIC
        min.beta <- result_GIC$beta
        min.coef <- fit_coef
      }
      z <- z + 1
      print(paste0('Finishing: ', 100 * round(z / total_num, 4), '%'))
    }
  }
  return(list(min.lambda = min.lambda, min.beta = min.beta))
}
