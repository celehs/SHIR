#' SHIR_fit: aggregate the locally derived summary data using the proposed method in SHIR paper.
#'
#' @param H_lst list of locally derived Hessian matrix
#' @param d_lst list of locally derived gradient vector
#' @param n_lst vector of the sample sizes at local sites
#' @param lambda_lst candidate set of the tuning parameter for mu,
#' (corresponds to sqrt(N)lambda in the SHIR paper)
#' @param lambda_g_lst candidate set of the tuning parameter for alpha,
#' (corresponds to sqrt(N)*lambdalambda_g in the SHIR paper)
#' @param tune the information criterion used for model selection:
#' options include 'AIC', 'BIC' and 'RIC', which put different weights on the degree of freedom.
#'
#' @return min.lambda: the selected tuning parameters with the minimum information criterion
#' @return min.beta: matrices of the fitted SHIR estimator
#' (the m-th column loads the fitted coefficient (beta) for the m-th site,
#' the first row is the fitted intercepts)
#' @export
#'
SHIR_fit <- function(H_lst, d_lst, n_lst,
                     lambda_lst, lambda_g_lst, tune = 'BIC'){

  options(warn = -1)
  intercept <- T

  M <- length(d_lst)
  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  n <- sum(n_lst)
  p <- length(H_lst[[1]][1, ]) - 1
  initial <- rep(0, M * (p + 1))

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

  X_all <- X_lst[[1]]
  Y_all <- Y_lst[[1]]
  for (i in 2:M){
    X_all <- cbind(X_all, - X_lst[[1]])
  }

  for (i in 2:M){
    X_i <- X_lst[[i]]
    for (j in 2:M){
      if (j == i){
        X_i <- cbind(X_i, X_lst[[i]])
      }else{
        X_i <- cbind(X_i, matrix(0, length(X_i[ ,1]), length(X_lst[[1]][1, ])))
      }
    }
    X_all <- rbind(X_all, X_i)
    Y_all <- c(Y_all, Y_lst[[i]])
  }

  if (intercept == T){
    indx_lst <- c(c(NA, 1:p),
                  rep(c(NA, (p + 1):(2 * p)), M - 1))
  }

  # tune with GIC
  min.lambda <- NULL
  min.GIC <- Inf
  min.coef <- NULL
  min.beta <- NULL

  total_num <- length(lambda_g_lst) * length(lambda_lst)
  t <- 0
  for (lambda in lambda_lst){
    for (lambda_g in lambda_g_lst){
      penscale <- function(x){
        if (x == 1){
          return(1)
        }else{
          return(lambda_g / lambda)
        }
      }
      invisible(utils::capture.output(
        fit_result <- grplasso::grplasso(X_all, y = Y_all, index = indx_lst, standardize = F,
                               center = FALSE, lambda = lambda, penscale = penscale, model = grplasso::LinReg(),
                               coef.init = initial)))
      fit_coef <- fit_result$coefficients
      fit_norm <- fit_result$norms.pen
      initial <- fit_coef
      result_GIC <- Cal_GIC_pool(d_lst, H_lst, X_all, fit_coef, n_lst, lambda_g / 2,
                                 intercept = T, type = tune)
      GIC <- result_GIC$GIC
      if (GIC < min.GIC){
        min.lambda <- c(lambda, lambda_g)
        min.GIC <- GIC
        min.beta <- result_GIC$beta
        min.coef <- fit_coef
      }
      t <- t + 1
      print(paste0('Finishing: ', 100 * round(t / total_num, 4), '%'))
    }
  }
  return(list(min.lambda = min.lambda, min.beta = min.beta))
}
