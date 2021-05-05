#' Fits lasso and derives summary statistics at each local site.
#'
#' @param Y Individual response vector.
#' @param X Individual covariates matrix.
#' @param lambda_lst Candidate set of tuning parameters.
#' Default value is range in glmnet.
#'
#' @return hessian: Derived Hessian matrix.
#' @return gradient: Derived gradient vector.
#' @return beta: Local lasso estimator.
#'
#' @export
Local_fit <- function(Y, X, lambda_lst = NULL){
  cv.result <- glmnet::cv.glmnet(X, Y, family = 'binomial', lambda = lambda_lst)
  lambda.cv <- cv.result$lambda.min
  model <- glmnet::glmnet(X, Y, family = 'binomial', lambda = lambda.cv)
  beta_fit <- c(as.vector(model$a0), as.vector(model$beta))

  n <- length(Y)
  X_all <- cbind(rep(1, n), X)
  pi_vec <- as.vector(1 / (1 + exp(- X_all %*% beta_fit)))
  grad <- t(X_all) %*% (Y - pi_vec)
  I_mat <- t(X_all) %*% diag(pi_vec * (1- pi_vec)) %*% X_all
  U <- I_mat %*% beta_fit + grad

  return(list(hessian = I_mat, gradient = U, beta = beta_fit))
}

