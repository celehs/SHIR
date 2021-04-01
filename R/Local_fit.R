#' Local_fit: Fit lasso and derive summary statistics at each local site.
#'
#' @param Y individual response vector
#' @param X individual covariates matrices
#' @param lambda_lst candidate set of tuning parameters, NULL for choose the default range in glmnet
#'
#' @return hessian: derived Hessian matrix
#' @return gradient: derived gradient vector
#' @return beta: local lasso estimator
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
  #pi_vec <- as.vector(X_all %*% beta_fit)
  grad <- t(X_all) %*% (Y - pi_vec)
  I_mat <- t(X_all) %*% diag(pi_vec * (1- pi_vec)) %*% X_all
  U <- I_mat %*% beta_fit + grad

  return(list(hessian = I_mat, gradient = U, beta = beta_fit))
}
