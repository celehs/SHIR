test_that("Local_fit works", {
  data(simu_data)
  X_lst <- simu_data$X
  Y_lst <- simu_data$Y
  n <- 400
  p <- 400
  M <- length(X_lst)
  length_lst <- c()
  I_lst <- vector('list', M)
  U_lst <- vector('list', M)
  beta_lst <- vector('list', M)

  for (m in 1:M){
    Y <- Y_lst[[m]]
    X <- X_lst[[m]]

    local_summary_m <- Local_fit(Y, X, lambda_lst = 0.1 * c(4:20) * sqrt((log(p)) / n))

    I_lst[[m]] <- local_summary_m$hessian
    U_lst[[m]] <- local_summary_m$gradient
    beta_lst[[m]] <- local_summary_m$beta
    length_lst <- c(length_lst, length(Y))
  }

  expect_equal(round(I_lst[[m]][1,1]), 83)
})


test_that("SHIR_fit works", {
  data(simu_data)
  X_lst <- simu_data$X
  Y_lst <- simu_data$Y
  n <- 400
  p <- 400
  M <- length(X_lst)
  length_lst <- c()
  I_lst <- vector('list', M)
  U_lst <- vector('list', M)
  beta_lst <- vector('list', M)

  for (m in 1:M){
    Y <- Y_lst[[m]]
    X <- X_lst[[m]]

    local_summary_m <- Local_fit(Y, X, lambda_lst = 0.1 * c(4:20) * sqrt((log(p)) / n))

    I_lst[[m]] <- local_summary_m$hessian
    U_lst[[m]] <- local_summary_m$gradient
    beta_lst[[m]] <- local_summary_m$beta
    length_lst <- c(length_lst, length(Y))
  }

  SHIR_train <- SHIR_fit(I_lst, U_lst, length_lst)

  expect_equal(round(SHIR_train$min.beta[1],2), -0.02)
})
