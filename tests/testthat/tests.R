test_that("Local_fit works", {
  data(simu_data)
  X_lst <- simu_data$X
  Y_lst <- simu_data$Y
  n <- 400
  p <- 400

  # number of sites:
  M <- length(X_lst)

  # vector for sample size of the local sites:
  length_lst <- c()

  # list for derived hessian matrices:
  I_lst <- vector('list', M)

  # list for derived gradient vectors:
  U_lst <- vector('list', M)

  # list for derived lasso estimators:
  beta_lst <- vector('list', M)

  for (m in 1:M){
    Y <- Y_lst[[m]]
    X <- X_lst[[m]]

    # Locally fit and derive the individual data at site m:

    local_summary_m <- Local_fit(Y, X, lambda_lst = 0.1 * c(3:30) * sqrt((log(p)) / n))

    I_lst[[m]] <- local_summary_m$hessian
    U_lst[[m]] <- local_summary_m$gradient
    beta_lst[[m]] <- local_summary_m$beta
  }

  expect_equal(round(I_lst[[m]][1,1],3), 78.315)
})


test_that("SHIR_fit works", {
  data(simu_data)
  X_lst <- simu_data$X
  Y_lst <- simu_data$Y
  n <- 400
  p <- 400

  # number of sites:
  M <- length(X_lst)

  # vector for sample size of the local sites:
  length_lst <- c()

  # list for derived hessian matrices:
  I_lst <- vector('list', M)

  # list for derived gradient vectors:
  U_lst <- vector('list', M)

  # list for derived lasso estimators:
  beta_lst <- vector('list', M)

  for (m in 1:M){
    Y <- Y_lst[[m]]
    X <- X_lst[[m]]

    # Locally fit and derive the individual data at site m:

    local_summary_m <- Local_fit(Y, X, lambda_lst = 0.1 * c(3:30) * sqrt((log(p)) / n))

    I_lst[[m]] <- local_summary_m$hessian
    U_lst[[m]] <- local_summary_m$gradient
    beta_lst[[m]] <- local_summary_m$beta
  }

  SHIR_train <- SHIR_fit(I_lst, U_lst, length_lst,
                         lambda_lst = sqrt(n * log(p)) * 0.3 * c(5:30),
                         lambda_g_lst = sqrt(n * (log(p) + M)) * 0.3 * (5:30),
                         tune = 'BIC')

  expect_equal(round(SHIR_train$min.beta[1],3), -0.021)
})
