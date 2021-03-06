% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Local_fit.R
\name{Local_fit}
\alias{Local_fit}
\title{Fits lasso and derives summary statistics at each local site.}
\usage{
Local_fit(Y, X, lambda_lst = NULL)
}
\arguments{
\item{Y}{Individual response vector.}

\item{X}{Individual covariates matrix.}

\item{lambda_lst}{Candidate set of tuning parameters.
Default value is range in glmnet.}
}
\value{
hessian: Derived Hessian matrix.

gradient: Derived gradient vector.

beta: Local lasso estimator.
}
\description{
Fits lasso and derives summary statistics at each local site.
}
