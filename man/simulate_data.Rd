% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slasso_functions.R
\name{simulate_data}
\alias{simulate_data}
\title{Simulate data through the function-on-function linear regression model}
\usage{
simulate_data(scenario, n_obs = 3000, type_x = "Bspline")
}
\arguments{
\item{scenario}{A  character strings indicating the scenario considered. It could be "Scenario I", "Scenario II", "Scenario III", and "Scenario IV".}

\item{n_obs}{Number of observations.}

\item{type_x}{Covariate generating mechanism, either Bspline or Brownian.}
}
\value{
A list containing the following arguments:

\code{X}: Covariate matrix, where  the rows  correspond to argument values and columns to replications.

\code{Y}: Response matrix, where  the rows  correspond to argument values and columns to replications.

\code{X_fd}: Coavariate functions.

\code{Y_fd}: Response functions.

\code{clus}: True cluster membership vector.
}
\description{
Generate synthetic data as in the simulation study of Centofanti et al. (2020).
}
\examples{
library(slasso)
data<-simulate_data("Scenario II",n_obs=150)
}
\references{
Centofanti, F., Fontana, M., Lepore, A., & Vantini, S. (2020).
Smooth LASSO Estimator for the Function-on-Function Linear Regression Model.
\emph{arXiv preprint arXiv:2007.00529}.
}
