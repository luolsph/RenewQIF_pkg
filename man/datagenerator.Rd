\name{datagenerator}
\alias{datagenerator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate data streams with correlated subjects within each single data batch.
}
\description{
  Generate input data for \code{renewqif}, may include two abnormal data batches.
}
\usage{
datagenerator(n, m, p, B, tempdatadir, type, beta0, beta1, beta2, cpt1, cpt2, 
  intercept, categorical, corst_x, rho_x, corst_y, rho_y, phi_y, seed)
}
\arguments{
  \item{n}{
  Sample size for each single data batch.
}
  \item{m}{
  Number of repeated measurements or cluster size.
}
  \item{p}{
  Number of covariates (including intercept). 
}
  \item{B}{Number of data batches/streams.}
  \item{tempdatadir}{Directory that stores the input data streams.}
  \item{type}{
  Type/Family of response variable, c("gaussian", "binomial", "poisson").
}
  \item{beta0}{
  Main underlying true regression coefficients, length equals to p.
}
  \item{beta1}{
  Coefficients for the first abnormal data batch, if it exists. Otherwise set beta1 = beta0.
}
  \item{beta2}{
  Coefficients for the second abnormal data batch, if it exists. Otherwise set beta2 = beta0.
}
  \item{cpt1}{
  Index for the first abnormal data batch.
}
  \item{cpt2}{
  Index for the second abnormal data batch.
  }
  \item{intercept}{
  Logical, if TRUE, beta0[1] will be the coefficient for intercept.
  }
  \item{categorical}{logical, if TRUE, the first covariate of X will be changed to categorical variable.}
  \item{corst_x}{Correlation structure of X, c("ind", "cs", "ar1").}
  \item{rho_x}{Correlation coefficient for the covariance of X.}
  \item{corst_y}{Correlation structure of y, c("ind", "cs", "ar1").}
  \item{rho_y}{Correlation coefficient for the covariance of y (outcome vector for each cluster).}
  \item{phi_y}{Dispersion parameter for generating y.}
  \item{seed}{random seed for simulation.}
}


\examples{


}

\keyword{ Cluster-correlated data streams }
\keyword{ Simulation }