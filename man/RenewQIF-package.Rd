\name{RenewQIF-package}
\alias{RenewQIF-package}
\alias{RenewQIF}
\docType{package}
\title{
\packageTitle{RenewQIF}
}
\description{
\packageDescription{RenewQIF}
}
\details{This package is used to conduct real-time regression analysis of streaming
clustered data with possible abnormal data batches.
}
\author{
\packageAuthor{RenewQIF}

Maintainer: \packageMaintainer{RenewQIF}
}
\references{
Real-time Regression Analysis of Streaming Clustered Data with Possible Abnormal 
Data Batches.
}

\keyword{Incremental statistical analysis} 
\keyword{Quadratic inference function}
\keyword{Quality control}
\keyword{Streaming datasets}


\examples{
NB <- 10000  # total sample size
B <- 100  # number of data batches
n <- round(NB/B) # data batch size
m <- 10  # cluster size

thred <- 0.05  # threshold of statistical significance level
QC <- FALSE  # specify whether to conduct the monitoring procedure

cpt1 <- round(B * 0.25)  
cpt2 <- round(B * 0.75)
beta0 <- c(1, 1)
beta1 <- beta2 <- beta0

p <- length(beta0)

type <- family <-"gaussian"
corst_x <- "cs"
rho_x <- 0.5

intercept <- TRUE
categorical <- FALSE

## True correlation matrix 
corst_y <- "cs"

rho_y <- 0.7
phi_y <- 1
seed <- 1
## Working matrix
corstr <- "exchangeable"

tempdatadir <- "~/Desktop/tempdata"

datagenerator(n, m, p, B, tempdatadir, type, beta0, beta1, beta2, cpt1, cpt2, 
	intercept = intercept, categorical = categorical, 
	corst_x = corst_x, rho_x = rho_x, corst_y = corst_y, rho_y = rho_y, phi_y = phi_y, 
	seed = seed)

out <- renewqif(B = B, tempdatadir = tempdatadir, family = family, intercept = intercept, corstr = corstr, 
    	QC = QC, thred = thred)

unlink(tempdatadir)

# Output the regression coefficient estimates without abnormal data batch
out$beta

}
