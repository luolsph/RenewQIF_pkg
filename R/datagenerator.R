library(MASS)
library(SimCorMultRes)
##########################################
## data generator
##########################################
sigmagenerator <- function(corst, p, rho = 0){
  if(corst == "ind") { Sigma <- diag(rep(1,p)) }
  if(corst == "cs")  { Sigma <- matrix(rho,p,p); diag(Sigma) <- 1 }
  if(corst == "AR-1") { Sigma <- rho^(abs(outer(1:p, 1:p, "-"))) }
  return(Sigma)
}

datagenerator <- function(n, m, p, B, tempdatadir, type, beta0, beta1, beta2, cpt1, cpt2, 
  intercept, categorical, corst_x, rho_x, corst_y, rho_y, phi_y, seed){
  
  if(length(n)!=1 & length(n)!=B){ stop("n must be a number or a vector of length B.") }
  if(!is.null(seed)){ set.seed(seed) }
  seed.list <- sample(1:1e8, B, replace=FALSE)
  if(length(m)!=1 & length(m)!=n){ stop("m must be a positive integer or a vector of length n.") }
 
  dir.create(tempdatadir)

  for(b in 1:B){
    set.seed(seed.list[b])

    # Generate covariance matrix for covariate matrix
    Sigma_x <- sigmagenerator(corst = corst_x, p = p, rho = rho_x)
  
    # number of visits is the same for every individual, equal to m
    if(length(m)==1){ X <- mvrnorm(n=n*m, mu=rep(0,p), Sigma=Sigma_x); 
                      id <- rep(((b-1)*n+1):(b*n), each=m) }

  
    if(categorical==TRUE){X[,2] <- rep(seq(0.1, 1.0, 0.1),n);}
    if(intercept==TRUE){ X[,1] <- 1; }                  
   
    #create correlations between visits of the same individual
    if(length(m)==1){
      Sigma_y <- phi_y * sigmagenerator(corst = corst_y, p = m, rho = rho_y)
      epsilon <- mvrnorm(n=n, mu=rep(0,m), Sigma=Sigma_y)
      epsilon <- c(t(epsilon))
    }
   
   if(length(m)==n){ stop("function for unbalanced design has not been implemented.") }
  
   if(b==cpt1){
        beta <- beta1
      }else if(b==cpt2){
        beta <- beta2
      }else{
        beta <- beta0
      }

      Xbeta <- drop(as.matrix(X) %*% beta)
  
      if(type=="gaussian"){ y <- Xbeta + epsilon }
      if(type=="binomial"){ 
      intercepts<-0
      y<-rbin(clsize = m, intercepts = intercepts, betas = beta, 
      xformula = ~X, cor.matrix = Sigma_y, link = "logit")$simdata[,1]
      }
   

   if(intercept==TRUE){X<-as.matrix(X[,-1])}
   save(id, y, X, file = paste(tempdatadir, "/", b, ".RData", sep=""))
  }
  

}

