
### Sequential updating in QIF
renewqif <- function (B, tempdatadir, family, intercept = FALSE, corstr = "independence", QC = FALSE, thred){
  
  if (!(family %in% c("gaussian", "binomial", "poisson"))) { stop("'family' value not supported") }
  if (!corstr %in% c("independence", "exchangeable", "AR-1", "CS+AR1","unstructured")) { stop("'corstr' value not supported") }

  tol <- 1e-4;
  maxit <- 200;
  time_load <- 0;
  Q_stat <- c();
  Q_new <- 0;
  
  p_vec <- c();
 
  N_accum <- 0;
  phi_new <- 0;
  use <- 1;
  cpt <- c()
  L <- 1;


  load(paste(tempdatadir, "/", 1, ".RData", sep = ""));

  ptm <- proc.time()

  beta_new <- coef(glm(y ~ X  ,family = family));
  
  if(intercept==TRUE){X <- cbind(1, X)}
  p <- ncol(X)
  
  ### initial inference statistics ### 
  if(corstr=="independence"){
    g_accum <-  matrix(rep(0, p),nrow = p);
    G_accum <-  matrix(rep(0, p * p), nrow = p);
    C_accum <-  matrix(rep(0, p * p), nrow = p);
  }else if(corstr=="CS+AR1"){
    g_accum <-  matrix(rep(0, 3 * p),nrow = 3 * p);
    G_accum <-  matrix(rep(0, 3 * p * p), nrow = 3 * p);
    C_accum <-  matrix(rep(0, 3 * p * 3 * p), nrow = 3 * p);
  }else{
    g_accum <- matrix(rep(0, 2 * p),nrow = 2 * p);
    G_accum <- matrix(rep(0, 2 * p * p), nrow = 2 * p);
    C_accum <- matrix(rep(0, 2 * p * 2 * p), nrow = 2 * p);
  }


  for(b in 1:B){
    
    beta_old <- beta_new
    g_old <- g_accum
    G_old <- G_accum
    C_old <- C_accum

    load <- proc.time()
    load(paste(tempdatadir,"/",b,".RData",sep=""))
    time_load <- time_load + (proc.time()-load)[3]

    if(intercept == TRUE){X1 <- cbind(1,X)}
    y1 <- y
    obs <- lapply(split(id, id), "length")
    # number of repeated measurements for each subject
    nobs1 <- as.numeric(obs)

    estimate <- increQIF(X1, y1, nobs1, family, corstr, beta_old, g_old, G_old, C_old, maxit, tol)
    
    beta_new <- estimate$beta;
    g_accum <- estimate$g_accum;
    G_accum <- estimate$G_accum;
    C_accum <- estimate$C_accum;
    phi_sub <- estimate$phi_sub;

    N_accum <- N_accum + nrow(X1);
    w <- (N_accum - nrow(X1)) / N_accum;
    phi_new <- phi_new * w + phi_sub * (1-w);   


    if(QC == TRUE && b >= 2){ 
      
      load(paste(tempdatadir,"/", L,".RData",sep = ""))

      if(intercept == TRUE){X0 <- cbind(1,X)}
      y0 <- y
      obs <- lapply(split(id, id), "length")
      nobs0 <- as.numeric(obs)

      Xc <- rbind(X0, X1)
      yc <- c(y0, y1)
      nobsc <- c(nobs0, nobs1)

      beta2 <- increQIF_sub(Xc, yc, nobsc, family, corstr, beta_old, maxit, tol)$beta_sub
    
      r_H0 <- increQIF_test(X0, y0, nobs0, family, corstr, beta2) 
      r_new <- increQIF_test(X1, y1, nobs1, family, corstr, beta2)

      Q_old <- r_H0$Q_H0
      Q_new <- r_new$Q_H0

      if(b == 2){
        Cb_old <- r_H0$Cb_sub;
        Cb_new <- r_new$Cb_sub;
        weight_m <- as.matrix(bdiag(Cb_old, Cb_new))
        df <-  as.numeric(rankMatrix(weight_m) - p)
      }

      Q_stat <- c(Q_stat, Q_old + Q_new);       
      pvalue <- pchisq(Q_old + Q_new, df = df, lower.tail = FALSE);   
      p_vec <- c(p_vec, pvalue)  # vector of unadjusted p-values
      
            if(p_vec[length(p_vec)] <= thred){
              use <- 0  
              cpt <- c(cpt, b)  # record the change point
            }else{
              use <- 1
              cpt <- c(cpt, 0)
              L <- b;  # update L only when the current data batch passes the test
            }

        if (use == 0) { # if change point is not detected       
          # update with data batch b      
          beta_new <- beta_old
          g_accum <- g_old
          G_accum <- G_old
          C_accum <- C_old

          N_accum <- N_accum - nrow(X1)

        }
          
    }
     
  }
  
  J_accum <- t(G_accum) %*% ginv(C_accum) %*% G_accum;
  varb <- sqrt(diag(ginv(J_accum)));

  out_beta <- cbind(Estimate = drop(beta_new), StdErr = varb, Zscore = beta_new / varb)

  
  time_total <- (proc.time()-ptm)[3]
  time_run <- time_total-time_load 
  out <- list()

  if(QC==TRUE){
      out$cpt <- cpt
      out$Q <- Q_stat
      out$N <- N_accum
  } 
  out$beta <- out_beta
  out$time <- time_total
  out$run <- time_run

  return(out)
}