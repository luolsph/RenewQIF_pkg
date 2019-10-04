// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

//[[Rcpp::export]]
List increQIF(arma::mat X, arma::vec y, arma::vec nobs, String family, String corstr,
	arma::vec beta_old, arma::vec g_accum, arma::mat G_accum,
	arma::mat C_accum, int maxit, double tol){
 
    int niter = 0;
    
    bool stop_flag = FALSE;  
    bool converged = FALSE;
    
    int N = X.n_rows;
    int p = X.n_cols;
    int n = nobs.n_elem;
   
    arma::vec g_sum; 
    arma::mat G_temp;
    arma::mat C_temp;

    arma::vec gb_new;
    arma::mat Gb_new;
    arma::mat Cb_new;

   
    //initialization by the beta estimated from previous data
    arma::vec beta_new = beta_old;

    arma::vec mu;
    arma::vec vu;
    arma::mat M0;
    arma::mat M1;
    arma::mat M2;
    arma::mat qif1dev;
    arma::mat qif2dev;

    while (!stop_flag){
    	niter += 1;
        //reset gb to 0 after each iteration
        if(corstr=="independence"){
        gb_new = zeros<vec>(p );
        Gb_new = zeros<mat>(p , p);
        Cb_new = zeros<mat>(p , p);
        }else if(corstr=="CS+AR1"){
        gb_new = zeros<vec>(p * 3);
        Gb_new = zeros<mat>(p * 3, p);
        Cb_new = zeros<mat>(p * 3, p * 3);
        }else{
        gb_new = zeros<vec>(p * 2);
        Gb_new = zeros<mat>(p * 2, p);
        Cb_new = zeros<mat>(p * 2, p * 2);
        }
        
        // update gb_new with beta_new over iterations
        arma::vec eta = X * beta_new;

        if(family == "gaussian") {
            mu = eta ;
            vu = ones<vec>(N) ;
        } else if(family == "binomial") {
            mu = exp(eta)/( 1 + exp(eta) ) ;
            vu = exp(eta)/(pow( (1 + exp(eta)), 2 )) ;
        } else if(family == "poisson") {
            mu = exp(eta) ;
            vu = exp(eta) ;
        } else{Rcpp::stop("Unknown distribution family\n");}
    
    
    int loc1 = 0;
    int loc2 = -1;
    for (int i = 0; i < n; i++){
        loc1 = loc2 + 1 ;
        loc2 = loc1 + nobs[i] -1;

        arma::mat Xi = X.rows(loc1, loc2);
        arma::vec yi = y.subvec(loc1, loc2);
        arma::vec mui = mu.subvec(loc1, loc2);
        arma::vec vui = vu.subvec(loc1, loc2);

        int mi = nobs[i] ;

        arma::mat Ai_half = diagmat(sqrt(vui)) ;
        arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;

        if (corstr == "independence") {
        	M0 = eye<mat>(mi, mi);
        	
        } else if (corstr == "exchangeable") {
        	M0 = eye<mat>(mi, mi);
        	// M1 is a matrix with 1 on off-diagonal elements
        	M1 = ones(mi, mi) - M0;  
        } else if (corstr == "AR-1") {
        	M0 = eye<mat>(mi, mi);
        	M1 = zeros<mat>(mi, mi);
        	for (int j = 0; j < mi; j++ ){
        		for (int k = 0; k < mi; k++){
        			if(abs(j-k)==1){ M1(j, k) = 1; }
        		}
        	}        	
        } else if (corstr == "unstructured"){
            int m = N/n;
            M0 = eye<mat>(m, m);
            M1 = zeros<mat>(m, m);
            int loc3 = 0;
            int loc4 = -1;
            for (int i = 0; i < n; i++){
                loc3 = loc4 +1;
                loc4 = loc3 + nobs[i] -1;
                arma::mat Xi = X.rows(loc3, loc4);
                arma::vec yi = y.subvec(loc3, loc4);
                arma::vec mui = mu.subvec(loc3, loc4);
                M1 += (yi - mui) * (yi - mui).t(); 
            }
            M1 = M1 / n;
        }else if(corstr=="CS+AR1"){
            M0 = eye<mat>(mi, mi);
            M1 = ones(mi, mi) - M0;  
            M2 = zeros<mat>(mi, mi);
            for (int j = 0; j < mi; j++ ){
                for (int k = 0; k < mi; k++){
                    if(abs(j-k)==1){ M2(j, k) = 1; }
                }
            }   
        }else {Rcpp::stop("Unknown correlation structure\n");}   
        
        if(corstr=="independence"){
        vec gi_new = Xi.t() * (yi - mui);
        gb_new += gi_new;
        Gb_new += Xi.t() * Ai_half * M0 * Ai_half * Xi ;
        Cb_new += gi_new * gi_new.t();
        }else if (corstr=="CS+AR1"){
        vec gi_new = join_cols( join_cols( Xi.t() * (yi - mui), 
            Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui)),
            Xi.t() * Ai_half * M2 * Ai_inv_half * (yi-mui) );
        gb_new += gi_new;
        Gb_new += join_cols( join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
            Xi.t() * Ai_half * M1 * Ai_half * Xi),
            Xi.t() * Ai_half * M2 * Ai_half * Xi) ;
        Cb_new += gi_new * gi_new.t();

        }else{
        vec gi_new = join_cols(Xi.t() * (yi - mui), Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
        gb_new += gi_new;
        Gb_new += join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
            Xi.t() * Ai_half * M1 * Ai_half * Xi) ;
        Cb_new += gi_new * gi_new.t();
        }

    }       

        g_sum = g_accum + G_accum * (beta_old - beta_new) + gb_new;

        G_temp = G_accum + Gb_new;
        C_temp = C_accum + Cb_new;

        qif1dev = G_temp.t() * pinv(C_temp) * g_sum ;

        qif2dev = G_temp.t() * pinv(C_temp) * G_temp ;
        

        vec d_beta = solve(qif2dev, qif1dev);

        double df_beta = as_scalar(qif1dev.t() * d_beta);

        beta_new += d_beta;

        if(fabs(df_beta) < tol) {converged = TRUE; stop_flag = TRUE;}
        if (niter > maxit) {stop_flag = TRUE;}
    }
//    if (converged==FALSE) {Rcpp::stop("algorithm reached 'maxit' but did not converge\n"); }
    
    vec res = (y - mu)/sqrt(vu);

    mat J_accum = G_accum.t() * pinv(C_accum) * G_accum;

    double phi_sub = as_scalar(res.t() * res / ( N - p ) );

    double Q = as_scalar(g_sum.t() * pinv(C_temp) * g_sum);

    double Q_adjust_old = as_scalar((g_sum - gb_new).t() * pinv(C_accum) * (g_sum - gb_new));

    double Q_new = as_scalar(gb_new.t() * pinv(Cb_new) * gb_new);

    double Q_old = as_scalar(g_accum.t() * pinv(C_accum) * g_accum);

   // double Q_test = Q_old + Q_sub - Q_accum_old;

return List::create(Named("beta") = beta_new,
                    Named("g_accum") = g_sum, 
                    Named("G_accum") = G_temp, 
                    Named("C_accum") = C_temp,
                    Named("gb_new") = gb_new, 
                    Named("Cb_new") = Cb_new,
                    Named("J_accum") = J_accum,
                    Named("Q_value") = Q,
                    Named("Q_old") = Q_old, 
                    Named("Q_adjust_old") = Q_adjust_old, 
                    Named("Q_new") = Q_new,
                    Named("phi_sub") = phi_sub          
                    );

}