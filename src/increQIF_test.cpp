// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

//[[Rcpp::export]]
List increQIF_test(arma::mat X, arma::vec y, arma::vec nobs, String family, String corstr,
	arma::vec beta_old){

	//int niter2= 0;

    //bool stop_flag2= FALSE;
    //bool converged2= FALSE;
    
    int N = X.n_rows;
    int p = X.n_cols;
    int n = nobs.n_elem;

    arma::vec gb_sub;
//    arma::mat Gb_sub;
    arma::mat Cb_sub;
   
    //initialization by the beta estimated from previous data

    //arma::vec beta_sub = beta_old;

    arma::vec mu;
    arma::vec vu;
    arma::mat M0;
    arma::mat M1;
    arma::mat M2;


        if(corstr=="independence"){
        gb_sub = zeros<vec>(p );
//        Gb_sub = zeros<mat>(p , p);
        Cb_sub = zeros<mat>(p , p);
        }else if(corstr=="CS+AR1"){
        gb_sub = zeros<vec>(p * 3);
//        Gb_sub = zeros<mat>(p * 3, p);
        Cb_sub = zeros<mat>(p * 3, p * 3);
        }else{
        gb_sub = zeros<vec>(p * 2);
//        Gb_sub = zeros<mat>(p * 2, p);
        Cb_sub = zeros<mat>(p * 2, p * 2);
        }
        
        // update gb_new with beta_new over iterations
        arma::vec eta2 = X * beta_old;

        if(family == "gaussian") {
            mu = eta2 ;
            vu = ones<vec>(N) ;
        } else if(family == "binomial") {
            mu = exp(eta2)/( 1 + exp(eta2) ) ;
            vu = exp(eta2)/(pow( (1 + exp(eta2)), 2 )) ;
        } else if(family == "poisson") {
            mu = exp(eta2) ;
            vu = exp(eta2) ;
        } else{Rcpp::stop("Unknown distribution family\n");}
    
    int loc1 = 0;
    int loc2 = -1;
    for (int i = 0; i < n; i++){
        loc1 = loc2 + 1 ;
        loc2 = loc1 + nobs[i] -1;

        arma::mat Xi = X.rows(loc1, loc2) ;
        arma::vec yi = y.subvec(loc1, loc2);
        arma::vec mui = mu.subvec(loc1, loc2);
        arma::vec vui = vu.subvec(loc1, loc2) ;

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
        vec gi_sub = Xi.t() * (yi - mui);
        gb_sub += gi_sub;
    //    Gb_sub += Xi.t() * Ai_half * M0 * Ai_half * Xi ;
        Cb_sub += gi_sub * gi_sub.t();

        }else if(corstr == "CS+AR1"){
        vec gi_sub = join_cols( join_cols( Xi.t() * (yi - mui), 
            Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui)),
            Xi.t() * Ai_half * M2 * Ai_inv_half * (yi-mui) );
        gb_sub += gi_sub;
    //    Gb_sub += join_cols( join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
    //        Xi.t() * Ai_half * M1 * Ai_half * Xi),
    //        Xi.t() * Ai_half * M2 * Ai_half * Xi) ;
        Cb_sub += gi_sub * gi_sub.t();
        
        }else{
        arma::vec gi_sub = join_cols(Xi.t() * (yi - mui), Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
        gb_sub += gi_sub;
    //    Gb_sub += join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
    //        Xi.t() * Ai_half * M1 * Ai_half * Xi) ;
        Cb_sub += gi_sub * gi_sub.t();
        }

    }       

    double Q_H0 = as_scalar(gb_sub.t() * pinv(Cb_sub) * gb_sub);

    return List::create(
                        Named("gb_sub") = gb_sub,
                    //    Named("Gb_sub") = Gb_sub,
                        Named("Cb_sub") = Cb_sub,
                        Named("Q_H0") = Q_H0
                      //  Named("phi_sub") = phi_sub      
                    );

}