// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>


using namespace Rcpp; using namespace arma;

// [[Rcpp::export]]

arma::mat inverse(arma::mat x) {
  return( inv(x) ) ;
}


// [[Rcpp::export]]
sp_mat multiply(sp_mat A, sp_mat B){ return A * B;
}
