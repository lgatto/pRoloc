
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// user includes

// declarations
// RcppExport is an alias to ‘extern ”C”‘ defined by Rcpp.
RcppExport SEXP C_setOfOnePerturbation(SEXP A, SEXP B, SEXP ind) ;


// definition

SEXP C_setOfOnePerturbation(SEXP A, SEXP B, SEXP ind){
  arma::mat Am = Rcpp::as< arma::mat >(A);
  arma::mat Bm = Rcpp::as< arma::mat >(B);
  int nbInd = Rcpp::as<int>(ind);
  
  arma::colvec  setOfOnePerturbations(nbInd);
  setOfOnePerturbations.zeros();
  
  for (int k=0; k < nbInd; k++) {
    setOfOnePerturbations(k) = arma::as_scalar(trans(Am.col(k)) * Bm * Am.col(k));
  }
 
  return Rcpp::wrap(setOfOnePerturbations);
}



