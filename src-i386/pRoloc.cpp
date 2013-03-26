
// includes from the plugin
#include <RcppArmadillo.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP file25fd75026a26( SEXP A, SEXP B, SEXP ind) ;
}

// definition

SEXP file25fd75026a26( SEXP A, SEXP B, SEXP ind ){
BEGIN_RCPP
arma::mat Am = Rcpp::as< arma::mat >(A);
                                    arma::mat Bm = Rcpp::as< arma::mat >(B);
                                    int nbInd = Rcpp::as<int>(ind);

                                    arma::colvec  setOfOnePerturbations(nbInd);
                                    setOfOnePerturbations.zeros();

                                    for (int k=0; k < nbInd; k++) {
                                    setOfOnePerturbations(k) = arma::as_scalar(trans(Am.col(k)) * Bm * Am.col(k));
                                    }

                                    return Rcpp::wrap(setOfOnePerturbations);
                                    
END_RCPP
}



