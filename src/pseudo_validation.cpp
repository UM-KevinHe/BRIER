// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <new>
#include <iterator>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
#include "utils.h"

using namespace Rcpp;
using namespace std;
using namespace arma;


//' Pseudo-validation score for summary-statistics models
//'
//' Compute a pseudo-validation score for each column of a coefficient matrix,
//' using a sparse LD reference (XtX) and a marginal correlation vector (r).
//' The score is defined as beta'r / sqrt(beta' XtX beta), and is used for
//' tuning lambda when validation data is unavailable.
//'
//' @param beta A p x L numeric matrix of coefficients across L lambda values.
//' @param XtX A p x p sparse LD reference matrix.
//' @param r A numeric vector (length p) of marginal correlations between
//'   genotype and outcome.
//'
//' @return A numeric vector of length L containing the pseudo-validation
//'   score for each lambda. Non-finite values (e.g. when beta is all zeros)
//'   are returned as NA.
// [[Rcpp::export]]
arma::vec pseudo_validation(
  const arma::mat& beta,
  const arma::sp_mat& XtX,
  const arma::vec& r
) {
  arma::vec numer = beta.t() * r;
  arma::mat XtXbeta = XtX * beta;
  arma::vec denom = arma::sqrt(arma::sum(beta % XtXbeta, 0).t());

  arma::vec out = numer / denom;
  out.elem(arma::find_nonfinite(out)).fill(NA_REAL);
  return out;
}