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


//' Standardize a design matrix
//'
//' Center and scale each column of X to have mean 0 and standard deviation 1.
//' Constant columns (sd = 0) are centered but not scaled.
//'
//' @param Xr An n x p numeric matrix.
//'
//' @return A list with the following elements:
//' \describe{
//'   \item{standardized}{The standardized n x p matrix.}
//'   \item{center}{A numeric vector of length p of column means.}
//'   \item{scale}{A numeric vector of length p of column standard deviations.}
//' }
//' @export
// [[Rcpp::export]]
List standardize_X(const NumericMatrix& Xr) {
  arma::mat X = as<arma::mat>(Xr);

  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat XX(n, p, arma::fill::zeros);
  NumericVector center(p);
  NumericVector scale(p);

  for (int j = 0; j < p; ++j) {
    double mean_j = arma::mean(X.col(j));
    arma::vec centered = X.col(j) - mean_j;
    double sd_j = std::sqrt(arma::mean(arma::square(centered)));

    XX.col(j) = centered;
    if (sd_j > 0) {
      XX.col(j) /= sd_j;
    }

    center[j] = mean_j;
    scale[j] = sd_j;
  }

  NumericMatrix XXr = wrap(XX);

  SEXP dimnames = Xr.attr("dimnames");
  if (!Rf_isNull(dimnames)) {
    XXr.attr("dimnames") = dimnames;

    List dn(dimnames);
    if (dn.size() == 2 && !Rf_isNull(dn[1])) {
      center.attr("names") = dn[1];
      scale.attr("names") = dn[1];
    }
  }

  return List::create(
    _["standardized"] = XXr,
    _["center"] = center,
    _["scale"] = scale
  );
}

//' Weighted standardization of a design matrix
//'
//' Center and scale each column of X using observation weights. Weights are
//' rescaled internally to sum to n for numerical stability.
//'
//' @param X An n x p numeric matrix.
//' @param wt A numeric vector of observation weights (length n).
//'
//' @return A list with the following elements:
//' \describe{
//'   \item{standardized}{The standardized n x p matrix.}
//'   \item{center}{A numeric vector of length p of weighted column means.}
//'   \item{scale}{A numeric vector of length p of weighted column standard deviations.}
//' }
//' @export
// [[Rcpp::export]]
List wstandardize_X(const arma::mat &X, const arma::vec &wt) {

  // Dimensions
  int n = X.n_rows;
  int p = X.n_cols;

  // Pre-allocate matrices and vectors
  arma::mat XX = arma::mat(n, p, fill::zeros);
  arma::vec c = arma::vec(p, fill::zeros);
  arma::vec s = arma::vec(p, fill::zeros);

  // rescale weights to avoid very small numbers
  arma::vec w = wt * n;
  double wsum = arma::sum(w);

  for (int j = 0; j < p; ++j) {

    // Weighted center
    double mean_j = arma::dot(w, X.col(j)) / wsum;
    c(j) = mean_j;
    arma::vec col_centered = X.col(j) - mean_j;

    // Weighted scale
    double var_j = arma::dot(w, col_centered % col_centered) / wsum;
    double sd_j = std::sqrt(var_j);

    s(j) = sd_j;
    if (sd_j > 0) {
      XX.col(j) = col_centered / sd_j;
    } 
  }

  List result = List::create(_["standardized"] = XX, _["center"] = c, _["scale"] = s);
  return result;
}


//' Unstandardize coefficients
//'
//' Given coefficients fitted on a standardized design matrix, recover the
//' coefficients on the original scale by reversing the centering and scaling.
//'
//' @param b0 A 1 x L matrix of intercepts (one per lambda).
//' @param b A p x L matrix of coefficients (one column per lambda).
//' @param scale A numeric vector of length p of column standard deviations.
//' @param center A numeric vector of length p of column means.
//' @param yy_center A numeric scalar (the response center, used for Gaussian models).
//' @param L An integer specifying the number of lambda values.
//' @param p An integer specifying the number of predictors.
//'
//' @return A (p+1) x L matrix where the first row is the unstandardized
//'   intercept and the remaining rows are the unstandardized coefficients.
// [[Rcpp::export]]
arma::mat unstand_beta(
  arma::mat &b0, arma::mat &b, arma::vec &scale, arma::vec &center, 
  double yy_center, int L, int p
) {

  arma::mat bb(p + 1, L, arma::fill::zeros);
  double sum_tmp;

  for (int l = 0; l < L; l++) {
    sum_tmp = 0;
    for (int j = 0; j < p; j++) {
      bb(j + 1, l) = b(j, l) / scale(j);
      sum_tmp += bb(j + 1, l) * center(j);
    }
    bb(0, l) = b0(l) + yy_center - sum_tmp;
  }

  return bb;
}