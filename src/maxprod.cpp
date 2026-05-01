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


//' Maximum scaled inner product
//'
//' Compute the maximum scaled |X'r| / m(j) across all variables, treating
//' variables with m(j) == 0 as excluded from the maximum. Used for selecting
//' the maximum lambda in penalized regression with individual-level data.
//'
//' @param X An n x p numeric matrix.
//' @param r A numeric vector of length n.
//' @param n An integer specifying the number of rows.
//' @param p An integer specifying the number of columns.
//' @param m A numeric vector of per-variable scaling factors (length p).
//'   Variables with m(j) == 0 are skipped.
//'
//' @return A numeric scalar: the maximum value of |X'r|(j) / m(j).
// [[Rcpp::export]]
double maxprod(
  const arma::mat &X, const arma::vec &r, 
  int n, int p, const arma::vec &m
) {
  double max_val = 0;
  double z = 0;

  for (int j = 0; j < p; j++) {
    if (m(j) != 0) {
      double U = 0.0;
      for (int i = 0; i < n; i++) {
        U += X(i, j) * r(i);
      }
      z = std::abs(U) / m(j);
      // Update the maximum value
      if (z > max_val) {
        max_val = z;
      }
    }
  }

  return(max_val);
}

//' Maximum Lambda Selection for Summary Statistics
//'
//' Calculate the maximum lambda for LASSO when each SNP is its own group,
//' treating any SNP with m(j) == 0 as unpenalized. Iteratively fits the
//' unpenalized variables, then computes the marginal score for each
//' penalized variable.
//'
//' @param XtY A numeric vector (length p) of marginal correlations between
//'   genotype and outcome.
//' @param ld_mat A p x p sparse matrix representing the LD structure.
//' @param m A numeric vector of per-SNP scaling factors (length p).
//'   m(j) == 0 indicates an unpenalized variable.
//' @param alpha A numeric scalar for the elastic net mixing parameter.
//' @param eps A numeric scalar specifying the convergence threshold for the
//'   unpenalized update.
//' @param max_iter An integer specifying the maximum iterations for the
//'   unpenalized update.
//'
//' @return A numeric scalar: the maximum lambda value among penalized SNPs.
// [[Rcpp::export]]
double maxlambda_summary(
  const arma::vec& XtY, const arma::sp_mat& ld_mat,
  const arma::vec& m, double alpha, double eps, int max_iter
) {
  int p = ld_mat.n_cols;

  // Identify unpenalized indices: m[i] == 0
  std::vector<int> unpenalized;
  for (int i = 0; i < p; i++) {
    if (m(i) == 0) unpenalized.push_back(i);
  }

  // Iteratively update unpenalized part
  arma::vec a = arma::zeros<arma::vec>(p);
  arma::vec b = arma::zeros<arma::vec>(p);
  for (int iter = 0; iter < max_iter; ++iter) {
    double maxChange = 0.0;
    for (int j : unpenalized) {
      double shift = XtY[j] - arma::dot(ld_mat.col(j), b);
      maxChange = std::max(maxChange, std::abs(shift));
      b(j) = a(j) + shift;
    }
    a = b;
    if (maxChange <= eps) break;
  }

  // Compute z for all SNPs
  arma::vec z = XtY - (ld_mat * b) + b;

  // Compute lambda_m for each penalized SNP
  arma::vec lambda_m(p, arma::fill::zeros);
  for (int j = 0; j < p; j++) {
    if (m(j) > 0) {
      lambda_m(j) = std::abs(z(j)) / (alpha * m(j));
    }
  }

  return lambda_m.max();
}