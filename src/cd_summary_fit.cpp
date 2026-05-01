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


// -- Coordinate descent update for summary statistics --

void cd_summary(
    const arma::vec& XtY, const arma::sp_mat& ld_mat,
    arma::mat& b, int j, arma::vec& a, int l,
    const std::string& penalty, double lam1, double lam2, double gamma,
    arma::vec& df, double& maxChange
) {

  double z = XtY(j) - arma::dot(ld_mat.col(j), b.col(l)) + b(j, l);
  // double z_norm = fabs(z);

  double len;
  if (penalty == "LASSO") len = lasso(z, lam1, lam2, 1);
  if (penalty == "MCP")   len = MCP(z, lam1, lam2, gamma, 1);
  if (penalty == "SCAD")  len = SCAD(z, lam1, lam2, gamma, 1);

  if (len != 0.0 || a(j) != 0.0) {
    // b(j, l) = len * z / z_norm;
    b(j, l) = len;
    double shift = b(j, l) - a(j);
    if (fabs(shift) > maxChange) maxChange = fabs(shift);
    if (len != 0.0) df[l] += len / fabs(z);
  }
}


//' Coordinate descent for BRIERs estimation procedure
//'
//' Coordinate-descent algorithm for penalized regression using summary
//' statistics and an LD reference panel.
//'
//' @param XtY A numeric vector of marginal correlations between genotype and outcome (length p).
//' @param ld_mat A p x p sparse matrix representing the linkage disequilibrium (LD) structure.
//' @param multiplier A numeric vector of per-variable penalty multipliers (length p).
//' @param penalty A string specifying the penalty: "LASSO", "MCP", or "SCAD".
//' @param lam A numeric vector of penalty values, sorted in decreasing order.
//' @param alpha A numeric scalar for the elastic net mixing parameter.
//' @param gamma A numeric scalar for MCP and SCAD penalties.
//' @param eps A numeric scalar for the convergence tolerance.
//' @param max_iter An integer specifying the maximum total iterations allowed.
//' @param dfmax An integer specifying the maximum number of selected variables.
//' @param user A logical flag. If TRUE, lambda is user-specified and iteration starts
//'   from the first value; if FALSE, the first lambda is skipped.
//'
//' @return A list with the following elements:
//' \describe{
//'   \item{beta}{A matrix of coefficients (p x length(lam)).}
//'   \item{iter}{An integer vector of iteration counts for each lambda.}
//'   \item{df}{A numeric vector of degrees of freedom for each lambda.}
//'   \item{dev}{A numeric vector of deviance values for each lambda.}
//' }
// [[Rcpp::export]]
List cd_summary_fit(
    const arma::vec& XtY, const arma::sp_mat& ld_mat,
    const arma::vec& multiplier, const std::string& penalty, const arma::vec& lam,
    double alpha, double gamma, double eps, int max_iter, int dfmax, bool user
) {

  int p = ld_mat.n_cols;
  int L = lam.n_elem;
  int tot_iter = 0;

  // Initialize
  arma::mat b(p, L, arma::fill::zeros);
  arma::ivec iter(L, arma::fill::zeros);
  arma::vec df(L, arma::fill::zeros);
  arma::vec dev(L, arma::fill::zeros);

  // Intermediate
  arma::vec a(p, arma::fill::zeros);
  arma::ivec e(p, arma::fill::zeros);
  int lstart = user ? 0 : 1;
  int nv, violations;
  double maxChange;

  for (int l = lstart; l < L; ++l) {
    R_CheckUserInterrupt();

    if (l != 0) {
      a = b.col(l - 1);
      nv = 0;
      for (int j = 0; j < p; j++) {
        if (a(j) != 0) nv++;
      }
      if (nv > dfmax || tot_iter == max_iter) {
        for (int ll = l; ll < L; ll++) { iter(ll) = NA_INTEGER; }
        break;
      }
    }

    while (tot_iter < max_iter) {
      // active set cycling
      while (tot_iter < max_iter) {
        iter(l) += 1;
        tot_iter++;
        df(l) = 0;
        maxChange = 0.0;

        for (int j = 0; j < p; j++) {
          if (e(j) == 1) {
            double l1 = alpha * lam(l) * multiplier(j);
            double l2 = (1 - alpha) * lam(l) * multiplier(j);
            cd_summary(XtY, ld_mat, b, j, a, l, penalty, l1, l2, gamma, df, maxChange);
          }
        }

        a = b.col(l);
        if (maxChange <= eps) break;
      }

      // KKT check for inactive variables
      violations = 0;
      for (int j = 0; j < p; j++) {
        if (e(j) == 0) {
          double l1 = alpha * lam(l) * multiplier(j);
          double l2 = (1 - alpha) * lam(l) * multiplier(j);
          cd_summary(XtY, ld_mat, b, j, a, l, penalty, l1, l2, gamma, df, maxChange);
          if (b(j, l) != 0) {
            e(j) = 1;
            violations++;
            a(j) = b(j, l);
          }
        }
      }

      if (violations == 0) {
        dev[l] = arma::as_scalar(a.t() * ld_mat * a) - 2.0 * arma::as_scalar(a.t() * XtY);
        break;
      }
    }
  }

  List result = List::create(
    _["beta"] = b,
    _["iter"] = iter,
    _["df"]   = df,
    _["dev"]  = dev
  );
  return result;
}