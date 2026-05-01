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


// -- Coordinate descent update for weighted Gaussian regression --

void cd_wgaussian(
  arma::mat &b, const arma::mat &X, arma::vec &r, const arma::vec &wt, int j, 
  const arma::vec &a, int n_obs, int l, 
  const std::string &penalty, double lam1, double lam2, double gamma,
  arma::vec &df, double &maxChange
  ){
  // n_obs: number of observations (n_obs = N1(1) given N1(0) = 0)
  // X: (n_obs)*p matrix; r: (n_obs)*1 vector
  // b: beta estimate, p*L matrix; a: old beta, p vector 

  // Compute U 
  double z = 0;
  z = wcrossprod(X, r, wt, n_obs, j);   
  z += a(j);

  double len = 0;
  if (penalty == "LASSO") len = lasso(z, lam1, lam2, 1);
  else if (penalty == "SCAD")  len = SCAD(z, lam1, lam2, gamma, 1); 
  else if (penalty == "MCP")   len = MCP(z, lam1, lam2, gamma, 1);
  
  if (len != 0 || a(j) != 0) {
    b(j, l) = len;
    double shift = b(j, l) - a(j);
    if (fabs(shift) > maxChange) {
        maxChange = fabs(shift);
    }
    r -= X.col(j) * shift;
  }

  // Update degrees of freedom (df)
  if (len != 0) {
    if (z != 0) df(l) += fabs(len / z);
  }
}


// -- Sequential strong rule screening for weighted Gaussian --

void ssr_wgaussian(
  arma::ivec &e2, const arma::mat &X, const arma::vec &r, const arma::vec &wt, int j, int n_obs, int l, 
  const std::string &penalty, const arma::vec &lam, double lam_max, double alpha, double gamma, const arma::vec &multiplier
  ) {
  double cutoff;
  double TOLERANCE = 1e-8;
  
  double z = 0;
  z =  wcrossprod(X, r, wt, n_obs, j);   

  if (l != 0) {
      if (penalty == "LASSO") cutoff = 2 * lam(l) - lam(l-1);
      if (penalty == "SCAD")  cutoff = lam(l) + gamma/(gamma-2)*(lam(l) - lam(l-1));
      if (penalty == "MCP")   cutoff = lam(l) + gamma/(gamma-1)*(lam(l) - lam(l-1));
  } else {
      if (penalty == "LASSO") cutoff = 2 * lam(l) - lam_max;
      if (penalty == "SCAD")  cutoff = lam(l) + gamma/(gamma-2)*(lam(l) - lam_max);
      if (penalty == "MCP")   cutoff = lam(l) + gamma/(gamma-1)*(lam(l) - lam_max);
  }

  if (fabs(z) + TOLERANCE > (cutoff * alpha * multiplier(j))) {
      e2(j) = 1; // not reject, in strong set
  } else {
      e2(j) = 0; // reject
  }
}


//' Coordinate descent for weighted Gaussian regression with sequential strong rule
//'
//' Coordinate-descent algorithm for penalised weighted Gaussian regression with
//' sequential strong rule (SSR) screening over a sequence of penalty values.
//'
//' @param X A numeric matrix of predictors (n x p).
//' @param y A numeric response vector of length n.
//' @param wt A numeric vector of observation weights (length n).
//' @param multiplier A numeric vector of per-variable penalty multipliers (length p).
//'   Variables with multiplier = 0 are treated as unpenalised.
//' @param penalty A string specifying the penalty: "LASSO", "MCP", or "SCAD".
//' @param lam A numeric vector of penalty values, sorted in decreasing order.
//' @param alpha A numeric scalar for the elastic net mixing parameter.
//' @param gamma A numeric scalar for MCP and SCAD penalties.
//' @param lam_max A numeric scalar for the maximum lambda value (used in SSR cutoff).
//' @param max_iter An integer specifying the maximum total iterations allowed.
//' @param eps A numeric scalar for the convergence tolerance.
//' @param dfmax An integer specifying the maximum number of selected variables.
//'
//' @return A list with the following elements:
//' \describe{
//'   \item{beta0}{A numeric vector of intercepts for each lambda.}
//'   \item{beta}{A matrix of coefficients (p x length(lam)).}
//'   \item{dev}{A numeric vector of weighted residual sums of squares.}
//'   \item{eeta}{A matrix of fitted linear predictors (n x length(lam)).}
//'   \item{df}{A numeric vector of degrees of freedom for each lambda.}
//'   \item{iter}{An integer vector of iteration counts for each lambda.}
//' }
// [[Rcpp::export]]
List cd_wgaussian_fit_ssr(
    const arma::mat &X, const arma::vec &y, const arma::vec &wt, const arma::vec &multiplier, 
    const std::string &penalty, const arma::vec &lam, double alpha, double gamma, 
    double lam_max, int max_iter, double eps, int dfmax
  ){

  // n_obs: number of observations
  // X: Data matrix: n_obs x p; y, r: n_obs*1 vector
  // b: beta estimate, p x L matrix; a: old beta, p vector
  
  int n_obs = y.n_elem;
  int L = lam.n_elem;
  int p = X.n_cols;
  int tot_iter = 0;

  // b0: intercept vector (length L); b: coefficient matrix (dimensions p x L)
  arma::vec b0(L, fill::zeros);
  arma::mat b(p, L, fill::zeros);
  arma::vec loss(L, fill::zeros);
  arma::mat Eta(n_obs, L, fill::zeros);
  arma::vec df(L, fill::zeros);
  arma::ivec iter(L, fill::zeros);

  double a0 = 0;            // Old intercept 
  arma::vec a(p, fill::zeros);  // Old beta coefficients

  arma::vec r = y;          // Residuals, initially equal to y

  int lstart = 0, nv, violations;
  double shift, l1, l2, maxChange;
  
  // Initialize screening variables
  arma::ivec e1(p, fill::zeros);  // Ever-active set indicator 
  arma::ivec e2(p, fill::zeros);  // strong set indicator

  double wrss = arma::dot(wt, r % r);
  // if (lam(0) == lam_max) {
  //   loss(0) = rss;
  //   lstart = 1; // Start from the second lambda value
  // }
  lstart = 0;
  double wsdy = sqrt(wrss);

  for (int l = lstart; l < L; l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      // Warm start: update old beta from previous lambda iteration
      a0 = b0(l-1);
      for (int j = 0; j < p; j++) {
        a(j) = b(j, l-1);
      }
      
      // Check dfmax conditions :
      nv = 0;
      for (int j = 0; j < p; j++) {
      // Check the first coefficient; adjust if needed.
      if ( a(j) != 0 ) {
        nv++;
      }
      }
      if (nv > dfmax || tot_iter == max_iter) {
      if (tot_iter == max_iter) {
        Rcpp::Rcout << "Algorithm has reached the maximum number of total iterations, stops..." << endl;
      } else {
        Rcpp::Rcout << "Algorithm has selected the maximum number of variables, stops..." << endl;
      }
      for (int ll = l; ll < L; ll++) {
        iter(ll) =  NA_INTEGER;
      }
      break;
      }
    }

    // Set SSR rule
    for (int j = 0; j < p; j++) {
      ssr_wgaussian(e2, X, r, wt, j, n_obs, l, penalty, lam, lam_max, alpha, gamma, multiplier);
    }

    while (tot_iter < max_iter) { // outer loop for coordinate descent algorithm
      while (tot_iter < max_iter) { // strong set loop
        while (tot_iter < max_iter) { // convergence loop
          tot_iter++;
          iter(l) += 1;
          df(l) = 0;
          maxChange = 0;
  
          // ---- Update Intercept ----
          shift = arma::dot(wt, r);        
          if (fabs(shift) > maxChange) {
            maxChange = fabs(shift);
          }
          b0(l) = shift + a0;
          for (int i = 0; i < n_obs; i++){
            r(i) -= shift;
          }
  
          // ---- Update Covariates: if multiplier == 0, unpenalized ----
          for (int j = 0; j < p; j++) {
            l1 = alpha * lam(l) * multiplier(j);
            l2 = (1 - alpha) * lam(l) * multiplier(j);
            if (e1(j) == 1) {
              cd_wgaussian(b, X, r, wt, j, a, n_obs, l, penalty, l1, l2, gamma, df, maxChange);
            }
          }
  
          // ---- Update Old Beta ----
          a0 = b0(l);
          for (int j = 0; j < p; j++) {
            a(j) = b(j, l);
          }
          
          // Check convergence
          if (maxChange <= eps * wsdy) {
              break;
          }
        }

        // ---- Check KKT Conditions for in strong set but not active set ----
        violations = 0;
        // double betachange;
        for (int j = 0; j < p; j++) {
          if (e1(j) == 0 && e2(j) == 1) {
            l1 = alpha * lam(l) * multiplier(j);
            l2 = (1 - alpha) * lam(l) * multiplier(j);
            double z = wcrossprod(X, r, wt, n_obs, j);
            if (std::fabs(z) > l1) {
              e1(j) = 1;
              e2(j) = 1;
              violations++;
            }
          }
        }
        if (violations == 0) break;
      }

      // ---- Check KKT Conditions not in strong set ----
      violations = 0;
      // double betachange;
      for (int j = 0; j < p; j++) {
        if (e2(j) == 0) {
          l1 = alpha * lam(l) * multiplier(j);
          l2 = (1 - alpha) * lam(l) * multiplier(j);
          double z = wcrossprod(X, r, wt, n_obs, j);
          if (std::fabs(z) > l1) {
            e1(j) = 1;
            e2(j) = 1;
            violations++;
          }
        }
      }

      if (violations == 0) {
        loss(l) = arma::dot(wt, r % r);
        Eta.col(l) = y - r;
        break;
      }
    } // end outer while
  } // end for (l in lambda sequence)

  return Rcpp::List::create(
    _["beta0"] = b0, 
    _["beta"] = b, 
    _["dev"] = loss, 
    _["eeta"] = Eta, 
    _["df"] = df, 
    _["iter"] = iter
  );
}