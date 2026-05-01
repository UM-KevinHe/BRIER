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

// -- Coordinate descent update for unweighted GLM --

void cd_glm(
  arma::mat &b, const arma::mat &X, arma::vec &r, const arma::vec &w, 
  arma::vec &eta, int j, const arma::vec &a, int n_obs, int l,  
  const std::string &penalty, double lam1, double lam2, double gamma,
  arma::vec &df, double &maxChange
) {

  double v = wsqsum(X, w, n_obs, j) / n_obs; 
  double u = wcrossprod(X, r, w, n_obs, j) / n_obs + v * a(j);
  
  double len = 0.0;
  if (penalty == "LASSO")  len = lasso(u, lam1, lam2, v);
  else if (penalty == "SCAD") len = SCAD(u, lam1, lam2, gamma, v); 
  else if (penalty == "MCP")  len = MCP(u, lam1, lam2, gamma, v);
  
  if (len != 0 || a(j) != 0) {
    b(j, l) = len;
    double shift = b(j, l) - a(j);
    // Update residual and eta:
    if (shift != 0) {
      r -= X.col(j) * shift;
      eta += X.col(j) * shift;
      if (fabs(shift) * sqrt(v) > maxChange) maxChange = fabs(shift) * sqrt(v);
    }
  }
  
  if (len != 0) df(l) += 1.0; 
}


// -- Sequential strong rule screening for GLM --

void ssr_glm(
  arma::ivec &e2, const arma::mat &X, const arma::vec &z, int j, int n_obs, int l, 
  const std::string &penalty, const arma::vec &lam, double lmax,
  double alpha, double gamma, const arma::vec &multiplier
){
  double cutoff;
  double TOLERANCE = 1e-8;

  if (l != 0) {
    if (penalty == "LASSO") cutoff = 2 * lam(l) - lam(l-1);
    if (penalty == "SCAD")  cutoff = lam(l) + gamma / (gamma-2) * (lam(l) - lam(l-1));
    if (penalty == "MCP")   cutoff = lam(l) + gamma / (gamma-1) * (lam(l) - lam(l-1));
  } else {
    if (penalty == "LASSO") cutoff = 2 * lam(l) - lmax;
    if (penalty == "SCAD")  cutoff = lam(l) + gamma / (gamma - 2) * (lam(l) - lmax);
    if (penalty == "MCP")   cutoff = lam(l) + gamma / (gamma - 1) * (lam(l) - lmax);
  }

  if (fabs(z(j)) + TOLERANCE > (cutoff * alpha * multiplier(j))) {
    e2(j) = 1; // not reject, in strong set
  } 
}


//' Coordinate descent for GLM with sequential strong rule
//'
//' Coordinate-descent algorithm for penalised generalised linear models (binomial
//' or Poisson) with sequential strong rule (SSR) screening over a sequence of
//' penalty values.
//'
//' @param X A numeric matrix of predictors (n x p).
//' @param y A numeric response vector of length n.
//' @param family A string: "binomial" or "poisson".
//' @param multiplier A numeric vector of per-variable penalty multipliers (length p).
//'   Variables with multiplier = 0 are treated as unpenalised.
//' @param penalty A string specifying the penalty: "LASSO", "MCP", or "SCAD".
//' @param lambda A numeric vector of penalty values, sorted in decreasing order.
//' @param alpha A numeric scalar for the elastic net mixing parameter.
//' @param gamma A numeric scalar for MCP and SCAD penalties.
//' @param lam_max A numeric scalar for the maximum lambda value (used in SSR cutoff).
//' @param max_iter An integer specifying the maximum total iterations allowed.
//' @param dfmax An integer specifying the maximum number of selected variables.
//' @param eps A numeric scalar for the convergence tolerance.
//' @param user A logical flag. If TRUE, lambda is user-specified and iteration starts
//'   from the first value; if FALSE, the first lambda is treated as the null model.
//' @param warn A logical flag. If TRUE, issues a warning when the model saturates.
//'
//' @return A list with the following elements:
//' \describe{
//'   \item{beta0}{A numeric vector of intercepts for each lambda.}
//'   \item{beta}{A matrix of coefficients (p x length(lambda)).}
//'   \item{dev}{A numeric vector of deviance values.}
//'   \item{eeta}{A matrix of fitted linear predictors (n x length(lambda)).}
//'   \item{iter}{An integer vector of iteration counts for each lambda.}
//' }
// [[Rcpp::export]]
Rcpp::List cd_glm_fit_ssr(
  const arma::mat& X, const arma::vec& y, std::string family, const arma::vec& multiplier,
  std::string penalty, const arma::vec& lambda, double alpha, double gamma,
  double lam_max, int max_iter, int dfmax, double eps, bool user = false, bool warn = true
) {

  const int n = X.n_rows;
  const int p = X.n_cols;
  const int L = lambda.n_elem;

  arma::vec b0(L, arma::fill::zeros);
  arma::mat b(p, L, arma::fill::zeros);   // store column per lambda
  arma::vec Dev(L, arma::fill::zeros);
  arma::mat Eta(n, L, arma::fill::zeros);
  arma::vec df(L, arma::fill::zeros);
  arma::ivec iter(L, arma::fill::zeros);

  // working
  double a0 = 0.0;
  arma::vec a(p, arma::fill::zeros);         // previous beta
  arma::vec r(n, arma::fill::zeros);
  arma::vec w(n, arma::fill::zeros);
  arma::vec s(n, arma::fill::zeros);
  arma::vec z(p, arma::fill::zeros);
  arma::vec eta(n, arma::fill::zeros);
  arma::ivec e1(p, arma::fill::zeros); // active set
  arma::ivec e2(p, arma::fill::zeros); // strong set

  int tot_iter = 0;

  // ---- Initialization ----
  const double ybar = arma::mean(y);
  double nullDev = 0.0;
  double TOLERANCE = 1e-8;

  if (family == "binomial") {
    const double yb = std::min(std::max(ybar, TOLERANCE), 1.0 - TOLERANCE);
    a0 = std::log(yb / (1.0 - yb));
    b0(0) = a0;
    for (int i = 0; i < n; ++i) {
      nullDev -= 2.0 * (y(i) * std::log(yb) + (1.0 - y(i)) * std::log(1.0 - yb));
    }
  } else if (family == "poisson") {
    const double yb = std::max(ybar, TOLERANCE);
    a0 = std::log(yb);
    b0(0) = a0;
    for (int i = 0; i < n; ++i) {
      if (y(i) != 0.0) nullDev += 2.0 * (y(i) * std::log(y(i) / yb) + yb - y(i));
      else nullDev += 2.0 * yb;
    }
  } else {
    stop("family must be 'binomial' or 'poisson'");
  }

  for (int i = 0; i < n; i++) s(i) = y(i) - ybar;
  eta.fill(a0);
  for (int j = 0; j < p; j++) z(j) = crossprod(X, s, n, j) / (double)n;

  int lstart = 0;
  if (!user) {
    lstart = 1;
    Dev(0) = nullDev;
    Eta.col(0) = eta;
  }

  // ---- path ----
  for (int l = lstart; l < L; ++l) {
    Rcpp::checkUserInterrupt();

    if (l != 0) {
      a0 = b0(l - 1);
      a = b.col(l - 1);
      // eta = Eta.col(l - 1);        
    } else {
      a0 = b0(0);
      a = b.col(0);
    }

    int nv = 0;
    for (int j = 0; j < p; ++j) if (a(j) != 0.0) nv++;
    if ((nv > dfmax) || (tot_iter >= max_iter)) {
      for (int ll = l; ll < L; ++ll) iter(ll) = NA_INTEGER;
      break;
    }

    // Define strong set through sequential strong rule
    const double lmax = arma::abs(z).max();
    for (int j = 0; j < p; j++) {
      ssr_glm(e2, X, z, j, n, l, penalty, lambda, lmax, alpha, gamma, multiplier);
    }

    // Outer loop: fit until no violations (strong then rest)
    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        while (tot_iter < max_iter) {

          iter(l)++;
          tot_iter++;
          df(l) = 0.0;

          Dev(l) = 0.0;

          // build working quantities and deviance
          if (family == "binomial") {
            for (int i = 0; i < n; ++i) {
              const double mu = p_binomial(eta(i));
              const double wi = std::max(mu * (1.0 - mu), 1e-4);
              w(i) = wi;
              s(i) = y(i) - mu;
              r(i) = s(i) / wi;

              Dev(l) -= 2.0 * ( y(i) * std::log(std::max(mu, TOLERANCE))
                              + (1.0 - y(i)) * std::log(std::max(1.0 - mu, TOLERANCE)) );
            }
          } else { // poisson
            for (int i = 0; i < n; ++i) {
              const double mu = std::max(std::exp(eta(i)), TOLERANCE);
              w(i) = mu;
              s(i) = y(i) - mu;
              r(i) = s(i) / mu;

              if (y(i) != 0.0) Dev(l) += 2.0 * (y(i) * std::log(std::max(y(i) / mu, TOLERANCE)) + mu - y(i));
              else  Dev(l) += 2 * (mu - y(i));
            }
          }

          // saturation check
          if (Dev(l) / nullDev < 0.01) {
            if (warn) Rcpp::warning("Model saturated; exiting...");
            for (int ll = l; ll < L; ++ll) iter(ll) = NA_INTEGER;
            tot_iter = max_iter;
            break;
          }

          // Intercept update
          const double u0 = arma::dot(w, r);     
          const double v0 = arma::sum(w);        
          b0(l) = u0 / v0 + a0;                  
          const double si = b0(l) - a0;

          if (si != 0.0) {
            r -= si;
            eta += si;
          }

          double maxChange = std::fabs(si) * v0 / (double)n;

          // covariate updates
          for (int j = 0; j < p; ++j) {
            const double lam1 = lambda(l) * multiplier(j) * alpha;
            const double lam2 = lambda(l) * multiplier(j) * (1.0 - alpha);

            if (e1(j) == 1) {
              cd_glm(b, X, r, w, eta, j, a, n, l, penalty, lam1, lam2, gamma, df, maxChange);
            }
          }

          // convergence check
          a0 = b0(l);
          a  = b.col(l);
          if (maxChange < eps) break;
        }

        if (tot_iter >= max_iter || iter(l) == NA_INTEGER) break;

        // ---- scan strong set violations (e2==1, e1==0): not active but in strong set ----
        int violations = 0;
        for (int j = 0; j < p; ++j) {
          if (e1(j) == 0 && e2(j) == 1) {
            z(j) = crossprod(X, s, n, j) / (double)n;
            const double lam1 = lambda(l) * multiplier(j) * alpha;
            if (std::fabs(z(j)) > lam1) {
              e1(j) = 1;
              e2(j) = 1;
              violations++;
            }
          }
        }
        if (violations == 0) break;
      }

      // ---- scan rest violations (e2==0) not in the strong set ----
      int violations = 0;
      for (int j = 0; j < p; ++j) {
        if (e2(j) == 0) {
          z(j) = crossprod(X, s, n, j) / (double)n;
          const double lam1 = lambda(l) * multiplier(j) * alpha;
          if (std::fabs(z(j)) > lam1) {
            e1(j) = 1;
            e2(j) = 1;
            violations++;
          }
        }
      }

      // if no violations, store eta and move to next lambda
      if (violations == 0) {
        Eta.col(l) = eta;
        break;
      }
    } // end outer while
  }   // end path

  return Rcpp::List::create(
    _["beta0"] = b0,
    _["beta"]  = b,        
    _["dev"]   = Dev,
    _["eeta"]   = Eta,      // linear predictor
    _["iter"]  = iter    // return number of iterations
  );
}


// -- Coordinate descent update for weighted GLM --

void cd_wglm(
  arma::mat &b, const arma::mat &X, arma::vec &r, const arma::vec &w, const arma::vec &wt,
  arma::vec &eta, int j, const arma::vec &a, int n_obs, int l,  
  const std::string &penalty, double lam1, double lam2, double gamma,
  arma::vec &df, double &maxChange
) {

  double v = wsqsum(X, w % wt, n_obs, j); 
  double u = wcrossprod(X, r, w % wt, n_obs, j) + v * a(j);
  
  double len = 0.0;
  if (penalty == "LASSO")  len = lasso(u, lam1, lam2, v);
  else if (penalty == "SCAD") len = SCAD(u, lam1, lam2, gamma, v); 
  else if (penalty == "MCP")  len = MCP(u, lam1, lam2, gamma, v);
  
  if (len != 0 || a(j) != 0) {
    b(j, l) = len;
    double shift = b(j, l) - a(j);
    // Update residual and eta:
    if (shift != 0) {
      r -= X.col(j) * shift;
      eta += X.col(j) * shift;
      if (fabs(shift) * sqrt(v) > maxChange) maxChange = fabs(shift) * sqrt(v);
    }
  }
  
  if (len != 0) df(l) += 1.0; 
}


//' Coordinate descent for weighted GLM with sequential strong rule
//'
//' Coordinate-descent algorithm for penalised weighted generalised linear models
//' (binomial or Poisson) with sequential strong rule (SSR) screening over a
//' sequence of penalty values.
//'
//' @param X A numeric matrix of predictors (n x p).
//' @param y A numeric response vector of length n.
//' @param wt A numeric vector of observation weights (length n).
//' @param family A string: "binomial" or "poisson".
//' @param multiplier A numeric vector of per-variable penalty multipliers (length p).
//' @param penalty A string specifying the penalty: "LASSO", "MCP", or "SCAD".
//' @param lambda A numeric vector of penalty values, sorted in decreasing order.
//' @param alpha A numeric scalar for the elastic net mixing parameter.
//' @param gamma A numeric scalar for MCP and SCAD penalties.
//' @param lam_max A numeric scalar for the maximum lambda value.
//' @param max_iter An integer specifying the maximum total iterations.
//' @param dfmax An integer specifying the maximum number of selected variables.
//' @param eps A numeric scalar for the convergence tolerance.
//' @param user A logical flag for user-specified lambda starting behaviour.
//' @param warn A logical flag for saturation warnings.
//'
//' @return A list with beta0, beta, dev, eeta, and iter.
// [[Rcpp::export]]
Rcpp::List cd_wglm_fit_ssr(
  const arma::mat& X, const arma::vec& y, const arma::vec &wt, std::string family, const arma::vec& multiplier,
  std::string penalty, const arma::vec& lambda, double alpha, double gamma,
  double lam_max, int max_iter, int dfmax, double eps, bool user = false, bool warn = true
) {

  const int n_obs = X.n_rows;
  const int p = X.n_cols;
  const int L = lambda.n_elem;

  arma::vec b0(L, arma::fill::zeros);
  arma::mat b(p, L, arma::fill::zeros);   // store column per lambda
  arma::vec Dev(L, arma::fill::zeros);
  arma::mat Eta(n_obs, L, arma::fill::zeros);
  arma::vec df(L, arma::fill::zeros);
  arma::ivec iter(L, arma::fill::zeros);

  // working
  double a0 = 0.0;
  arma::vec a(p, arma::fill::zeros);         // previous beta
  arma::vec r(n_obs, arma::fill::zeros);
  arma::vec w(n_obs, arma::fill::zeros);
  arma::vec s(n_obs, arma::fill::zeros);
  arma::vec z(p, arma::fill::zeros);
  arma::vec eta(n_obs, arma::fill::zeros);
  arma::ivec e1(p, arma::fill::zeros); // active set
  arma::ivec e2(p, arma::fill::zeros); // strong set

  int tot_iter = 0;
  double TOLERANCE = 1e-8;

  // ---- Initialization ----
  const double ybar = arma::dot(wt, y);
  double nullDev = 0.0;

  if (family == "binomial") {
    const double yb = std::min(std::max(ybar, TOLERANCE), 1.0 - TOLERANCE);
    a0 = std::log(yb / (1.0 - yb));
    b0(0) = a0;
    for (int i = 0; i < n_obs; ++i) {
      nullDev -= 2.0 * wt(i) * (y(i) * std::log(yb) + (1.0 - y(i)) * std::log(1.0 - yb));
    }
  } else if (family == "poisson") {
    const double yb = std::max(ybar, TOLERANCE);
    a0 = std::log(yb);
    b0(0) = a0;
    for (int i = 0; i < n_obs; ++i) {
      if (y(i) != 0.0) nullDev += 2.0 * wt(i) * (y(i) * std::log(y(i) / yb) + yb - y(i));
      else nullDev += 2.0 * wt(i) * yb;
    }
  } else {
    Rcpp::stop("family must be 'binomial' or 'poisson'");
  }

  for (int i = 0; i < n_obs; i++) s(i) = y(i) - ybar;
  eta.fill(a0);
  for (int j = 0; j < p; j++) z(j) = wcrossprod(X, s, wt, n_obs, j);

  int lstart = 0;
  if (!user) {
    lstart = 1;
    Dev(0) = nullDev;
    Eta.col(0) = eta;
  }

  // ---- path ----
  for (int l = lstart; l < L; ++l) {
    Rcpp::checkUserInterrupt();

    if (l != 0) {
      a0 = b0(l - 1);
      a = b.col(l - 1);
      // eta = Eta.col(l - 1);        
    } else {
      a0 = b0(0);
      a = b.col(0);
    }

    int nv = 0;
    for (int j = 0; j < p; ++j) if (a(j) != 0.0) nv++;
    if ((nv > dfmax) || (tot_iter >= max_iter)) {
      for (int ll = l; ll < L; ++ll) iter(ll) = NA_INTEGER;
      break;
    }

    // Define strong set through sequential strong rule
    const double lmax = arma::abs(z).max();
    for (int j = 0; j < p; j++) {
      ssr_glm(e2, X, z, j, n_obs, l, penalty, lambda, lmax, alpha, gamma, multiplier);
    }

    // Outer loop: fit until no violations (strong then rest)
    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        while (tot_iter < max_iter) {

          iter(l)++;
          tot_iter++;
          df(l) = 0.0;

          Dev(l) = 0.0;

          // build working quantities and deviance
          if (family == "binomial") {
            for (int i = 0; i < n_obs; ++i) {
              const double mu = p_binomial(eta(i));
              const double wi = std::max(mu * (1.0 - mu), TOLERANCE);
              w(i) = wi;
              s(i) = y(i) - mu;
              r(i) = s(i) / wi;

              Dev(l) -= 2.0 * wt(i) * ( y(i) * std::log(std::max(mu, TOLERANCE))
                              + (1.0 - y(i)) * std::log(std::max(1.0 - mu, TOLERANCE)) );
            }
          } else { // poisson
            for (int i = 0; i < n_obs; ++i) {
              const double mu = std::max(std::exp(eta(i)), TOLERANCE);
              w(i) = mu;
              s(i) = y(i) - mu;
              r(i) = s(i) / mu;

              if (y(i) != 0.0) Dev(l) += 2.0 * wt(i) * (y(i) * std::log(std::max(y(i) / mu, TOLERANCE)) + mu - y(i));
              else  Dev(l) += 2 * wt(i) * (mu - y(i));
            }
          }

          // saturation check
          if (Dev(l) / nullDev < 0.01) {
            if (warn) Rcpp::warning("Model saturated; exiting...");
            for (int ll = l; ll < L; ++ll) iter(ll) = NA_INTEGER;
            tot_iter = max_iter;
            break;
          }

          // Intercept update
          const double u0 = arma::dot(w % wt, r);     
          const double v0 = arma::sum(w % wt);        
          b0(l) = u0 / v0 + a0;                  
          const double si = b0(l) - a0;

          if (si != 0.0) {
            r -= si;
            eta += si;
          }

          double maxChange = std::fabs(si) * v0;

          // covariate updates
          for (int j = 0; j < p; ++j) {
            const double lam1 = lambda(l) * multiplier(j) * alpha;
            const double lam2 = lambda(l) * multiplier(j) * (1.0 - alpha);

            if (e1(j) == 1) {
              cd_wglm(b, X, r, w, wt, eta, j, a, n_obs, l, penalty, lam1, lam2, gamma, df, maxChange);
            }
          }

          // convergence check
          a0 = b0(l);
          a  = b.col(l);
          if (maxChange < eps) break;
        }

        if (tot_iter >= max_iter || iter(l) == NA_INTEGER) break;

        // ---- scan strong set violations (e2==1, e1==0): not active but in strong set ----
        int violations = 0;
        for (int j = 0; j < p; ++j) {
          if (e1(j) == 0 && e2(j) == 1) {
            z(j) = wcrossprod(X, s, wt, n_obs, j);
            const double lam1 = lambda(l) * multiplier(j) * alpha;
            if (std::fabs(z(j)) > lam1) {
              e1(j) = 1;
              e2(j) = 1;
              violations++;
            }
          }
        }
        if (violations == 0) break;
      }

      // ---- scan rest violations (e2==0) not in the strong set ----
      int violations = 0;
      for (int j = 0; j < p; ++j) {
        if (e2(j) == 0) {
          z(j) = wcrossprod(X, s, wt, n_obs, j);
          const double lam1 = lambda(l) * multiplier(j) * alpha;
          if (std::fabs(z(j)) > lam1) {
            e1(j) = 1;
            e2(j) = 1;
            violations++;
          }
        }
      }

      // if no violations, store eta and move to next lambda
      if (violations == 0) {
        Eta.col(l) = eta;
        break;
      }
    } // end outer while
  }   // end path

  return Rcpp::List::create(
    _["beta0"] = b0,
    _["beta"]  = b,        
    _["dev"]   = Dev,
    _["eeta"]   = Eta,      // linear predictor
    _["iter"]  = iter    // return number of iterations
  );
}