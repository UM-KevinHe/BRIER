#ifndef BRIER_UTILS_H
#define BRIER_UTILS_H

#include <new>
#include <iterator>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>

// Simple soft-thresholding (scalar, no ridge/scale)
inline double Soft_thresh(double z, double l) {
  if (z >  l) return z - l;
  if (z < -l) return z + l;
  return 0.0;
}

// -- Thresholding operators --

inline double lasso(double z, double l1, double l2, double v) {
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return 0;
  else return s * (fabs(z) - l1) / (v * (1 + l2));
}

inline double MCP(double z, double l1, double l2, double gamma, double v) {
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return 0;
  else if (fabs(z) <= gamma * l1 * (1 + l2)) return s * (fabs(z) - l1) / (v * (1 + l2 - 1 / gamma));
  else return z / (v * (1 + l2));
}

inline double SCAD(double z, double l1, double l2, double gamma, double v) {
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return 0;
  else if (fabs(z) <= (l1 * (1 + l2) + l1)) return s * (fabs(z) - l1) / (v * (1 + l2));
  else if (fabs(z) <= gamma * l1 * (1 + l2)) return s * (fabs(z) - gamma * l1 / (gamma - 1)) / (v * (1 - 1 / (gamma - 1) + l2));
  else return z / (v * (1 + l2));
}

// -- Weighted linear algebra --

inline double wsqsum(const arma::mat &Z, const arma::vec &w, int n_obs, int j) {
  double val = 0;
  for (int i = 0; i < n_obs; i++) {
    val += w(i) * Z(i, j) * Z(i, j);
  }
  return val;
}

inline double wcrossprod(const arma::mat &Z, const arma::vec &r, const arma::vec &w, int n_obs, int j) {
  double val = 0;
  for (int i = 0; i < n_obs; i++) {
    val += w(i) * Z(i, j) * r(i);
  }
  return val;
}

inline double crossprod(const arma::mat &Z, const arma::vec &r, int n_obs, int j) {
  double val = 0;
  for (int i = 0; i < n_obs; i++) {
    val += Z(i, j) * r(i);
  }
  return val;
}

// -- Utilities --

inline int sum_rejections(const arma::ivec &x, int n_obs) {
  int val = 0;
  for (int i = 0; i < n_obs; i++) val += x(i);
  return val;
}

inline double p_binomial(double eta) {
  if (eta > 10) return 1;
  else if (eta < -10) return 0;
  else return exp(eta) / (1 + exp(eta));
}

// -- Deviance --

inline double calc_deviance(const arma::vec &y, const arma::vec &mu, const std::string &family) {
  double dev = 0.0;
  double TOLERANCE = 1e-8;
  int n_obs = y.n_elem;

  if (family == "binomial") {
    for (int i = 0; i < n_obs; ++i) {
      double m = std::min(std::max(mu(i), TOLERANCE), 1.0 - TOLERANCE);
      dev -= 2.0 * (y(i) * std::log(m) + (1.0 - y(i)) * std::log(1.0 - m));
    }
  } else if (family == "poisson") {
    for (int i = 0; i < n_obs; ++i) {
      double m = std::max(mu(i), TOLERANCE);
      if (y(i) > 0) dev += 2.0 * (y(i) * std::log(y(i) / m) - (y(i) - m));
      else dev += 2.0 * m;
    }
  } else if (family == "gaussian") {
    for (int i = 0; i < n_obs; ++i) {
      double r = y(i) - mu(i);
      dev += r * r;
    }
  }
  return dev;
}

inline double calc_wdeviance(const arma::vec &y, const arma::vec &mu, const arma::vec &wt, const std::string &family) {
  double dev = 0.0;
  double TOLERANCE = 1e-8;
  int n_obs = y.n_elem;

  if (family == "binomial") {
    for (int i = 0; i < n_obs; ++i) {
      double m = std::min(std::max(mu(i), TOLERANCE), 1.0 - TOLERANCE);
      dev -= 2.0 * wt(i) * (y(i) * std::log(m) + (1.0 - y(i)) * std::log(1.0 - m));
    }
  } else if (family == "poisson") {
    for (int i = 0; i < n_obs; ++i) {
      double m = std::max(mu(i), TOLERANCE);
      if (y(i) > 0) dev += 2.0 * wt(i) * (y(i) * std::log(y(i) / m) - (y(i) - m));
      else dev += 2.0 * m;
    }
  } else if (family == "gaussian") {
    for (int i = 0; i < n_obs; ++i) {
      double r = y(i) - mu(i);
      dev += r * r * wt(i);
    }
  }
  return dev;
}

#endif // BRIER_UTILS_H