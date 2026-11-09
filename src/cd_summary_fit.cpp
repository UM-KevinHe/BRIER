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


// -- Input checks shared by both fitters --

// The running gradient g = XtY - R b is updated with COLUMN j of R, so it
// tracks R b. The coordinate update reads (R b)_j from g, which is the
// gradient of 1/2 b'Rb only when R is symmetric, so symmetry is asserted.
//
// Choice of check: one O(nnz) merge over the CSC arrays with O(p) extra
// memory, and NO transposed copy. arma::approx_equal(ld_mat, ld_mat.t()) was
// rejected because ld_mat.t() materialises a second sparse matrix, which for
// chromosome LD (about 1.4e9 nonzeros) is about 22 GB per forked worker.
// Sampling was rejected because one asymmetric pair is enough to break it.
//
// Method: entries (i, j) with i > j are visited in column-major order; for a
// fixed column i their mirrors (j, i) appear in column i in increasing j,
// which is exactly the order in which they are visited. One cursor per column
// therefore matches every strictly lower entry with its strictly upper mirror
// in a single pass. An entry without a stored mirror is compared against 0.
// Tolerance: |R_ij - R_ji| <= 1e-12 * max(1, |R_ij|, |R_ji|); an LD matrix
// has entries bounded by 1, so this is in effect an absolute 1e-12, above the
// rounding noise of a symmetric LD build. It guards against gross asymmetry
// (one triangle stored, a missing mirror).

inline bool sym_close(double x, double y, double tol) {
  double s = std::max(1.0, std::max(fabs(x), fabs(y)));
  return fabs(x - y) <= tol * s;
}

void stop_asym(arma::uword i, arma::uword j) {
  stop(
    "ld_mat is not symmetric: R[%d, %d] != R[%d, %d]. Symmetrize it, e.g. (R + t(R)) / 2.",
    (long long) i + 1, (long long) j + 1, (long long) j + 1, (long long) i + 1
  );
}

void check_symmetric_ld(const arma::sp_mat& ld_mat) {
  const double sym_tol = 1e-12;
  const arma::uword p = ld_mat.n_cols;
  ld_mat.sync();
  // raw CSC arrays, not the column iterator: the merge keeps one independent cursor per column
  const arma::uword* cp = ld_mat.col_ptrs;
  const arma::uword* ri = ld_mat.row_indices;
  const double* vv = ld_mat.values;

  // cursor(i): next not yet matched entry of column i (starts at the top)
  std::vector<arma::uword> cursor(cp, cp + p);

  for (arma::uword j = 0; j < p; j++) {
    for (arma::uword k = cp[j]; k < cp[j + 1]; k++) {
      arma::uword i = ri[k];
      if (i <= j) continue;  // strictly lower entries only
      // upper entries of column i above row j have no lower mirror
      while (cursor[i] < cp[i + 1] && ri[cursor[i]] < j) {
        if (!sym_close(vv[cursor[i]], 0.0, sym_tol)) stop_asym(ri[cursor[i]], i);
        cursor[i]++;
      }
      if (cursor[i] < cp[i + 1] && ri[cursor[i]] == j) {
        if (!sym_close(vv[k], vv[cursor[i]], sym_tol)) stop_asym(i, j);
        cursor[i]++;
      } else {
        if (!sym_close(vv[k], 0.0, sym_tol)) stop_asym(i, j);
      }
    }
  }
  // leftover strictly upper entries have no lower mirror
  for (arma::uword i = 0; i < p; i++) {
    while (cursor[i] < cp[i + 1] && ri[cursor[i]] < i) {
      if (!sym_close(vv[cursor[i]], 0.0, sym_tol)) stop_asym(ri[cursor[i]], i);
      cursor[i]++;
    }
  }
}


//' Check an LD matrix once for the BRIERs fitters
//'
//' The O(nnz) LD checks of cd_summary_fit and cd_summary_fit_ssr (square,
//' finite, symmetric), exported so that a driver can validate each LD (or each
//' block) once, before mclapply over eta, and then call the fitters with
//' check_ld = FALSE.
//'
//' @param ld_mat A p x p sparse matrix representing the linkage disequilibrium (LD) structure.
//' @return Nothing; stops with an error if a check fails.
// [[Rcpp::export]]
void check_ld_summary(const arma::sp_mat& ld_mat) {
  if (ld_mat.n_rows != ld_mat.n_cols) {
    stop("ld_mat must be square (p x p).");
  }
  if (!ld_mat.is_finite()) {
    stop(
      "ld_mat contains non-finite values. This usually indicates "
      "zero-variance SNPs in the LD reference panel or division-by-zero "
      "in correlation computation. Check calLD() inputs and tau threshold."
    );
  }
  check_symmetric_ld(ld_mat);
}

// Same order as before (XtY, ld_mat, multiplier), plus the length checks:
// g = XtY is written with the unchecked g[it.row()] for rows up to p - 1, so a
// short XtY must stop here, not write out of bounds. Square is tested once,
// inside check_ld_summary when check_ld is TRUE.
void check_inputs_summary(
    const arma::vec& XtY, const arma::sp_mat& ld_mat, const arma::vec& multiplier,
    bool check_ld
) {
  if (!XtY.is_finite()) {
    stop("XtY contains non-finite values (NaN or Inf). Check preprocessing.");
  }
  if (XtY.n_elem != ld_mat.n_cols) {
    stop("length(XtY) must equal ncol(ld_mat).");
  }
  if (check_ld) {
    check_ld_summary(ld_mat);
  } else if (ld_mat.n_rows != ld_mat.n_cols) {
    stop("ld_mat must be square (p x p).");
  }
  if (!multiplier.is_finite()) {
    stop("multiplier contains non-finite values.");
  }
  if (multiplier.n_elem != ld_mat.n_cols) {
    stop("length(multiplier) must equal ncol(ld_mat).");
  }
  ld_mat.sync();  // CSC form is final from here on; cd_summary uses the no_sync iterators
}


// -- Small helpers --

// penalty string to a named code, once per call (unknown: len stays 0)
enum penalty_type { PEN_LASSO = 0, PEN_MCP = 1, PEN_SCAD = 2, PEN_UNKNOWN = -1 };

int penalty_code(const std::string& penalty) {
  if (penalty == "LASSO") return PEN_LASSO;
  if (penalty == "MCP")   return PEN_MCP;
  if (penalty == "SCAD")  return PEN_SCAD;
  return PEN_UNKNOWN;
}

// ever-active set as an index list in increasing j, the sweep order of e
void active_index(const arma::ivec& e, std::vector<int>& act) {
  act.clear();
  for (int j = 0; j < (int) e.n_elem; j++) {
    if (e(j) == 1) act.push_back(j);
  }
}

// g = XtY - R b.col(l), recomputed exactly from the ever-active columns only
// (b.col(l) is 0 elsewhere): cost is the sum of their nnz, not nnz(R)
void refresh_gradient(
    const arma::vec& XtY, const arma::sp_mat& ld_mat, const arma::mat& b, int l,
    const std::vector<int>& act, arma::vec& g
) {
  g = XtY;
  for (int k = 0; k < (int) act.size(); k++) {
    int j = act[k];
    double bj = b(j, l);
    if (bj != 0.0) {
      arma::sp_mat::const_col_iterator it_end = ld_mat.end_col_no_sync(j);  // built once, not per element
      for (arma::sp_mat::const_col_iterator it = ld_mat.begin_col_no_sync(j); it != it_end; ++it) {
        g[it.row()] -= bj * (*it);
      }
    }
  }
}


// -- Coordinate descent update for summary statistics, running gradient --
//
// z_j is read from the running gradient g = XtY - R b, z_j = g_j + b_j, instead
// of being recomputed as XtY(j) - arma::dot(ld_mat.col(j), b.col(l)) + b(j, l).
// Unit diagonal assumed (the "+ b_j" and the denominator 1). a(j) is set to
// the new b(j, l) here, so a == b.col(l) between calls and the sweep needs no
// full copy (a(j) is only read at j, once per sweep).

void cd_summary(
    const arma::sp_mat& ld_mat, arma::vec& g,
    arma::mat& b, int j, arma::vec& a, int l,
    int pen, double lam1, double lam2, double gamma,
    arma::vec& df, double& maxChange
) {

  double z = g(j) + b(j, l);

  double len = 0.0;
  if (pen == PEN_LASSO) len = lasso(z, lam1, lam2, 1);
  if (pen == PEN_MCP)   len = MCP(z, lam1, lam2, gamma, 1);
  if (pen == PEN_SCAD)  len = SCAD(z, lam1, lam2, gamma, 1);

  if (len != 0.0 || a(j) != 0.0) {
    b(j, l) = len;
    double shift = b(j, l) - a(j);
    if (fabs(shift) > maxChange) maxChange = fabs(shift);
    if (len != 0.0) df[l] += len / fabs(z);

    // running gradient g -= shift * R[, j]; no_sync: ld_mat synced in check_inputs_summary
    if (shift != 0.0) {
      arma::sp_mat::const_col_iterator it_end = ld_mat.end_col_no_sync(j);  // built once, not per element
      for (arma::sp_mat::const_col_iterator it = ld_mat.begin_col_no_sync(j); it != it_end; ++it) {
        g[it.row()] -= shift * (*it);
      }
    }
    a(j) = b(j, l);
  }
}


// Coordinate descent with a running gradient g = XtY - R b that is updated
// only when a coefficient moves. Sweep order, in-pass KKT updates (an
// entering b_j is written and g updated at once, so later z's in the same
// pass see it), ever-active set, warm starts, stopping rule and the dfmax /
// max_iter semantics are those of the original fitter, so for a symmetric
// ld_mat the iterates agree with it up to rounding. g is recomputed exactly
// from the ever-active columns once at the start of every lambda (after the
// warm start) to bound rounding drift, and dev is computed in O(p) from the
// converged gradient: R b = XtY - g, so b'Rb - 2 b'XtY = -b'(XtY + g).

//' Coordinate descent for BRIERs estimation procedure
//'
//' Coordinate-descent algorithm for penalized regression using summary
//' statistics and an LD reference panel, with a running gradient and no
//' screening.
//'
//' @param XtY A numeric vector of marginal correlations between genotype and outcome (length p).
//' @param ld_mat A p x p sparse matrix representing the linkage disequilibrium (LD) structure.
//'   Must be symmetric (checked unless check_ld = FALSE; then it must have passed
//'   check_ld_summary()) and is assumed to have a unit diagonal.
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
//' @param check_ld A logical flag. If TRUE (default), run the O(nnz) LD checks
//'   (finite, symmetric); FALSE skips them when ld_mat was validated once with
//'   check_ld_summary().
//'
//' @return A list with the following elements:
//' \describe{
//'   \item{beta}{A matrix of coefficients (p x length(lam)).}
//'   \item{iter}{An integer vector of iteration counts for each lambda (NA if not fitted).}
//'   \item{df}{A numeric vector of degrees of freedom for each lambda.}
//'   \item{dev}{A numeric vector of deviance values for each lambda.}
//' }
// [[Rcpp::export]]
List cd_summary_fit(
    const arma::vec& XtY, const arma::sp_mat& ld_mat,
    const arma::vec& multiplier, const std::string& penalty, const arma::vec& lam,
    double alpha, double gamma, double eps, int max_iter, int dfmax, bool user,
    bool check_ld = true
) {

  check_inputs_summary(XtY, ld_mat, multiplier, check_ld);

  int p = ld_mat.n_cols;
  int L = lam.n_elem;
  int tot_iter = 0;
  int pen = penalty_code(penalty);

  // Initialize
  arma::mat b(p, L, arma::fill::zeros);
  arma::ivec iter(L, arma::fill::zeros);
  arma::vec df(L, arma::fill::zeros);
  arma::vec dev(L, arma::fill::zeros);

  // Intermediate
  arma::vec a(p, arma::fill::zeros);
  arma::ivec e(p, arma::fill::zeros);
  std::vector<int> act;                // indices with e(j) == 1, increasing
  arma::vec g(p, arma::fill::zeros);   // running gradient XtY - R b.col(l)
  int lstart = user ? 0 : 1;
  int nv, violations;
  double maxChange;

  for (int l = lstart; l < L; ++l) {
    Rcpp::checkUserInterrupt();

    if (l != 0) {
      nv = 0;
      for (int j = 0; j < p; j++) {
        if (a(j) != 0) nv++;
      }
      if (nv > dfmax || tot_iter == max_iter) {
        for (int ll = l; ll < L; ll++) { iter(ll) = NA_INTEGER; }
        break;
      }

      // warm start
      b.col(l) = b.col(l - 1);
      a = b.col(l - 1);
    }

    refresh_gradient(XtY, ld_mat, b, l, act, g);

    while (tot_iter < max_iter) {
      // active set cycling
      while (tot_iter < max_iter) {
        iter(l) += 1;
        tot_iter++;
        df(l) = 0;
        maxChange = 0.0;

        for (int k = 0; k < (int) act.size(); k++) {
          int j = act[k];
          double l1 = alpha * lam(l) * multiplier(j);
          double l2 = (1 - alpha) * lam(l) * multiplier(j);
          cd_summary(ld_mat, g, b, j, a, l, pen, l1, l2, gamma, df, maxChange);
        }

        if (maxChange <= eps) break;
      }

      // KKT check for inactive variables (in-pass update)
      violations = 0;
      for (int j = 0; j < p; j++) {
        if (e(j) == 0) {
          double l1 = alpha * lam(l) * multiplier(j);
          double l2 = (1 - alpha) * lam(l) * multiplier(j);
          cd_summary(ld_mat, g, b, j, a, l, pen, l1, l2, gamma, df, maxChange);
          if (b(j, l) != 0) {
            e(j) = 1;
            violations++;
          }
        }
      }
      if (violations > 0) active_index(e, act);

      if (violations == 0) {
        // NOTE: also reached when the convergence loop stopped at max_iter with
        // maxChange > eps; dev is then set on an unconverged column.
        dev[l] = -arma::dot(a, XtY + g);
        break;
      }
    }
  }

  List result = List::create(
    _["beta"] = b,
    _["iter"] = iter_to_R(iter),
    _["df"]   = df,
    _["dev"]  = dev
  );
  return result;
}


// -- Sequential strong rule screening for summary statistics --
//
// Same cutoffs as BRIERi's ssr_wgaussian, with c_j = g(j), the gradient at the
// warm start b_hat(lam_{l-1}) (the ridge term vanishes at b_j = 0). For l == 0
// the previous lambda is lam_max; with user = FALSE the first fitted lambda is
// l = 1 and its previous lambda is lam(0). A variable with multiplier 0 is
// always kept by construction: the right side is 0 and fabs(z) + TOLERANCE > 0.

void ssr_summary(
    arma::ivec& e2, const arma::vec& g, int j, int l,
    int pen, const arma::vec& lam, double lam_max,
    double alpha, double gamma, const arma::vec& multiplier
) {
  double cutoff = 0.0;
  double TOLERANCE = 1e-8;

  double z = g(j);

  if (l != 0) {
    if (pen == PEN_LASSO) cutoff = 2 * lam(l) - lam(l - 1);
    if (pen == PEN_SCAD)  cutoff = lam(l) + gamma / (gamma - 2) * (lam(l) - lam(l - 1));
    if (pen == PEN_MCP)   cutoff = lam(l) + gamma / (gamma - 1) * (lam(l) - lam(l - 1));
  } else {
    if (pen == PEN_LASSO) cutoff = 2 * lam(l) - lam_max;
    if (pen == PEN_SCAD)  cutoff = lam(l) + gamma / (gamma - 2) * (lam(l) - lam_max);
    if (pen == PEN_MCP)   cutoff = lam(l) + gamma / (gamma - 1) * (lam(l) - lam_max);
  }

  if (fabs(z) + TOLERANCE > (cutoff * alpha * multiplier(j))) {
    e2(j) = 1; // not reject, in strong set
  } else {
    e2(j) = 0; // reject
  }
}


// cd_summary_fit plus BRIERi-style sequential strong rule (SSR) screening.
// At the start of each lambda the strong set e2 is computed from the gradient
// at the warm start. The loop nesting is that of BRIERi's
// cd_wgaussian_fit_ssr: cycle the ever-active set e to convergence; KKT check
// over the strong set not yet active (violators enter, repeat); then KKT
// check over every variable outside the strong set (violators enter, repeat)
// until there is no violation. Exactness (KKT stationarity to eps) comes from
// these checks; the strong rule itself is a heuristic.
//
// First fitted lambda: the rule assumes b = 0 is the solution at the previous
// lambda. That fails when some multiplier is 0 (the unpenalised coordinates
// are not 0 at lambda.max, and with user = FALSE lam(0) is never fitted, so
// b.col(0) stays 0). In that case EVERY variable is kept in the strong set at
// the first fitted lambda; the iterates are not changed.
//
// KKT checks UPDATE b IN PASS (as cd_summary_fit), they do not only flag (as
// BRIERi). With the running gradient the test is O(1) per variable. It also
// makes this function a strict generalisation of cd_summary_fit: if the strong
// set were every variable, the iterates, iter, df and dev would be identical.

//' Coordinate descent for BRIERs estimation procedure with sequential strong rule
//'
//' Coordinate-descent algorithm for penalized regression using summary
//' statistics and an LD reference panel, with a running gradient and
//' sequential strong rule (SSR) screening over a sequence of penalty values.
//'
//' @param XtY A numeric vector of marginal correlations between genotype and outcome (length p).
//' @param ld_mat A p x p sparse matrix representing the linkage disequilibrium (LD) structure.
//'   Must be symmetric (checked unless check_ld = FALSE; then it must have passed
//'   check_ld_summary()) and is assumed to have a unit diagonal.
//' @param multiplier A numeric vector of per-variable penalty multipliers (length p).
//'   Variables with multiplier = 0 are unpenalised and always in the strong set.
//' @param penalty A string specifying the penalty: "LASSO", "MCP", or "SCAD".
//' @param lam A numeric vector of penalty values, sorted in decreasing order.
//' @param alpha A numeric scalar for the elastic net mixing parameter.
//' @param gamma A numeric scalar for MCP and SCAD penalties.
//' @param lam_max A numeric scalar: lambda.max from
//'   maxlambda_summary(XtY, ld_mat, multiplier, alpha, eps, max_iter). The
//'   caller must pass the real value, computed even when lam is user-specified;
//'   it is the previous lambda for l = 0. A negative or non-finite value (such
//'   as BRIERi's user-lambda sentinel -1) is rejected.
//' @param eps A numeric scalar for the convergence tolerance.
//' @param max_iter An integer specifying the maximum total iterations allowed.
//' @param dfmax An integer specifying the maximum number of selected variables.
//' @param user A logical flag. If TRUE, lambda is user-specified and iteration starts
//'   from the first value; if FALSE, the first lambda is skipped.
//' @param check_ld A logical flag. If TRUE (default), run the O(nnz) LD checks
//'   (finite, symmetric); FALSE skips them when ld_mat was validated once with
//'   check_ld_summary().
//'
//' @return A list with the following elements:
//' \describe{
//'   \item{beta}{A matrix of coefficients (p x length(lam)).}
//'   \item{iter}{An integer vector of iteration counts for each lambda (NA if not fitted).}
//'   \item{df}{A numeric vector of degrees of freedom for each lambda.}
//'   \item{dev}{A numeric vector of deviance values for each lambda.}
//' }
// [[Rcpp::export]]
List cd_summary_fit_ssr(
    const arma::vec& XtY, const arma::sp_mat& ld_mat,
    const arma::vec& multiplier, const std::string& penalty, const arma::vec& lam,
    double alpha, double gamma, double lam_max,
    double eps, int max_iter, int dfmax, bool user,
    bool check_ld = true
) {

  check_inputs_summary(XtY, ld_mat, multiplier, check_ld);
  if (!std::isfinite(lam_max) || lam_max < 0) {
    stop(
      "lam_max must be a finite, non-negative lambda.max from maxlambda_summary(), "
      "computed even for a user lambda (the -1 sentinel is not accepted)."
    );
  }

  int p = ld_mat.n_cols;
  int L = lam.n_elem;
  int tot_iter = 0;
  int pen = penalty_code(penalty);
  bool any_unpen = arma::any(multiplier == 0);

  // Initialize
  arma::mat b(p, L, arma::fill::zeros);
  arma::ivec iter(L, arma::fill::zeros);
  arma::vec df(L, arma::fill::zeros);
  arma::vec dev(L, arma::fill::zeros);

  // Intermediate
  arma::vec a(p, arma::fill::zeros);
  arma::ivec e(p, arma::fill::zeros);    // ever-active set
  std::vector<int> act;                  // indices with e(j) == 1, increasing
  arma::ivec e2(p, arma::fill::zeros);   // strong set
  arma::vec g(p, arma::fill::zeros);     // running gradient XtY - R b.col(l)
  int lstart = user ? 0 : 1;
  int nv, violations;
  double maxChange;

  for (int l = lstart; l < L; ++l) {
    Rcpp::checkUserInterrupt();

    if (l != 0) {
      nv = 0;
      for (int j = 0; j < p; j++) {
        if (a(j) != 0) nv++;
      }
      if (nv > dfmax || tot_iter == max_iter) {
        for (int ll = l; ll < L; ll++) { iter(ll) = NA_INTEGER; }
        break;
      }

      // warm start
      b.col(l) = b.col(l - 1);
      a = b.col(l - 1);
    }

    refresh_gradient(XtY, ld_mat, b, l, act, g);

    // Set SSR rule (first fitted lambda with unpenalised variables: keep all)
    if (l == lstart && any_unpen) {
      e2.ones();
    } else {
      for (int j = 0; j < p; j++) {
        ssr_summary(e2, g, j, l, pen, lam, lam_max, alpha, gamma, multiplier);
      }
    }

    while (tot_iter < max_iter) {      // outer loop
      while (tot_iter < max_iter) {    // strong set loop
        // active set cycling
        while (tot_iter < max_iter) {  // convergence loop
          iter(l) += 1;
          tot_iter++;
          df(l) = 0;
          maxChange = 0.0;

          for (int k = 0; k < (int) act.size(); k++) {
            int j = act[k];
            double l1 = alpha * lam(l) * multiplier(j);
            double l2 = (1 - alpha) * lam(l) * multiplier(j);
            cd_summary(ld_mat, g, b, j, a, l, pen, l1, l2, gamma, df, maxChange);
          }

          if (maxChange <= eps) break;
        }

        // KKT check: in the strong set, not yet active (in-pass update)
        violations = 0;
        for (int j = 0; j < p; j++) {
          if (e(j) == 0 && e2(j) == 1) {
            double l1 = alpha * lam(l) * multiplier(j);
            double l2 = (1 - alpha) * lam(l) * multiplier(j);
            cd_summary(ld_mat, g, b, j, a, l, pen, l1, l2, gamma, df, maxChange);
            if (b(j, l) != 0) {
              e(j) = 1;
              violations++;
            }
          }
        }
        if (violations > 0) active_index(e, act);
        if (violations == 0) break;
      }
      // max_iter reached inside the strong set loop with pending violators:
      // leave this lambda unconverged (dev stays 0), as cd_summary_fit does
      if (violations > 0) break;

      // KKT check: outside the strong set (in-pass update). Ever-active
      // variables are cycled above, so only e(j) == 0 is checked here; BRIERi
      // also rechecks ever-active ones that fell out of e2, which only costs
      // an extra round.
      violations = 0;
      for (int j = 0; j < p; j++) {
        if (e(j) == 0 && e2(j) == 0) {
          double l1 = alpha * lam(l) * multiplier(j);
          double l2 = (1 - alpha) * lam(l) * multiplier(j);
          cd_summary(ld_mat, g, b, j, a, l, pen, l1, l2, gamma, df, maxChange);
          if (b(j, l) != 0) {
            e(j) = 1;
            e2(j) = 1;
            violations++;
          }
        }
      }
      if (violations > 0) active_index(e, act);

      if (violations == 0) {
        // NOTE: also reached when the convergence loop stopped at max_iter with
        // maxChange > eps; dev is then set on an unconverged column, as in
        // cd_summary_fit
        dev[l] = -arma::dot(a, XtY + g);
        break;
      }
    }
  }

  List result = List::create(
    _["beta"] = b,
    _["iter"] = iter_to_R(iter),
    _["df"]   = df,
    _["dev"]  = dev
  );
  return result;
}
