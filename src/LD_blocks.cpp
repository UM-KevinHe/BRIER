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

//' LD Block Structured SNP grouping
//'
//' Calculate LD-block group boundaries for a list of SNPs given a pre-computed
//' block definition (Berisa et al.). The inputs must be sorted by chromosome and
//' position.
//'
//' @param chr Numeric vector of chromosome IDs (length n).
//' @param pos Numeric vector of base-pair positions (length n).
//' @param LDB Numeric matrix (m x 3) of (chr, start, stop) for each LD block.
//' @return Integer vector of length K+1 giving the 0-based index in (chr, pos)
//'         where each block starts, including the end marker at n.
//' @export
// [[Rcpp::export]]
arma::ivec LD_blocks(const arma::vec& chr, const arma::vec& pos, const arma::mat& LDB) {

  int n = pos.n_elem;
  int m = (int)LDB.n_rows;
  std::vector<int> boundaries;
  boundaries.reserve(n + 1);

  // first block always starts at index 0
  boundaries.push_back(0);

  // locate the block for the first SNP
  int j = 0;
  while (j < m && pos(0) >= LDB(j, 2)) {
    ++j;
  }
  double current_chr = chr(0);

  // scan the rest of the SNPs
  for (int i = 1; i < n; i++) {
    // chromosome change -> new block
    if (chr(i) != current_chr) {
      boundaries.push_back(i);
      current_chr = chr[i];
      // find new block for this SNP
      while (j < m && !(LDB(j, 0) == current_chr && pos(i) < LDB(j, 2))) {
        j++;
      }
    }
    // past the end of current block -> new block
    if (j < m && pos(i) >= LDB(j, 2)) {
      boundaries.push_back(i);
      while (j < m && pos(i) >= LDB(j, 2)) {
        j++;
      }
    }
  }

  // append end marker
  boundaries.push_back(n);

  // copy into Armadillo vector
  arma::ivec out(boundaries.size());
  for (size_t k = 0; k < boundaries.size(); ++k) {
    out(k) = boundaries[k];
  }
  return out;
}


//' Block-Structured LD Calculation
//'
//' Calculate LD matrix within blocks using soft-thresholding.
//'
//' @param X N x p standardized genotype matrix.
//' @param blk Vector of block boundaries of length K+1 (0-based).
//' @param tau Soft-threshold parameter.
//' @return p x p sparse LD matrix (arma::sp_mat).
//' @export
// [[Rcpp::export]]
arma::sp_mat LD_sigma(const arma::mat& X, const arma::uvec& blk, double tau) {

  int N = X.n_rows;
  int p = X.n_cols;
  int K = blk.n_elem - 1;  // number of blocks

  // Coordinate lists
  std::vector<arma::uword> I, J;
  std::vector<double> V;
  I.reserve(p * 10);
  J.reserve(p * 10);
  V.reserve(p * 10);

  // Loop over blocks
  for (int k = 0; k < K; k++) {
    int start = blk(k);
    int end   = blk(k+1);    // one‐past last
    int sz    = end - start;

    // Compute blockwise LD = (X_block^T * X_block) / N
    arma::mat sub = X.cols(start, end - 1);
    arma::mat mat = sub.t() * sub / double(N);

    // Diagonal entries = 1
    for (int i = 0; i < sz; i++) {
      int ii = start + i;
      I.push_back(ii);
      J.push_back(ii);
      V.push_back(1.0);
    }

    // Off‐diagonals: soft‐thresholded
    for (int i = 0; i < sz; ++i) {
      for (int j = i + 1; j < sz; ++j) {
        double val = Soft_thresh(mat(i, j), tau);
        if (val != 0.0) {
          int ii = start + i;
          int jj = start + j;
          // symmetric entries
          I.push_back(ii); J.push_back(jj); V.push_back(val);
          I.push_back(jj); J.push_back(ii); V.push_back(val);
        }
      }
    }
  }

  // Build arma::imat of coordinates (2×nnz)
  arma::umat coords(2, I.size());
  for (size_t idx = 0; idx < I.size(); ++idx) {      // previously: int idx, and I(idx)
    coords(0, idx) = I[idx];                         // previously: static_cast<arma::uword>(I(idx));
    coords(1, idx) = J[idx];
  }

  arma::vec vals(V.size());
  for (size_t i = 0; i < V.size(); ++i) {
    vals(i) = V[i];
  }

  // Construct sparse matrix
  arma::sp_mat S_out(coords, vals, p, p);
  return S_out;
}

