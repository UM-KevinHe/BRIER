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
//' block definition (Berisa et al.). The SNPs must be grouped by chromosome
//' (each chromosome in one contiguous run, in any chromosome order) and sorted by
//' position within each chromosome; \code{calLD} checks this. The rows of
//' \code{LDB} may be in any order.
//'
//' A SNP belongs to the first block of its own chromosome whose stop is above
//' its position; a SNP beyond the last block of its chromosome belongs to that
//' last block, and all SNPs of a chromosome with no row in \code{LDB} form one
//' block. Before 1.7.8 the first SNP's block was looked up without its
//' chromosome, so a SNP set that did not start on the first chromosome of
//' \code{LDB} was cut at that chromosome's boundaries on its first chromosome.
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
  if (n == 0) {
    arma::ivec out0(1);
    out0(0) = 0;
    return out0;
  }

  // the block stops of each chromosome, sorted
  std::map<double, std::vector<double> > stops;
  for (int r = 0; r < m; r++) {
    stops[LDB(r, 0)].push_back(LDB(r, 2));
  }
  for (auto& kv : stops) {
    std::sort(kv.second.begin(), kv.second.end());
  }

  // block of SNP i within its chromosome: the first stop above its position,
  // the last block when it lies beyond every stop, -1 when the chromosome has
  // no block at all
  auto block_of = [&](int i) -> long {
    auto it = stops.find(chr(i));
    if (it == stops.end()) return -1;
    const std::vector<double>& s = it->second;
    long k = (long)(std::upper_bound(s.begin(), s.end(), pos(i)) - s.begin());
    if (k >= (long)s.size()) k = (long)s.size() - 1;
    return k;
  };

  long prev_blk = block_of(0);
  for (int i = 1; i < n; i++) {
    long b = block_of(i);
    // a new block at every chromosome change and every change of block
    if (chr(i) != chr(i - 1) || b != prev_blk) {
      boundaries.push_back(i);
    }
    prev_blk = b;
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

