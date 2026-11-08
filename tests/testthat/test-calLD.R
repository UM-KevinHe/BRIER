# Tests for calLD with LD blocks.
#
# The block boundaries returned by LD_blocks are column indices, so they must be
# computed on the columns that survive the removal of constant columns. A
# constant SNP early in the panel used to shift every later block onto the wrong
# columns. The right answer is known by construction: calLD on a panel with a
# constant column must equal calLD on the same panel with that column removed.


# -- Fixtures ----------------------------------------------------------------
#
# Two chromosomes, two blocks on chromosome 1 and one on chromosome 2, with
# correlated SNPs inside each block so a misaligned block shows up as wrong
# entries, not only as a different sparsity pattern.

ld_panel <- function(n = 200, seed = 1) {
  set.seed(seed)
  blk_sizes <- c(4, 3, 5)
  X <- do.call(cbind, lapply(blk_sizes, function(k) {
    z <- stats::rnorm(n)
    sapply(seq_len(k), function(j) {
      as.numeric(stats::rbinom(n, 2, stats::plogis(0.8 * z + stats::rnorm(n, sd = 0.5))))
    })
  }))
  SNP.info <- data.frame(
    CHR = c(rep(1, 7), rep(2, 5)),
    BP  = c(100, 200, 300, 400, 1100, 1200, 1300, 100, 200, 300, 400, 500)
  )
  LDB <- cbind(c(1, 1, 2), c(0, 1000, 0), c(1000, 2000, 1000))
  list(X = X, SNP.info = SNP.info, LDB = LDB)
}


# -- Tests -------------------------------------------------------------------

test_that("calLD with blocks: a constant column leaves the other blocks intact", {
  d <- ld_panel()
  X_const <- d$X
  X_const[, 2] <- 1                       # constant SNP inside the first block

  ld_ref <- calLD(d$X[, -2], SNP.info = d$SNP.info[-2, ], LDB = d$LDB)
  ld_new <- suppressWarnings(calLD(X_const, SNP.info = d$SNP.info, LDB = d$LDB))

  expect_identical(ld_new$nz, c(1L, 3:12))
  expect_equal(as.matrix(ld_new$XtX), as.matrix(ld_ref$XtX), tolerance = 1e-12)
  expect_identical(ld_new$blk, ld_ref$blk)
})

test_that("calLD with blocks: a constant column at the start of a chromosome", {
  d <- ld_panel(seed = 2)
  X_const <- d$X
  X_const[, 8] <- 0                       # first SNP of chromosome 2

  ld_ref <- calLD(d$X[, -8], SNP.info = d$SNP.info[-8, ], LDB = d$LDB)
  ld_new <- suppressWarnings(calLD(X_const, SNP.info = d$SNP.info, LDB = d$LDB))

  expect_equal(as.matrix(ld_new$XtX), as.matrix(ld_ref$XtX), tolerance = 1e-12)
  expect_identical(ld_new$blk, ld_ref$blk)
})

test_that("calLD with blocks: no cross-block entries", {
  d <- ld_panel(seed = 3)
  ld <- calLD(d$X, SNP.info = d$SNP.info, LDB = d$LDB)
  R <- as.matrix(ld$XtX)
  blk_id <- rep(seq_along(c(4, 3, 5)), c(4, 3, 5))
  expect_true(all(R[outer(blk_id, blk_id, "!=")] == 0))
  expect_equal(unname(diag(R)), rep(1, 12))
})
