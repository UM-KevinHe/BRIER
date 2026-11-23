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


# -- LD_blocks: the block of the first SNP, chromosome order, edges (1.7.8) ---
#
# Before 1.7.8 the first SNP's block was looked up without its chromosome, so a
# SNP set that did not start on the first chromosome of LDB was cut at that
# chromosome's boundaries (found on chromosome 7 of real data, 2026-11-23).

ldb_two <- cbind(
  c(1, 1, 7, 7, 7),
  c(0, 1000, 0, 500, 1500),
  c(1000, 2000, 500, 1500, 3000)
)

test_that("LD_blocks: a set that starts on another chromosome gets its own blocks", {
  b <- as.integer(LD_blocks(c(7, 7, 7, 7, 7), c(100, 200, 600, 700, 1600), ldb_two))
  expect_equal(b, c(0L, 2L, 4L, 5L))
})

test_that("LD_blocks: the same answer as before for a set that starts on chromosome 1", {
  d <- ld_panel()
  b <- as.integer(LD_blocks(d$SNP.info$CHR, d$SNP.info$BP, d$LDB))
  expect_equal(b, c(0L, 4L, 7L, 12L))
})

test_that("LD_blocks: chromosome order and LDB row order do not matter", {
  b1 <- as.integer(LD_blocks(c(7, 7, 1, 1), c(100, 600, 100, 1100), ldb_two))
  expect_equal(b1, c(0L, 1L, 2L, 3L, 4L))
  b2 <- as.integer(LD_blocks(c(7, 7, 1, 1), c(100, 600, 100, 1100), ldb_two[c(5, 2, 4, 1, 3), ]))
  expect_equal(b2, b1)
})

test_that("LD_blocks: beyond the last block, and a chromosome without blocks", {
  # chr 1 SNP at 2500 lies beyond chr 1's last stop: it stays in that block, and chr 7 is unaffected
  b <- as.integer(LD_blocks(c(1, 1, 7, 7), c(1500, 2500, 100, 600), ldb_two))
  expect_equal(b, c(0L, 2L, 3L, 4L))
  # chr 3 has no row: its SNPs are one block
  b3 <- as.integer(LD_blocks(c(3, 3, 3, 7), c(100, 900, 5000, 100), ldb_two))
  expect_equal(b3, c(0L, 3L, 4L))
})

test_that("calLD: chromosome 7 alone with the whole bundled block file equals its own rows", {
  set.seed(7)
  n <- 150
  ldb <- utils::read.table(getLDB(ancestry = "AFR", build = "hg38"), header = TRUE)
  ldb[, 1] <- as.numeric(sub("^chr", "", ldb[, 1]))
  ldb <- as.matrix(ldb)
  b7 <- ldb[ldb[, 1] == 7, , drop = FALSE]
  bp <- sort(c(b7[1:4, 2] + 1000, b7[1:4, 2] + 5000))
  X <- matrix(stats::rbinom(n * length(bp), 2, 0.3), n, length(bp))
  info <- data.frame(CHR = 7, BP = bp)
  full <- calLD(X, SNP.info = info, LDB = ldb, tau = 0)
  own <- calLD(X, SNP.info = info, LDB = b7, tau = 0)
  expect_equal(as.integer(full$blk), c(0L, 2L, 4L, 6L, 8L))
  expect_equal(as.matrix(full$XtX), as.matrix(own$XtX))
})

test_that("calLD: SNP.info that is not grouped and sorted is an error", {
  d <- ld_panel()
  i <- c(2, 1, 3:12)
  expect_error(calLD(d$X[, i], SNP.info = d$SNP.info[i, ], LDB = d$LDB), "sorted by position")
  j <- c(1:4, 8:12, 5:7)
  expect_error(calLD(d$X[, j], SNP.info = d$SNP.info[j, ], LDB = d$LDB), "grouped by chromosome")
})
