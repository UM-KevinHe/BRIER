# Tests for convertLD() and convertLDref().
#
# convertLD() only removes entries, so its answer is known by construction: the
# block-wise (chromosome-wise) restriction of the full in-sample LD must equal what
# calLD() computes directly with those blocks (one block per chromosome). The
# PRS-CS converter is checked on a tiny HDF5 panel written here, when rhdf5 is
# installed: allele alignment flips signs, missing SNPs are reported, not filled.


# -- Fixtures ----------------------------------------------------------------

conv_panel <- function(n = 300, seed = 2) {
  set.seed(seed)
  blk_sizes <- c(4, 3, 5)
  z0 <- stats::rnorm(n)
  X <- do.call(cbind, lapply(blk_sizes, function(k) {
    z <- 0.6 * z0 + stats::rnorm(n)        # shared factor: LD also BETWEEN blocks
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


# -- convertLD ----------------------------------------------------------------

test_that("convertLD to blocks equals calLD with the same blocks", {
  d <- conv_panel()
  full <- calLD(d$X)$XtX
  expect_gt(max(abs(full[1:4, 5:7])), 0.01)   # the fixture has LD across blocks
  got <- convertLD(full, d$SNP.info, to = "block", LDB = d$LDB)
  ref <- calLD(d$X, SNP.info = d$SNP.info, LDB = d$LDB)
  expect_s3_class(got, "LD")
  expect_equal(as.matrix(got$XtX), as.matrix(ref$XtX), tolerance = 1e-12, ignore_attr = TRUE)
  expect_identical(as.integer(got$blk), as.integer(ref$blk))
  expect_identical(got$nz, seq_len(12))
})

test_that("convertLD to chromosomes equals calLD with one block per chromosome", {
  d <- conv_panel()
  full <- calLD(d$X)$XtX
  got <- convertLD(full, d$SNP.info, to = "chr")
  ref <- calLD(d$X, SNP.info = d$SNP.info, LDB = cbind(c(1, 2), 0, 1e10))
  expect_equal(as.matrix(got$XtX), as.matrix(ref$XtX), tolerance = 1e-12, ignore_attr = TRUE)
  expect_identical(as.integer(got$blk), c(0L, 7L, 12L))
  expect_true(all(as.matrix(got$XtX)[1:7, 8:12] == 0))
  expect_equal(as.matrix(got$XtX)[1:7, 1:7], as.matrix(full)[1:7, 1:7], ignore_attr = TRUE)
})

test_that("convertLD accepts a base matrix and a dsCMatrix, and keeps dimnames", {
  d <- conv_panel()
  full <- as.matrix(calLD(d$X)$XtX)
  dimnames(full) <- list(paste0("rs", 1:12), paste0("rs", 1:12))
  a <- convertLD(full, d$SNP.info, to = "block", LDB = d$LDB)
  b <- convertLD(Matrix::forceSymmetric(Matrix::Matrix(full, sparse = TRUE)), d$SNP.info, to = "block", LDB = d$LDB)
  expect_s4_class(a$XtX, "dgCMatrix")
  expect_equal(as.matrix(a$XtX), as.matrix(b$XtX), ignore_attr = TRUE)
  expect_identical(rownames(a$XtX), paste0("rs", 1:12))
  # block to chromosome is a no-op on a matrix that is already block-wise
  again <- convertLD(a$XtX, d$SNP.info, to = "chr")
  expect_equal(as.matrix(again$XtX), as.matrix(a$XtX), ignore_attr = TRUE)
})

test_that("convertLD rejects unsorted SNPs, a missing LDB and a size mismatch", {
  d <- conv_panel()
  full <- calLD(d$X)$XtX
  bad <- d$SNP.info[c(2, 1, 3:12), ]
  expect_error(convertLD(full, bad, to = "chr"), "sorted")
  expect_error(convertLD(full, d$SNP.info, to = "block"), "LDB")
  expect_error(convertLD(full[1:11, 1:11], d$SNP.info, to = "chr"), "rows")
})


# -- convertLDref("prscs") ---------------------------------------------------

test_that("convertLDref reads a PRS-CS panel, aligns signs and reports missing SNPs", {
  skip_if_not_installed("rhdf5")
  dir <- file.path(tempdir(), "ldblk_test")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  h5 <- file.path(dir, "ldblk_test_chr1.hdf5")
  if (file.exists(h5)) { file.remove(h5) }
  rhdf5::h5createFile(h5)
  b1 <- matrix(c(1, 0.5, 0.2, 0.5, 1, 0.3, 0.2, 0.3, 1), 3)
  b2 <- matrix(c(1, -0.4, -0.4, 1), 2)
  rhdf5::h5createGroup(h5, "blk_1")
  rhdf5::h5write(b1, h5, "blk_1/ldblk")
  rhdf5::h5write(c("rs1", "rs2", "rs3"), h5, "blk_1/snplist")
  rhdf5::h5createGroup(h5, "blk_2")
  rhdf5::h5write(b2, h5, "blk_2/ldblk")
  rhdf5::h5write(c("rs4", "rs5"), h5, "blk_2/snplist")
  rhdf5::h5closeAll()
  utils::write.table(
    data.frame(CHR = 1, SNP = paste0("rs", 1:5), BP = 1:5 * 100,
               A1 = c("A", "C", "G", "T", "A"), A2 = c("G", "T", "A", "C", "C"), MAF = 0.2),
    file.path(dir, "snpinfo_test_hm3"), quote = FALSE, row.names = FALSE, sep = "\t"
  )
  s <- data.frame(
    SNP = c("rs1", "rs2", "rs3", "rs9", "rs4", "rs5"),
    CHR = 1,
    REF = c("G", "C", "A", "A", "C", "A"),     # rs2 swapped (panel A1 = our REF); rs5 mismatched
    ALT = c("A", "T", "G", "G", "T", "G")
  )
  r <- convertLDref(dir, s)
  expect_identical(r$keep, c(1L, 2L, 3L, 5L))
  expect_identical(sort(r$missing$SNP), c("rs5", "rs9"))
  expect_identical(r$missing$reason[r$missing$SNP == "rs9"], "not in the panel")
  expect_identical(r$missing$reason[r$missing$SNP == "rs5"], "allele mismatch")
  expect_identical(r$sign.flipped, 1L)
  R <- as.matrix(r$XtX)
  expect_equal(R["rs1", "rs2"], -0.5)       # sign flipped once
  expect_equal(R["rs1", "rs3"], 0.2)
  expect_equal(R["rs2", "rs3"], -0.3)
  expect_equal(R["rs4", "rs4"], 1)
  expect_equal(R["rs1", "rs4"], 0)          # different PRS-CS blocks: nothing between them
  expect_identical(r$blocks, 2L)
})
