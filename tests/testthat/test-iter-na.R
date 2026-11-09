# Tests for the iteration counts returned by the C++ fitters.
#
# The package compiles with -DARMA_64BIT_WORD (src/Makevars), so arma::ivec is a
# 64-bit integer vector and RcppArmadillo used to return it to R as double: the
# NA_INTEGER that marks a lambda cut off by dfmax or max_iter arrived as
# -2147483648, is.na() missed it, and BRIERi.eta / BRIERs.eta kept the unfitted
# lambdas as all-zero columns. The fitters now return a true integer NA.


# -- Fixtures ----------------------------------------------------------------

iter_case <- function(n = 200, p = 30, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  b <- c(rep(0.6, 8), rep(0, p - 8))
  lp <- drop(X %*% b)
  y <- lp + stats::rnorm(n)
  Xs <- scale(X)
  R <- Matrix::Matrix(stats::cor(X), sparse = TRUE)
  R <- Matrix::forceSymmetric(R)
  list(
    X = X, Xs = Xs, y = y, y_bin = stats::rbinom(n, 1, stats::plogis(lp)),
    R = methods::as(R, "generalMatrix"), XtY = drop(stats::cor(X, y)),
    n = n, p = p
  )
}

expect_true_na_tail <- function(iter) {
  expect_type(iter, "integer")
  expect_true(anyNA(iter))
  expect_false(any(iter < 0, na.rm = TRUE))
  # once a lambda is cut off, every later one is too
  first_na <- which(is.na(iter))[1]
  expect_true(all(is.na(iter[first_na:length(iter)])))
}

quiet <- function(expr) {
  utils::capture.output(res <- suppressWarnings(expr))
  res
}


# -- Tests -------------------------------------------------------------------

test_that("the C++ fitters return a true NA for lambdas cut off by dfmax", {
  d <- iter_case()
  m <- rep(1, d$p)
  wt <- rep(1, d$n) / d$n
  yc <- d$y - mean(d$y)
  lam_i <- exp(seq(log(1), log(0.01), length.out = 20))

  # summary statistics
  lam_max <- BRIER:::maxlambda_summary(d$XtY, d$R, m, 1, 1e-4, 1e6)
  lam_s <- exp(seq(log(lam_max), log(0.01 * lam_max), length.out = 20))
  f <- BRIER:::cd_summary_fit(d$XtY, d$R, m, "LASSO", lam_s, 1, 3, 1e-4, 1e6, 2, FALSE)
  expect_true_na_tail(f$iter)
  f <- BRIER:::cd_summary_fit_ssr(d$XtY, d$R, m, "LASSO", lam_s, 1, 3, lam_max, 1e-4, 1e6, 2, FALSE)
  expect_true_na_tail(f$iter)

  # weighted gaussian
  f <- quiet(BRIER:::cd_wgaussian_fit_ssr(d$Xs, yc, wt, m, "LASSO", lam_i, 1, 3, lam_i[1], 1e6, 1e-4, 2))
  expect_true_na_tail(f$iter)
  f <- quiet(BRIER:::cd_wgaussian_fit(d$Xs, yc, wt, m, "LASSO", lam_i, 1, 3, 1e6, 1e-4, 2))
  expect_true_na_tail(f$iter)

  # weighted and unweighted glm
  lam_b <- exp(seq(log(0.3), log(0.003), length.out = 20))
  f <- quiet(BRIER:::cd_wglm_fit_ssr(d$Xs, d$y_bin, wt, "binomial", m, "LASSO", lam_b, 1, 3, lam_b[1], 1e6, 2, 1e-4))
  expect_true_na_tail(f$iter)
  f <- quiet(BRIER:::cd_wglm_fit(d$Xs, d$y_bin, wt, "binomial", m, "LASSO", lam_b, 1, 3, 1e6, 2, 1e-4))
  expect_true_na_tail(f$iter)
  f <- quiet(BRIER:::cd_glm_fit_ssr(d$Xs, d$y_bin, "binomial", m, "LASSO", lam_b, 1, 3, lam_b[1], 1e6, 2, 1e-4))
  expect_true_na_tail(f$iter)
})

test_that("the C++ fitters return a true NA for lambdas cut off by max_iter", {
  d <- iter_case(seed = 2)
  m <- rep(1, d$p)
  lam_max <- BRIER:::maxlambda_summary(d$XtY, d$R, m, 1, 1e-4, 1e6)
  lam_s <- exp(seq(log(lam_max), log(0.01 * lam_max), length.out = 20))
  f <- BRIER:::cd_summary_fit(d$XtY, d$R, m, "LASSO", lam_s, 1, 3, 1e-10, 5, d$p, FALSE)
  expect_true_na_tail(f$iter)
  f <- BRIER:::cd_summary_fit_ssr(d$XtY, d$R, m, "LASSO", lam_s, 1, 3, lam_max, 1e-10, 5, d$p, FALSE)
  expect_true_na_tail(f$iter)
})

test_that("BRIERi.eta and BRIERs.eta drop the lambdas cut off by nvar.max", {
  d <- iter_case(seed = 3)
  for (ssr in c(TRUE, FALSE)) {
    fi <- quiet(BRIERi.eta(d$X, d$y, nlambda = 30, nvar.max = 2, ssr = ssr))
    expect_lt(length(fi$lambda), 30)
    expect_false(anyNA(fi$iter))
    expect_false(any(fi$iter < 0))
    expect_identical(ncol(fi$beta), length(fi$lambda))

    fb <- quiet(BRIERi.eta(d$X, d$y_bin, family = "binomial", nlambda = 30, nvar.max = 2, ssr = ssr))
    expect_lt(length(fb$lambda), 30)
    expect_false(anyNA(fb$iter))
    expect_identical(ncol(fb$beta), length(fb$lambda))

    fs <- BRIERs.eta(d$R, d$XtY, nlambda = 30, nvar.max = 2, ssr = ssr)
    expect_lt(length(fs$lambda), 30)
    expect_false(anyNA(fs$iter))
    expect_false(any(fs$iter < 0))
    expect_identical(ncol(fs$beta), length(fs$lambda))
  }
})
