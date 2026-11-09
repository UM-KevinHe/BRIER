# Tests for the BRIERs fitters: cd_summary_fit (running gradient, no screening)
# and cd_summary_fit_ssr (the same plus the sequential strong rule).
#
# The strong rule is a heuristic and the KKT check repairs it, so for the LASSO
# both fitters must reach the same solution: identical supports, and objective
# and R b equal up to the convergence tolerance. With every strong rule cutoff
# negative nothing is screened out, and the ssr fitter then runs exactly the
# operations of the plain one.


# -- Fixtures ----------------------------------------------------------------
#
# Three LD blocks of correlated SNPs on one chromosome, block LD from calLD, and
# marginal correlations from a sparse true model on the same genotypes, so XtY is
# consistent with the LD.

ssr_summary_case <- function(n = 400, seed = 1) {
  set.seed(seed)
  blk_sizes <- c(20, 15, 25)
  X <- do.call(cbind, lapply(blk_sizes, function(k) {
    z <- stats::rnorm(n)
    sapply(seq_len(k), function(j) {
      as.numeric(stats::rbinom(n, 2, stats::plogis(0.8 * z + stats::rnorm(n, sd = 0.7))))
    })
  }))
  p <- ncol(X)
  SNP.info <- data.frame(
    CHR = 1,
    BP  = c(
      seq(100, by = 10, length.out = 20),
      seq(1100, by = 10, length.out = 15),
      seq(2100, by = 10, length.out = 25)
    )
  )
  LDB <- cbind(c(1, 1, 1), c(0, 1000, 2000), c(1000, 2000, 3000))
  ld <- calLD(X, SNP.info = SNP.info, LDB = LDB)

  b <- rep(0, p)
  b[c(3, 25, 40, 41)] <- c(0.5, -0.4, 0.3, 0.3)
  y <- drop(scale(X) %*% b) + stats::rnorm(n)
  list(R = ld$XtX, XtY = drop(stats::cor(X, y)), p = p)
}

# LASSO objective of the BRIERs fitters, one value per lambda (multiplier 1)
lasso_objective <- function(R, XtY, beta, lam) {
  Rb <- as.matrix(R %*% beta)
  0.5 * colSums(beta * Rb) - drop(crossprod(beta, XtY)) + lam * colSums(abs(beta))
}


# -- Tests -------------------------------------------------------------------

test_that("cd_summary_fit and cd_summary_fit_ssr reach the same LASSO solution", {
  d <- ssr_summary_case()
  m <- rep(1, d$p)
  eps <- 1e-10
  lam_max <- BRIER:::maxlambda_summary(d$XtY, d$R, m, 1, eps, 1e6)
  lam <- exp(seq(log(lam_max), log(0.01 * lam_max), length.out = 30))

  f0 <- BRIER:::cd_summary_fit(d$XtY, d$R, m, "LASSO", lam, 1, 3, eps, 1e6, d$p, FALSE)
  f1 <- BRIER:::cd_summary_fit_ssr(d$XtY, d$R, m, "LASSO", lam, 1, 3, lam_max, eps, 1e6, d$p, FALSE)

  expect_false(anyNA(f0$iter))
  expect_false(anyNA(f1$iter))
  expect_identical(f0$beta != 0, f1$beta != 0)
  expect_true(any(f0$beta != 0))
  expect_lt(
    max(abs(lasso_objective(d$R, d$XtY, f0$beta, lam) - lasso_objective(d$R, d$XtY, f1$beta, lam))),
    1e-8
  )
  expect_lt(max(abs(as.matrix(d$R %*% (f0$beta - f1$beta)))), 1e-8)
})

test_that("cd_summary_fit_ssr with nothing screened out equals cd_summary_fit", {
  d <- ssr_summary_case(seed = 2)
  m <- rep(1, d$p)
  lam_max <- BRIER:::maxlambda_summary(d$XtY, d$R, m, 1, 1e-10, 1e6)
  # user lambda: a huge lam_max makes the first cutoff, 2 lam(1) - lam_max,
  # hugely negative, and a ratio of 0.4 between lambdas keeps every later cutoff,
  # 2 lam(l) - lam(l - 1), negative too
  lam <- 0.9 * lam_max * 0.4^(0:7)

  f0 <- BRIER:::cd_summary_fit(d$XtY, d$R, m, "LASSO", lam, 1, 3, 1e-10, 1e6, d$p, TRUE)
  f1 <- BRIER:::cd_summary_fit_ssr(d$XtY, d$R, m, "LASSO", lam, 1, 3, 1e10, 1e-10, 1e6, d$p, TRUE)

  expect_true(any(f0$beta != 0))
  expect_identical(f1$beta, f0$beta)
  expect_identical(f1$iter, f0$iter)
  expect_identical(f1$df, f0$df)
  expect_identical(f1$dev, f0$dev)
})

test_that("BRIERs.eta gives the same LASSO coefficients with ssr TRUE and FALSE", {
  d <- ssr_summary_case(seed = 3)
  B <- matrix(stats::rnorm(d$p, sd = 0.05), ncol = 1)
  XtY.external <- as.matrix(d$R %*% B)

  for (eta in c(0, 1)) {
    ext <- if (eta == 0) NULL else XtY.external
    f_ssr <- BRIERs.eta(d$R, d$XtY, eta = eta, XtY.external = ext, nlambda = 30, ssr = TRUE)
    f_all <- BRIERs.eta(d$R, d$XtY, eta = eta, XtY.external = ext, nlambda = 30, ssr = FALSE)
    expect_identical(f_ssr$lambda, f_all$lambda)
    expect_identical(f_ssr$beta != 0, f_all$beta != 0)
    expect_lt(max(abs(f_ssr$beta - f_all$beta)), 1e-6)
  }

  # and through the BRIERs wrapper, serial
  ss <- data.frame(corr = d$XtY)
  g_ssr <- BRIERs(ss, d$R, eta.list = list(c(0, 1)), beta.external = B,
                  nlambda = 30, parallel = FALSE, ssr = TRUE)
  g_all <- BRIERs(ss, d$R, eta.list = list(c(0, 1)), beta.external = B,
                  nlambda = 30, parallel = FALSE, ssr = FALSE)
  for (k in seq_along(g_ssr$res)) {
    expect_lt(max(abs(g_ssr$res[[k]]$beta - g_all$res[[k]]$beta)), 1e-6)
  }

  # the switch really reaches the fitter: with the strong rule fitter disabled,
  # ssr = FALSE still fits and the default does not
  local_mocked_bindings(cd_summary_fit_ssr = function(...) stop("strong rule fitter called"))
  expect_error(
    BRIERs(ss, d$R, eta.list = list(0), beta.external = B, nlambda = 10, parallel = FALSE),
    "strong rule fitter called"
  )
  expect_s3_class(
    BRIERs(ss, d$R, eta.list = list(0), beta.external = B, nlambda = 10,
           parallel = FALSE, ssr = FALSE),
    "BRIER"
  )
})

test_that("BRIERs.eta with a user lambda uses the strong rule too", {
  d <- ssr_summary_case(seed = 4)
  f_def <- BRIERs.eta(d$R, d$XtY, nlambda = 20)
  lam <- f_def$lambda[5:15]
  f_ssr <- BRIERs.eta(d$R, d$XtY, lambda = lam, ssr = TRUE)
  f_all <- BRIERs.eta(d$R, d$XtY, lambda = lam, ssr = FALSE)
  expect_identical(f_ssr$beta != 0, f_all$beta != 0)
  expect_lt(max(abs(f_ssr$beta - f_all$beta)), 1e-6)
})

test_that("the BRIERs fitters reject an asymmetric LD, a wrong-length XtY and a negative lam_max", {
  d <- ssr_summary_case()
  m <- rep(1, d$p)
  lam <- c(0.1, 0.05)

  R_asym <- d$R
  R_asym[1, 2] <- R_asym[1, 2] + 0.01
  expect_error(
    BRIER:::cd_summary_fit(d$XtY, R_asym, m, "LASSO", lam, 1, 3, 1e-4, 1e6, d$p, TRUE),
    "not symmetric"
  )
  expect_error(
    BRIER:::cd_summary_fit_ssr(d$XtY, R_asym, m, "LASSO", lam, 1, 3, 0.2, 1e-4, 1e6, d$p, TRUE),
    "not symmetric"
  )
  expect_error(BRIERs.eta(R_asym, d$XtY, nlambda = 5), "not symmetric")

  expect_error(
    BRIER:::cd_summary_fit(d$XtY[-1], d$R, m, "LASSO", lam, 1, 3, 1e-4, 1e6, d$p, TRUE),
    "length\\(XtY\\)"
  )
  expect_error(
    BRIER:::cd_summary_fit_ssr(d$XtY[-1], d$R, m, "LASSO", lam, 1, 3, 0.2, 1e-4, 1e6, d$p, TRUE),
    "length\\(XtY\\)"
  )

  expect_error(
    BRIER:::cd_summary_fit_ssr(d$XtY, d$R, m, "LASSO", lam, 1, 3, -1, 1e-4, 1e6, d$p, TRUE),
    "lam_max"
  )
  expect_error(BRIERs.eta(d$R, d$XtY, nlambda = 5, ssr = NA), "ssr must be")
})
