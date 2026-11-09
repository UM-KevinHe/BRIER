# Tests for the BRIERi switch between the strong rule fitters
# (cd_wgaussian_fit_ssr, cd_wglm_fit_ssr) and their non-screening counterparts
# (cd_wgaussian_fit, cd_wglm_fit). For the LASSO both must reach the same
# solution along the path: identical supports and coefficients up to the
# convergence tolerance.


# -- Fixtures ----------------------------------------------------------------

ssr_individual_case <- function(n = 300, p = 40, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  b <- c(0.8, -0.6, 0.5, 0.4, -0.3, rep(0, p - 5))
  lp <- drop(X %*% b)
  list(
    X = X,
    y_gauss = lp + stats::rnorm(n),
    y_bin = stats::rbinom(n, 1, stats::plogis(lp)),
    # an external prediction from a perturbed true model, on the response scale
    ext_gauss = drop(X %*% (b + stats::rnorm(p, sd = 0.1))),
    ext_bin = stats::plogis(drop(X %*% (b + stats::rnorm(p, sd = 0.1))))
  )
}

fit_i <- function(X, y, ...) {
  suppressWarnings(suppressMessages(BRIERi.eta(X, y, nlambda = 30, ...)))
}

expect_same_path <- function(f_ssr, f_all) {
  expect_identical(f_ssr$lambda, f_all$lambda)
  expect_identical(f_ssr$beta != 0, f_all$beta != 0)
  expect_lt(max(abs(f_ssr$beta - f_all$beta)), 1e-6)
}


# -- Tests -------------------------------------------------------------------

test_that("BRIERi.eta gaussian LASSO: ssr TRUE and FALSE give the same path", {
  d <- ssr_individual_case()
  expect_same_path(
    fit_i(d$X, d$y_gauss, ssr = TRUE),
    fit_i(d$X, d$y_gauss, ssr = FALSE)
  )
  expect_same_path(
    fit_i(d$X, d$y_gauss, eta = 1, y.external = d$ext_gauss, ssr = TRUE),
    fit_i(d$X, d$y_gauss, eta = 1, y.external = d$ext_gauss, ssr = FALSE)
  )
})

test_that("BRIERi.eta binomial LASSO: ssr TRUE and FALSE give the same path", {
  d <- ssr_individual_case(seed = 2)
  expect_same_path(
    fit_i(d$X, d$y_bin, family = "binomial", ssr = TRUE),
    fit_i(d$X, d$y_bin, family = "binomial", ssr = FALSE)
  )
  expect_same_path(
    fit_i(d$X, d$y_bin, family = "binomial", eta = 1, y.external = d$ext_bin, ssr = TRUE),
    fit_i(d$X, d$y_bin, family = "binomial", eta = 1, y.external = d$ext_bin, ssr = FALSE)
  )
})

test_that("ssr is passed from BRIERi to BRIERi.eta", {
  d <- ssr_individual_case(seed = 3)
  B <- matrix(c(0, rep(0.1, ncol(d$X))), ncol = 1)
  g_ssr <- suppressWarnings(BRIERi(d$X, d$y_gauss, eta.list = list(c(0, 1)),
                                   beta.external = B, nlambda = 30,
                                   parallel = FALSE, ssr = TRUE))
  g_all <- suppressWarnings(BRIERi(d$X, d$y_gauss, eta.list = list(c(0, 1)),
                                   beta.external = B, nlambda = 30,
                                   parallel = FALSE, ssr = FALSE))
  for (k in seq_along(g_ssr$res)) {
    expect_same_path(g_ssr$res[[k]], g_all$res[[k]])
  }
  expect_error(fit_i(d$X, d$y_gauss, ssr = "yes"), "ssr must be")

  # the switch really reaches the fitter: with the strong rule fitter disabled,
  # ssr = FALSE still fits (through BRIERi and, via ..., through BRIERi.cv) and
  # ssr = TRUE does not
  local_mocked_bindings(cd_wgaussian_fit_ssr = function(...) stop("strong rule fitter called"))
  expect_error(
    BRIERi(d$X, d$y_gauss, eta.list = list(0), nlambda = 10, parallel = FALSE),
    "strong rule fitter called"
  )
  expect_s3_class(
    BRIERi(d$X, d$y_gauss, eta.list = list(0), nlambda = 10, parallel = FALSE, ssr = FALSE),
    "BRIER"
  )
  expect_s3_class(
    suppressWarnings(BRIERi.cv(d$X, d$y_gauss, eta.list = list(0), nlambda = 10,
                               nfolds = 3, parallel = FALSE, ssr = FALSE)),
    "BRIER"
  )
})
