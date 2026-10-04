## df.method: the degrees of freedom the information criteria use (R/df.R).

df_case <- function(n = 300, p = 12, seed = 3) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  b <- c(1, -0.8, 0.6, 0.4, rep(0, p - 4))
  list(X = X, y = drop(X %*% b) + rnorm(n), n = n, p = p)
}

fit_eta0 <- function(d, y = d$y, ...) {
  suppressWarnings(suppressMessages(BRIERi.eta(d$X, y, eta = 0, ...)))
}

test_that("the default is the active count, and nothing extra is stored", {
  d <- df_case()
  f <- fit_eta0(d, nlambda = 20)
  expect_identical(f$df.method, "active")
  expect_null(f$df.eff)
})

test_that("for the LASSO the divergence df IS the active count", {
  d <- df_case()
  f <- fit_eta0(d, nlambda = 20, df.method = "divergence")
  expect_equal(f$df.eff, as.numeric(f$k), tolerance = 1e-8)
  fb <- suppressWarnings(suppressMessages(BRIERi.eta(
    d$X, as.numeric(d$y > 0), eta = 0, family = "binomial", nlambda = 20,
    df.method = "divergence")))
  expect_equal(fb$df.eff, as.numeric(fb$k), tolerance = 1e-8)
})

test_that("the LQA df lies below the active count for the LASSO", {
  d <- df_case()
  f <- fit_eta0(d, nlambda = 20, df.method = "lqa")
  act <- f$k > 1
  expect_true(all(f$df.eff[act] < f$k[act]))
})

## THE CHECK THAT MATTERS: the divergence df equals the Stein divergence of the fit,
## sum_i d yhat_i / d y_i, measured by perturbing each y_i and refitting at the same
## lambda. It tests the whole formula, including the lambda scaling, the MCP curvature
## and the ridge term.
stein_df <- function(d, lam, ...) {
  fit_at <- function(y) {
    f <- fit_eta0(d, y = y, lambda = lam, ...)
    drop(cbind(1, d$X) %*% f$beta[, 1])
  }
  h <- 1e-5
  y0 <- d$y; f0 <- fit_at(y0)
  sum(vapply(seq_len(d$n), function(i) {
    y1 <- y0; y1[i] <- y1[i] + h
    (fit_at(y1)[i] - f0[i]) / h
  }, 0))
}

test_that("MCP: the divergence df matches the measured Stein divergence", {
  d <- df_case(n = 150, p = 8)
  lam <- 0.08
  f <- fit_eta0(d, lambda = lam, penalty = "MCP", gamma = 3, df.method = "divergence",
                eps = 1e-10)
  expect_equal(f$df.eff, stein_df(d, lam, penalty = "MCP", gamma = 3, eps = 1e-10),
               tolerance = 1e-3)
  expect_gt(f$df.eff, f$k)    # nonconvex: more df than the active count
})

test_that("elastic net: the ridge term enters the divergence df", {
  d <- df_case(n = 150, p = 8)
  lam <- 0.08
  f <- fit_eta0(d, lambda = lam, alpha = 0.5, df.method = "divergence", eps = 1e-10)
  expect_equal(f$df.eff, stein_df(d, lam, alpha = 0.5, eps = 1e-10), tolerance = 1e-3)
  expect_lt(f$df.eff, f$k)    # ridge shrinkage: fewer df than the active count
})

test_that("the selection uses the stored df", {
  d <- df_case()
  fit <- suppressWarnings(suppressMessages(BRIERi(
    d$X, d$y, eta.list = 0, beta.external = rep(0, d$p + 1), penalty = "MCP",
    df.method = "divergence", parallel = FALSE, ncores = 1)))
  r <- fit$res[[1]]
  expect_false(is.null(r$df.eff))
  s <- suppressMessages(BRIERi.selection(fit, criteria = "BIC"))
  expect_equal(as.numeric(s$lambda.min.index),
               as.numeric(which.min(log(r$deviance) + log(d$n) * r$df.eff / d$n)))
})

test_that("BRIERs: divergence = the active count for the LASSO; MCP exceeds it", {
  d <- df_case()
  ys <- as.numeric(scale(d$y)); Xs <- scale(d$X)
  R <- Matrix::Matrix(crossprod(Xs) / (d$n - 1), sparse = TRUE)
  r <- as.numeric(crossprod(Xs, ys)) / (d$n - 1)
  fl <- suppressWarnings(suppressMessages(BRIERs.eta(R, r, nlambda = 20,
                                                    df.method = "divergence")))
  expect_equal(fl$df.eff, as.numeric(fl$k), tolerance = 1e-8)
  fm <- suppressWarnings(suppressMessages(BRIERs.eta(R, r, nlambda = 20, penalty = "MCP",
                                                    df.method = "divergence")))
  act <- fm$k > 0
  expect_true(all(fm$df.eff[act] >= fm$k[act] - 1e-8))
})

test_that("the sparse LD trace is solved block by block and equals the dense trace", {
  set.seed(7)
  blk <- function(k) { A <- matrix(rnorm(k * 40), 40, k); stats::cor(A) }
  R <- Matrix::bdiag(blk(5), blk(8), blk(3), blk(6))
  R <- methods::as(methods::as(R, "generalMatrix"), "CsparseMatrix")
  curv <- runif(ncol(R), 0, 0.3)
  expect_identical(as.integer(.block_ends(R)), c(5L, 13L, 16L, 22L))
  expect_equal(.df_trace(R, curv), .df_trace_dense(as.matrix(R), curv), tolerance = 1e-10)
  Rs <- Matrix::forceSymmetric(R)            # upper-triangle storage reads the same
  expect_identical(as.integer(.block_ends(Rs)), c(5L, 13L, 16L, 22L))
})
