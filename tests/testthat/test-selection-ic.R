## The information criteria of BRIERi.selection, on the scale of the deviance they
## penalise (v1.4.4), and for a Gaussian outcome free of the units of y (v1.4.5). BRIERi.eta stores a MEAN deviance (the weights sum to one), so the penalty
## is divided by n. Before v1.4.4 it was not, and every criterion chose the empty model
## on data with real signal.

make_ic_case <- function(family = "gaussian", n = 2000, p = 50, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  b <- c(rep(if (family == "gaussian") 0.1 else 0.25, 10), rep(0, p - 10))
  lp <- drop(X %*% b)
  y <- if (family == "gaussian") lp + rnorm(n) else rbinom(n, 1, plogis(lp))
  fit <- suppressWarnings(suppressMessages(BRIERi(
    X, y, family = family, eta.list = 0, beta.external = rep(0, p + 1),
    parallel = FALSE, ncores = 1
  )))
  list(fit = fit, n = n)
}

nonzero_at <- function(sel) {
  sum(sel$res[[sel$eta.min.index]]$beta[-1, sel$lambda.min.index] != 0)
}

test_that("BIC on a Gaussian outcome recovers real signal, not the empty model", {
  d <- make_ic_case("gaussian")
  s <- suppressMessages(BRIERi.selection(d$fit, criteria = "BIC"))
  expect_gte(nonzero_at(s), 5)
  expect_lte(nonzero_at(s), 25)
})

test_that("Gaussian BIC / AIC are the profile-likelihood forms, n log(RSS/n) + kappa df", {
  d <- make_ic_case("gaussian")
  r <- d$fit$res[[1]]
  df <- r$k / (1 + sum(r$eta))
  rss <- d$n * r$deviance
  s <- suppressMessages(BRIERi.selection(d$fit, criteria = "BIC"))
  expect_equal(as.numeric(s$lambda.min.index),
               as.numeric(which.min(d$n * log(rss / d$n) + log(d$n) * df)))
  a <- suppressMessages(BRIERi.selection(d$fit, criteria = "AIC"))
  expect_equal(as.numeric(a$lambda.min.index),
               as.numeric(which.min(d$n * log(rss / d$n) + 2 * df)))
})

test_that("Gaussian BIC does not move with the units of y", {
  # v1.4.4 penalised RSS / n directly: y x 10 kept 48 of 50, y x 0.1 the empty model.
  set.seed(1)
  n <- 2000; p <- 50
  X <- matrix(rnorm(n * p), n, p)
  y <- drop(X %*% c(rep(0.1, 10), rep(0, p - 10))) + rnorm(n)
  pick <- function(yy) {
    f <- suppressWarnings(suppressMessages(BRIERi(
      X, yy, eta.list = 0, beta.external = rep(0, p + 1), parallel = FALSE, ncores = 1)))
    nonzero_at(suppressMessages(BRIERi.selection(f, criteria = "BIC")))
  }
  k1 <- pick(y)
  expect_identical(pick(0.1 * y), k1)
  expect_identical(pick(10 * y), k1)
  expect_gte(k1, 5)
})

test_that("binary BIC is the total-deviance BIC (dispersion 1)", {
  d <- make_ic_case("binomial")
  r <- d$fit$res[[1]]
  df <- r$k / (1 + sum(r$eta))
  s <- suppressMessages(BRIERi.selection(d$fit, criteria = "BIC"))
  expect_equal(as.numeric(s$lambda.min.index),
               as.numeric(which.min(d$n * r$deviance + log(d$n) * df)))
})

test_that("AIC penalises less than BIC, so it keeps at least as many effects", {
  d <- make_ic_case("gaussian")
  a <- suppressMessages(BRIERi.selection(d$fit, criteria = "AIC"))
  b <- suppressMessages(BRIERi.selection(d$fit, criteria = "BIC"))
  expect_gte(nonzero_at(a), nonzero_at(b))
})

test_that("BIC on a binary outcome recovers real signal too", {
  d <- make_ic_case("binomial")
  s <- suppressMessages(BRIERi.selection(d$fit, criteria = "BIC"))
  expect_gte(nonzero_at(s), 3)
})
