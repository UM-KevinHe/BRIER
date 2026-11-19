## The validation deviance of the selection routines is the MEAN per validation observation
## (v1.7.7), the scale evalMetric(), calcDev() and BRIERi.cv() already used. Until v1.7.6
## binomial.dev and poisson.dev were TOTALS over the validation rows, so a selection's
## measure (about n times larger) could not be compared with evalMetric() of the same model
## on the same data, or with an external model scored by evalMetric().

make_val_case <- function(family, n = 600, n.val = 300, p = 30, seed = 1) {
  set.seed(seed)
  b <- c(rep(0.4, 5), rep(0, p - 5))
  draw <- function(m) {
    X <- matrix(rnorm(m * p), m, p)
    lp <- drop(X %*% b)
    y <- switch(family,
      gaussian = lp + rnorm(m),
      binomial = rbinom(m, 1, plogis(lp - 0.5)),
      poisson = rpois(m, exp(0.2 * lp))
    )
    list(X = X, y = y)
  }
  tr <- draw(n)
  va <- draw(n.val)
  fit <- suppressWarnings(suppressMessages(BRIERi(
    tr$X, tr$y, family = family, eta.list = c(0, 1), beta.external = c(0, b),
    parallel = FALSE, ncores = 1
  )))
  list(fit = fit, va = va)
}

same_as_evalMetric <- function(family, criteria) {
  d <- make_val_case(family)
  sel <- suppressMessages(BRIERi.selection(
    d$fit, criteria = criteria, X.val = d$va$X, y.val = d$va$y
  ))
  pred <- predict(sel, X = d$va$X, type = "response")
  expect_equal(
    sel$eta.lambda$measure.min[sel$eta.min.index],
    evalMetric(pred = as.numeric(pred), y = d$va$y, criteria = criteria),
    tolerance = 1e-8
  )
}

test_that("binomial.dev of a selection equals evalMetric of the selected model", {
  same_as_evalMetric("binomial", "binomial.dev")
})

test_that("poisson.dev of a selection equals evalMetric of the selected model", {
  same_as_evalMetric("poisson", "poisson.dev")
})

test_that("gaussian.mspe of a selection equals evalMetric of the selected model", {
  same_as_evalMetric("gaussian", "gaussian.mspe")
})

test_that("the validation deviance does not grow with the validation sample size", {
  d <- make_val_case("binomial", n.val = 300)
  half <- seq_len(150)
  s.all <- suppressMessages(BRIERi.selection(
    d$fit, criteria = "binomial.dev", X.val = d$va$X, y.val = d$va$y
  ))
  s.half <- suppressMessages(BRIERi.selection(
    d$fit, criteria = "binomial.dev", X.val = d$va$X[half, ], y.val = d$va$y[half]
  ))
  ratio <- min(s.all$eta.lambda$measure.min) / min(s.half$eta.lambda$measure.min)
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.4)
})
