## predict() / coef() at a non-selected eta, and subsetEta() (v1.5.0).
##
## Before v1.5.0, predict(sel, X, eta = 0) used lambda.min.index, the best lambda of the
## SELECTED eta, at eta = 0: a different model from the target-only model the selection
## had actually chosen a lambda for. The tests check three things, in this order:
##   1. every call that worked before gives the same answer (backward compatibility);
##   2. a non-selected eta now gets its own best lambda;
##   3. subsetEta() restricts fits and selections correctly.

make_sel_case <- function(seed = 7, n = 600, p = 30) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  b <- c(rep(0.15, 8), rep(0, p - 8))
  y <- drop(X %*% b) + rnorm(n)
  Xv <- matrix(rnorm(300 * p), 300, p)
  yv <- drop(Xv %*% b) + rnorm(300)
  be <- c(0, b + rnorm(p, sd = 0.03))
  fit <- suppressWarnings(suppressMessages(BRIERi(
    X, y, family = "gaussian", beta.external = be,
    eta.list = c(0, 0.1, 0.5, 1, 3.5938, 10), parallel = FALSE, ncores = 1
  )))
  sel <- suppressMessages(BRIERi.selection(fit, criteria = "gaussian.mspe", X.val = Xv, y.val = yv))
  list(X = X, y = y, Xv = Xv, yv = yv, fit = fit, sel = sel)
}

eta_pred <- function(obj, i, j, X) {
  predict.BRIER.eta(obj$res[[i]], X = X, which = j)
}

## ---- 1. backward compatibility ----------------------------------------------------

test_that("the default prediction is still the selected eta at its selected lambda", {
  d <- make_sel_case()
  s <- d$sel
  expect_equal(predict(s, d$Xv), eta_pred(s, s$eta.min.index, s$lambda.min.index, d$Xv))
  expect_equal(coef(s), coef.BRIER.eta(s$res[[s$eta.min.index]], which = s$lambda.min.index))
})

test_that("an explicit which.lambda is always honoured, at any eta", {
  d <- make_sel_case()
  s <- d$sel
  for (i in seq_len(nrow(s$eta.grid))) {
    expect_equal(predict(s, d$Xv, which.eta = i, which.lambda = 5), eta_pred(s, i, 5, d$Xv))
    expect_equal(coef(s, which.eta = i, which.lambda = 5),
                 coef.BRIER.eta(s$res[[i]], which = 5))
  }
})

test_that("the selected eta given explicitly is unchanged (lambda.min.index)", {
  d <- make_sel_case()
  s <- d$sel
  expect_equal(predict(s, d$Xv, which.eta = s$eta.min.index), predict(s, d$Xv))
  expect_equal(predict(s, d$Xv, eta = s$eta.min), predict(s, d$Xv))
})

test_that("every selection sets lambda.min.index from eta.lambda (the invariant the fix relies on)", {
  d <- make_sel_case()
  for (s in list(
    d$sel,
    suppressMessages(BRIERi.selection(d$fit, criteria = "Cp")),
    suppressMessages(BRIERi.selection(d$fit, criteria = "gaussian.rsq", X.val = d$Xv, y.val = d$yv))
  )) {
    el <- s$eta.lambda
    expect_equal(s$lambda.min.index, el$lambda.min.index[el$eta.index == s$eta.min.index])
  }
})

test_that("multiple etas still use each eta's best lambda (unchanged branch)", {
  d <- make_sel_case()
  s <- d$sel
  el <- s$eta.lambda
  got <- predict(s, d$Xv, which.eta = c(1, 3))
  expect_equal(unname(got[[1]]), eta_pred(s, 1, el$lambda.min.index[1], d$Xv))
  expect_equal(unname(got[[2]]), eta_pred(s, 3, el$lambda.min.index[3], d$Xv))
})

test_that("an untuned fit still needs a lambda", {
  d <- make_sel_case()
  expect_error(predict(d$fit, d$Xv, which.eta = 2), "No lambda selected")
  expect_equal(predict(d$fit, d$Xv, which.eta = 2, which.lambda = 4), eta_pred(d$fit, 2, 4, d$Xv))
})

## ---- 2. the fix --------------------------------------------------------------------

test_that("a non-selected eta is predicted at ITS OWN best lambda", {
  d <- make_sel_case()
  s <- d$sel
  el <- s$eta.lambda
  for (i in seq_len(nrow(s$eta.grid))) {
    want <- eta_pred(s, i, el$lambda.min.index[el$eta.index == i], d$Xv)
    expect_equal(predict(s, d$Xv, which.eta = i), want)
    expect_equal(predict(s, d$Xv, eta = s$eta.grid[i, 1]), want)
    expect_equal(coef(s, which.eta = i),
                 coef.BRIER.eta(s$res[[i]], which = el$lambda.min.index[el$eta.index == i]))
  }
})

test_that("the eta = 0 prediction is the one predict() over the grid gives for eta = 0", {
  d <- make_sel_case()
  s <- d$sel
  multi <- predict(s, d$Xv, eta = s$eta.grid)
  expect_equal(predict(s, d$Xv, eta = 0), unname(multi[[1]]))
})

test_that("typed eta values match the nearest fitted eta within 1% only", {
  d <- make_sel_case()
  s <- d$sel
  expect_equal(predict(s, d$Xv, eta = 3.59), predict(s, d$Xv, which.eta = 5))
  expect_error(predict(s, d$Xv, eta = 7), "does not match any row")
  expect_error(predict(s, d$Xv, eta = 0.005), "does not match any row")
  expect_error(predict(s, d$Xv, eta = 0.11), "does not match any row")
})

## ---- 3. subsetEta() --------------------------------------------------------------

test_that("one eta on a selection: pinned to it, with its own best lambda", {
  d <- make_sel_case()
  s <- d$sel
  t0 <- subsetEta(s, eta = 0)
  expect_s3_class(t0, "BRIER.selection")
  expect_equal(nrow(t0$eta.grid), 1)
  expect_equal(unname(t0$eta.min), 0)
  expect_equal(t0$lambda.min.index, s$eta.lambda$lambda.min.index[1])
  expect_equal(predict(t0, d$Xv), predict(s, d$Xv, eta = 0))
  expect_equal(t0$eta.subset, 1L)
})

test_that("a vector of etas: the best (eta, lambda) WITHIN them", {
  d <- make_sel_case()
  s <- d$sel
  keep <- c(1, 2, 3)
  sub <- subsetEta(s, which.eta = keep)
  best <- keep[which.min(s$eta.lambda$measure.min[keep])]
  expect_equal(sub$eta.subset[sub$eta.min.index], best)
  expect_equal(predict(sub, d$Xv), predict(s, d$Xv, which.eta = best))
  expect_equal(subsetEta(s, eta = s$eta.grid[keep, 1])$eta.subset, keep)
})

test_that("keeping every eta reproduces the selection", {
  d <- make_sel_case()
  s <- d$sel
  all <- subsetEta(s, which.eta = seq_len(nrow(s$eta.grid)))
  expect_equal(predict(all, d$Xv), predict(s, d$Xv))
  expect_equal(coef(all), coef(s))
  expect_equal(all$eta.min.index, s$eta.min.index)
})

test_that("an untuned fit is subset and can then be selected as usual", {
  d <- make_sel_case()
  sub <- subsetEta(d$fit, eta = c(0, 1, 10))
  expect_s3_class(sub, "BRIER")
  expect_false(inherits(sub, "BRIER.selection"))
  expect_null(sub$eta.lambda)
  expect_equal(length(sub$res), 3)
  expect_equal(sub$eta.list, c(0, 1, 10))
  s2 <- suppressMessages(BRIERi.selection(sub, criteria = "gaussian.mspe", X.val = d$Xv, y.val = d$yv))
  full <- d$sel$eta.lambda[d$sel$eta.lambda$eta.index %in% c(1, 4, 6), ]
  expect_equal(s2$eta.lambda$measure.min, full$measure.min)
  expect_equal(s2$eta.lambda$lambda.min.index, full$lambda.min.index)
})

test_that("subsetEta refuses what was not fitted and malformed requests", {
  d <- make_sel_case()
  expect_error(subsetEta(d$sel, eta = 2), "does not match any row")
  expect_error(subsetEta(d$sel, eta = 0, which.eta = 1), "not both")
  expect_error(subsetEta(d$sel), "Give the eta values")
  expect_error(subsetEta(d$sel, which.eta = 99), "between 1 and")
  expect_error(subsetEta(list(), eta = 0), "class 'BRIER'")
})

test_that("subsetEta works on BRIERs and on BRIERi.cv selections", {
  d <- make_sel_case()
  n <- nrow(d$X)
  Xs <- scale(d$X)
  XtX <- Matrix::Matrix(crossprod(Xs) / (n - 1), sparse = TRUE)
  ss <- data.frame(corr = drop(crossprod(Xs, scale(d$y))) / (n - 1))
  fs <- suppressWarnings(suppressMessages(BRIERs(
    ss, XtX, family = "gaussian", beta.external = d$fit$beta.external[-1],
    eta.list = c(0, 1, 10), parallel = FALSE, ncores = 1
  )))
  ssel <- suppressMessages(BRIERs.selection(fs, criteria = "gaussian.mspe",
                                            X.val = scale(d$Xv), y.val = as.numeric(scale(d$yv))))
  t0 <- subsetEta(ssel, eta = 0)
  expect_equal(predict(t0, scale(d$Xv)), predict(ssel, scale(d$Xv), eta = 0))
  expect_equal(t0$lambda.min.index, ssel$eta.lambda$lambda.min.index[1])

  cv <- suppressWarnings(suppressMessages(BRIERi.cv(
    d$X, d$y, family = "gaussian", beta.external = d$fit$beta.external,
    eta.list = c(0, 1), nfolds = 3, parallel = FALSE, ncores = 1
  )))
  c0 <- subsetEta(cv, eta = 0)
  expect_s3_class(c0, "BRIER.cv")
  expect_equal(c0$lambda.min.index, cv$eta.lambda$lambda.min.index[1])
  expect_equal(predict(c0, d$Xv), predict(cv, d$Xv, eta = 0))
})
