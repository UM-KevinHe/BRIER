# =============================================================================
# screenExternals: the cascade, and the property that makes it a cascade.
#
# The rules under test come from the ten-trait application, so the tests assert
# the SHAPE of the rule (exclusivity, ordering, the clips) on panels built so
# the right answer is known, rather than re-deriving the application's numbers.
# =============================================================================

# A panel where the truth is known by construction: three models carrying real
# signal, one pure noise, one NEAR-copy of the first at the same scale, one
# RESCALED copy of the first, one all zero.
#
# The near-copy and the rescaled copy are different tests and both are needed.
# dedupExternals measures cosine similarity, which is scale invariant, so it
# calls BOTH duplicates of signal1. The SCREEN does not: its calibration leg
# d2 = 2b'r - b'Rb is scale SENSITIVE, so a 2.5x rescaling overshoots and is
# rejected before dedup ever sees it. Only the same-scale near-copy reaches
# de-duplication, which is what makes it the right fixture for that test.
make_panel <- function(p = 60, seed = 11) {
  set.seed(seed)
  b1 <- c(rnorm(10, 0, 0.3), rep(0, p - 10))
  b2 <- c(rep(0, 5), rnorm(10, 0, 0.3), rep(0, p - 15))
  b3 <- c(rnorm(6, 0, 0.25), rep(0, p - 6))
  noise <- c(rep(0, p - 8), rnorm(8, 0, 0.02))
  nearcopy <- b1 + c(rnorm(10, 0, 0.004), rep(0, p - 10))
  cbind(signal1 = b1, signal2 = b2, signal3 = b3, noise = noise,
        nearcopy = nearcopy, rescaled = 2.5 * b1, empty = rep(0, p))
}

make_summary <- function(p = 60, seed = 11) {
  B <- make_panel(p, seed)
  set.seed(seed + 1)
  n <- 400
  X <- matrix(rnorm(n * p), n, p)
  R <- crossprod(X) / n
  y <- drop(X %*% B[, "signal1"]) + rnorm(n, 0, 0.5)
  r <- as.numeric(crossprod(X, y)) / n
  r <- r / stats::sd(y)
  list(B = B, XtX = R, corr = r, n = 5000)
}

test_that("the two screens are exclusive, and neither silently wins", {
  d <- make_summary()
  expect_error(
    screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = d$n,
                    fit = structure(list(XtY = 1), class = "BRIER"),
                    X.val = matrix(0, 4, 60), y.val = rep(0, 4)),
    "alternatives, not stages"
  )
  expect_error(screenExternals(d$B), "No screen can run")
})

test_that("the validation-free screen removes what carries nothing", {
  d <- make_summary()
  s <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = d$n,
                       family = "gaussian", dedup.cor = NULL)
  e <- s$externals
  expect_identical(s$screen, "validation-free")
  # an all-zero model is degenerate before any screen runs
  expect_true(e$degenerate[e$label == "empty"])
  expect_false(e$keep[e$label == "empty"])
  # the model the outcome was built from must survive
  expect_true(e$keep[e$label == "signal1"])
  # both legs are reported, and keep is exactly their conjunction
  expect_true(all(c("d2", "r2", "inner_product", "prs_variance") %in% names(e)))
  expect_equal(e$keep,
               !e$degenerate & is.finite(e$d2) & e$d2 > 0 &
                 is.finite(e$r2) & e$r2 > s$floor)
})

test_that("the no-signal floor moves with n, and is 0.5 for a binomial outcome", {
  d <- make_summary()
  lo <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = 100,
                        family = "gaussian", dedup.cor = NULL)
  hi <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = 1e6,
                        family = "gaussian", dedup.cor = NULL)
  expect_gt(lo$floor, hi$floor)
  expect_gte(sum(lo$externals$keep), sum(hi$externals$keep))
  bn <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = 5000,
                        family = "binomial", dedup.cor = NULL)
  expect_equal(bn$floor, 0.5)
})

test_that("the validation-free screen refuses to guess n", {
  d <- make_summary()
  expect_error(screenExternals(d$B, corr = d$corr, XtX = d$XtX), "sample size")
})

test_that("the calibration leg is scale sensitive, and rejects a rescaled copy", {
  # d2 = 2b'r - b'Rb is measured at the PUBLISHED scale, so a score that is
  # right in shape but wrong in scale fails it even though its scale-free R^2
  # is identical to the original's. That is the chosen rule, not an accident:
  # the alternative is a sign-only leg that cannot see mis-scaling at all.
  d <- make_summary()
  s <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = d$n,
                       family = "gaussian", dedup.cor = NULL)
  e <- s$externals
  sig <- e[e$label == "signal1", ]; res <- e[e$label == "rescaled", ]
  expect_equal(res$r2, sig$r2, tolerance = 1e-8)   # scale free: identical
  expect_lt(res$d2, sig$d2)                        # at scale: worse
  expect_false(res$keep)
})

test_that("de-duplication runs in RANK order, so the better copy survives", {
  # nearcopy is signal1 plus a tiny same-scale perturbation, so cosine
  # similarity is ~1 and both clear the screen. Exactly one may survive dedup,
  # and WHICH one is the point: it must be the higher-ranked of the pair, not
  # whichever came first in the panel.
  d <- make_summary()
  s <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = d$n,
                       family = "gaussian", dedup.cor = 0.9)
  e <- s$externals
  expect_true(e$keep[e$label == "signal1"])
  expect_true(e$keep[e$label == "nearcopy"])       # both reach dedup
  kept <- e$label[e$dedup]
  expect_false(all(c("signal1", "nearcopy") %in% kept))
  expect_gt(s$n.deduplicated, 0L)

  both <- e[e$label %in% c("signal1", "nearcopy"), ]
  winner <- both$label[which.max(both$r2)]
  expect_true(winner %in% kept)
})

test_that("rank is assigned over the final set, in order, best first", {
  d <- make_summary()
  s <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = d$n,
                       family = "gaussian")
  e <- s$externals
  expect_identical(which(!is.na(e$rank)), sort(s$kept))
  expect_identical(e$rank[s$kept], seq_along(s$kept))
  expect_identical(e$dedup, !is.na(e$rank))
  # ordered by the screen's own statistic, decreasing
  expect_false(is.unsorted(rev(e$r2[s$kept])))
  expect_identical(s$rank.statistic, "r2")
})

test_that("top.k cuts the ranked survivors and nothing else", {
  d <- make_summary()
  full <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = d$n,
                          family = "gaussian")
  cut <- screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = d$n,
                         family = "gaussian", top.k = 2)
  expect_lte(length(cut$kept), 2L)
  expect_identical(cut$kept, utils::head(full$kept, length(cut$kept)))
  # the screen itself is untouched: top.k is a study-design choice, not a screen
  expect_identical(cut$survivors, full$survivors)
})

test_that("an empty panel and a labels mismatch are both refused", {
  d <- make_summary()
  expect_error(screenExternals(d$B[, 0, drop = FALSE], corr = d$corr,
                               XtX = d$XtX, n = d$n), "empty")
  expect_error(screenExternals(d$B, corr = d$corr, XtX = d$XtX, n = d$n,
                               labels = c("a", "b")), "length 2")
  expect_identical(ncol(d$B), 7L)
})


# -- the cross-validated screen -----------------------------------------------

make_cv_case <- function(p = 40, n = 300, seed = 21) {
  set.seed(seed)
  B <- cbind(good = c(rnorm(8, 0, 0.4), rep(0, p - 8)),
             bad  = c(rep(0, p - 6), rnorm(6, 0, 0.5)),
             zero = rep(0, p))
  X <- matrix(rnorm(n * p), n, p)
  y <- drop(X %*% B[, "good"]) + rnorm(n, 0, 0.4)
  y <- as.numeric(scale(y))
  Xv <- matrix(rnorm(n * p), n, p)
  yv <- drop(Xv %*% B[, "good"]) + rnorm(n, 0, 0.4)
  yv <- as.numeric(scale(yv))
  list(B = B, X = X, y = y, X.val = Xv, y.val = yv, p = p)
}

test_that("the cross-validated screen runs, reports both legs, and ranks on loss", {
  d <- make_cv_case()
  base <- suppressWarnings(suppressMessages(BRIERi(
    d$X, d$y, family = "gaussian", eta.list = 0,
    beta.external = rbind(0, d$B[, "zero", drop = FALSE]),
    nlambda = 10, parallel = FALSE, ncores = 1
  )))
  s <- suppressWarnings(suppressMessages(screenExternals(
    d$B, fit = base, X.val = d$X.val, y.val = d$y.val,
    family = "gaussian", nfolds = 3, seed = 1, dedup.cor = NULL
  )))
  expect_identical(s$screen, "cv")
  expect_identical(s$nfolds, 3L)
  expect_identical(s$rank.statistic, "ratio_loss")
  expect_true(all(c("ratio_loss", "ratio_acc") %in% names(s$externals)))
  # the model the outcome was built from beats the noise model on the loss leg
  e <- s$externals
  expect_gt(e$ratio_loss[e$label == "good"], e$ratio_loss[e$label == "bad"])
  # an all-zero model is degenerate here too
  expect_true(e$degenerate[e$label == "zero"])
})

test_that("keep is the conjunction of both legs at the bar", {
  d <- make_cv_case()
  base <- suppressWarnings(suppressMessages(BRIERi(
    d$X, d$y, family = "gaussian", eta.list = 0,
    beta.external = rbind(0, d$B[, "zero", drop = FALSE]),
    nlambda = 10, parallel = FALSE, ncores = 1
  )))
  s <- suppressWarnings(suppressMessages(screenExternals(
    d$B, fit = base, X.val = d$X.val, y.val = d$y.val,
    family = "gaussian", nfolds = 3, seed = 1, dedup.cor = NULL
  )))
  e <- s$externals
  expect_equal(e$keep,
               !e$degenerate &
                 is.finite(e$ratio_loss) & e$ratio_loss >= 1 &
                 is.finite(e$ratio_acc) & e$ratio_acc >= 1)
  # raising the bar can only remove models, never add them
  hard <- suppressWarnings(suppressMessages(screenExternals(
    d$B, fit = base, X.val = d$X.val, y.val = d$y.val,
    family = "gaussian", nfolds = 3, seed = 1, bar = 1.5, dedup.cor = NULL
  )))
  expect_true(all(which(hard$externals$keep) %in% which(e$keep)))
})

test_that("nfolds is an argument and 3 is only its default", {
  d <- make_cv_case()
  base <- suppressWarnings(suppressMessages(BRIERi(
    d$X, d$y, family = "gaussian", eta.list = 0,
    beta.external = rbind(0, d$B[, "zero", drop = FALSE]),
    nlambda = 10, parallel = FALSE, ncores = 1
  )))
  five <- suppressWarnings(suppressMessages(screenExternals(
    d$B, fit = base, X.val = d$X.val, y.val = d$y.val,
    family = "gaussian", nfolds = 5, seed = 1, dedup.cor = NULL
  )))
  expect_identical(five$nfolds, 5L)
  expect_error(
    screenExternals(d$B, fit = base, X.val = d$X.val, y.val = d$y.val,
                    family = "gaussian", nfolds = 1),
    "at least 2"
  )
  expect_error(
    screenExternals(d$B, fit = base, X.val = d$X.val, y.val = d$y.val,
                    family = "gaussian", nfolds = nrow(d$X.val) + 1),
    "samples"
  )
})
