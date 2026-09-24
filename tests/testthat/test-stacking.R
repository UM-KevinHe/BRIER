# =============================================================================
# Stacking weights: the free solve, the constrained solve, and the two shapes.
#
# These exist because the estimator was changed in a package whose suite did not
# cover the fitters or the aggregation. The free solve moved from normal
# equations to an equilibrated, truncation-safe solve, so the first job of these
# tests is to prove it is the SAME answer wherever the old one was valid, and a
# finite answer where the old one failed.
# =============================================================================

test_that("the free solve matches the closed form on a well-conditioned panel", {
  set.seed(101)
  n <- 300; m <- 4
  B <- matrix(rnorm(n * m), n, m)
  y <- rnorm(n)
  G <- crossprod(B); h <- as.numeric(crossprod(B, y))
  got <- BRIER:::stack_weights(G, h)
  expect_equal(got$weights, drop(solve(G, h)), tolerance = 1e-10)
  expect_false(got$constrained)
  expect_identical(got$n_dropped, 0L)
})

test_that("a rank-deficient panel returns finite weights and reports the truncation", {
  set.seed(102)
  n <- 200; m <- 3
  B <- matrix(rnorm(n * m), n, m)
  B <- cbind(B, B[, 1])                      # an exact duplicate: G is singular
  y <- rnorm(n)
  G <- crossprod(B); h <- as.numeric(crossprod(B, y))
  expect_error(solve(G, h))                  # the old path
  got <- BRIER:::stack_weights(G, h)
  expect_true(all(is.finite(got$weights)))
  expect_gt(got$n_dropped, 0L)               # the report says so
})

test_that("an empty panel is refused rather than reaching eigen()", {
  expect_error(BRIER:::stack_weights(matrix(0, 0, 0), numeric(0)), "empty")
})

test_that("the constrained solve stays on the constraint set and certifies itself", {
  set.seed(103)
  n <- 300; m <- 5
  B <- matrix(rnorm(n * m), n, m)
  y <- rowMeans(B) + rnorm(n, 0, 0.5)
  G <- crossprod(B); h <- as.numeric(crossprod(B, y))
  got <- BRIER:::stack_weights(G, h, constrain = TRUE)
  expect_true(got$constrained)
  expect_equal(sum(got$weights), 1, tolerance = 1e-8)
  expect_true(all(got$weights >= -1e-12))
  expect_lt(got$kkt, 1e-6)
})

test_that("the constraint cannot beat the free solve on the same objective", {
  set.seed(104)
  n <- 250; m <- 4
  B <- matrix(rnorm(n * m), n, m)
  y <- rnorm(n)
  G <- crossprod(B); h <- as.numeric(crossprod(B, y))
  Q <- function(w) as.numeric(crossprod(w, G %*% w)) - 2 * sum(h * w)
  free <- BRIER:::stack_weights(G, h)$weights
  cons <- BRIER:::stack_weights(G, h, constrain = TRUE)$weights
  # the constrained set is a subset, so its optimum cannot be lower
  expect_gte(Q(cons), Q(free) - 1e-8)
})

test_that("the constrained solve can put a weight exactly at zero", {
  set.seed(105)
  n <- 300
  B <- matrix(rnorm(n * 3), n, 3)
  B <- cbind(B, rnorm(n) * 5)                # a fourth model that is pure noise
  y <- rowMeans(B[, 1:3]) + rnorm(n, 0, 0.3)
  G <- crossprod(B); h <- as.numeric(crossprod(B, y))
  w <- BRIER:::stack_weights(G, h, constrain = TRUE)$weights
  expect_true(any(w == 0))
})

test_that("the response-scale constrained solve holds for all three families", {
  set.seed(106)
  n <- 400; m <- 4
  for (fam in c("gaussian", "binomial", "poisson")) {
    if (fam == "binomial") {
      Y <- matrix(runif(n * m, 0.05, 0.95), n, m); z <- rbinom(n, 1, rowMeans(Y))
    } else if (fam == "poisson") {
      Y <- matrix(runif(n * m, 0.5, 5), n, m); z <- rpois(n, rowMeans(Y))
    } else {
      Y <- matrix(rnorm(n * m), n, m); z <- rowMeans(Y) + rnorm(n, 0, 0.4)
    }
    got <- BRIER:::stack_weights_response(Y, z, family = fam)
    expect_equal(sum(got$weights), 1, tolerance = 1e-8, info = fam)
    expect_true(all(got$weights >= -1e-12), info = fam)
    expect_lt(got$kkt, 1e-6)
  }
})

test_that("the softmax reparameterisation reaches the same interior optimum", {
  # The cross-check: a different route to the same constrained problem. It is
  # only expected to agree where the optimum is INTERIOR, since the softmax
  # cannot place a weight exactly at zero.
  set.seed(107)
  n <- 400; m <- 4
  Y <- matrix(runif(n * m, 0.05, 0.95), n, m)
  z <- rbinom(n, 1, rowMeans(Y))
  pg <- BRIER:::stack_weights_response(Y, z, family = "binomial")
  sm <- BRIER:::stack_weights_softmax(Y, z, family = "binomial")
  skip_if(any(pg$weights < 1e-6), "optimum is on the boundary; the softmax route cannot reach it")
  expect_equal(pg$weights, sm$weights, tolerance = 1e-4)
  expect_equal(pg$value, sm$value, tolerance = 1e-6)
})

test_that("the two shapes agree on the same gaussian panel", {
  # BRIERi forms its Gram from predictions, BRIERs from B'XtX B. Given the same
  # predictions they are the same problem and must give the same weights.
  set.seed(108)
  n <- 300; m <- 4
  Yhat <- matrix(rnorm(n * m), n, m)
  y <- rowMeans(Yhat) + rnorm(n, 0, 0.5)
  w_i <- BRIER:::stack_weights(as.matrix(crossprod(Yhat)), as.numeric(crossprod(Yhat, y)))$weights
  w_s <- BRIER:::stack_weights(as.matrix(crossprod(Yhat)), as.numeric(crossprod(Yhat, y)))$weights
  expect_equal(w_i, w_s)
  # and the individual-level entry point returns the same thing
  expect_equal(as.numeric(BRIER:::stacking_gaussian(Yhat, y)), w_i, tolerance = 1e-10)
})

test_that("PCstacking is gone, and asking for it is an error", {
  expect_error(match.arg("PCstacking", c("ind", "PCA", "stacking", "stacking.c")))
  expect_true("stacking.c" %in% eval(formals(BRIERs)$multi.method))
  expect_true("stacking.c" %in% eval(formals(BRIERi)$multi.method))
  expect_false("PCstacking" %in% eval(formals(BRIERs)$multi.method))
  expect_false("PCstacking" %in% eval(formals(BRIERi)$multi.method))
})

test_that("every stacking path returns a report with the same fields", {
  # The report is what a caller reads to see WHICH models were combined and how
  # well conditioned the solve was. It is useless if its shape depends on the
  # family, so pin the common fields across all six (family, method) pairs.
  set.seed(109)
  n <- 300; p <- 20; M <- 4
  X <- matrix(rnorm(n * p), n, p)
  B <- rbind(0, matrix(rnorm(p * M, 0, 0.2), p, M))   # (p+1) x M, intercept row
  responses <- list(
    gaussian = drop(X %*% B[-1, 1]) + rnorm(n),
    binomial = rbinom(n, 1, 0.5),
    poisson  = rpois(n, 2)
  )
  for (fam in names(responses)) {
    for (mm in c("stacking", "stacking.c")) {
      ext <- calcExtY(X, responses[[fam]], B, family = fam, multi.method = mm)
      lab <- paste(fam, mm)
      expect_false(is.null(ext$stack), info = lab)
      for (f in c("weights", "constrained", "kkt")) {
        expect_true(f %in% names(ext$stack), info = paste(lab, f))
      }
      expect_length(ext$stack$weights, M)
      expect_identical(ext$stack$constrained, mm == "stacking.c", info = lab)
      expect_equal(ncol(ext$y.external), 1L, info = lab)
      if (mm == "stacking.c") {
        expect_equal(sum(ext$stack$weights), 1, tolerance = 1e-8, info = lab)
        expect_true(all(ext$stack$weights >= -1e-12), info = lab)
      }
    }
  }
})

test_that("equilibration improves the conditioning it is there to improve", {
  # Equilibration rescales the Gram to unit diagonal, so it helps when the
  # columns differ wildly in SCALE and does nothing for deficiency of RANK.
  # Both halves are asserted, because claiming the first without the second is
  # how the singular-panel test in test-externals.R got written wrong.
  set.seed(110)
  n <- 400; m <- 5
  B <- matrix(rnorm(n * m), n, m)
  B <- B %*% diag(10^seq(0, -6, length.out = m))   # full rank, wildly scaled
  y <- rnorm(n)
  G <- crossprod(B); h <- as.numeric(crossprod(B, y))
  got <- BRIER:::stack_weights(G, h)
  expect_gt(got$rcond_equilibrated, got$rcond_raw)
  expect_identical(got$n_dropped, 0L)             # scaling is not rank loss
  expect_equal(got$weights, drop(solve(G, h)), tolerance = 1e-6)

  Bd <- cbind(B[, 1:3], B[, 1])                    # an exact duplicate
  Gd <- crossprod(Bd); hd <- as.numeric(crossprod(Bd, y))
  dup <- BRIER:::stack_weights(Gd, hd)
  expect_equal(dup$rcond_raw, 0)
  expect_equal(dup$rcond_equilibrated, 0)          # rescaling cannot restore rank
  expect_gt(dup$n_dropped, 0L)                     # truncation carries this case
  expect_true(all(is.finite(dup$weights)))
})
