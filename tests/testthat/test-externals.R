# Tests for the external-panel preparation: de-duplication, which every
# multi-source fit runs first, and the principal-component reduction behind
# multi.method = "PCstacking".
#
# Both steps read only beta.external, so they are tested directly on small
# hand-built panels where the right answer is known by construction, and then
# once through each fitter to confirm the wiring: that M, the eta grid and the
# per-source settings all follow the panel that survived.


# -- Fixtures ----------------------------------------------------------------
#
# b.sig and b.alt are informative and orthogonal to each other; b.dup is a pure
# rescaling of b.sig, so it is the same model and carries no new information.

ext_panel <- function(p = 25) {
  b.sig <- c(rep(1.2, 5), rep(0, p - 5))
  b.alt <- c(rep(0, p - 5), rep(1, 5))
  list(
    p     = p,
    sig   = b.sig,
    alt   = b.alt,
    dup   = b.sig * 0.5,
    near  = b.sig + 0.05 * b.alt + c(0.01, rep(0, p - 1)),
    # Nearly the negation of b.sig, at a cosine of about -0.99. The extra block
    # keeps it outside the span of b.sig and b.alt, so only the similarity rule
    # can catch it.
    neg   = -b.sig + c(rep(0, 10), rep(0.15, 5), rep(0, p - 15)),
    zero  = rep(0, p)
  )
}

# Three models at cosine 0.94 along a chain: A to B and B to C are both above a
# 0.9 threshold, but A to C is 0.94^2 = 0.884, below it. Each step introduces a
# direction the previous pair does not span.
similarity_chain <- function(p = 25, rho = 0.94) {
  blk <- function(i) { v <- numeric(p); v[i] <- 1 / sqrt(5); v }
  u <- blk(1:5); v <- blk(6:10); w <- blk(11:15)
  s <- sqrt(1 - rho^2)
  A <- u
  B <- rho * u + s * v
  C <- rho * B + s * w
  cbind(A = A, B = B, C = C)
}

# Four distinct models sharing one dominant direction: full rank, so nothing is
# de-duplicated, but the variance concentrates in the leading directions and a
# principal-component reduction has something to truncate.
correlated_panel <- function(p = 25, seed = 7) {
  set.seed(seed)
  common <- c(rep(1.2, 5), rep(0, p - 5))
  B <- sapply(1:4, function(k) common * runif(1, 0.6, 1.4) + rnorm(p, sd = 0.6))
  colnames(B) <- paste0("m", 1:4)
  B
}


# -- dedupExternals ----------------------------------------------------------

test_that("a rescaled copy is dropped and the first occurrence is kept", {
  e <- ext_panel()
  d <- dedupExternals(cbind(a = e$sig, b = e$dup, c = e$alt), cor.min = 0.9)

  expect_true(d$applied)
  expect_equal(d$keep, c(1L, 3L))
  expect_equal(ncol(d$beta.external), 2)
  expect_equal(as.numeric(d$beta.external[, 1]), e$sig)
  expect_equal(d$dropped$index, 2L)
  expect_equal(d$dropped$duplicate.of, "a")
})

test_that("which copy survives is decided by position, not by magnitude", {
  e <- ext_panel()
  # The larger model second: the smaller one still survives, because it came
  # first. The rule must not depend on the coefficients' scale.
  d <- dedupExternals(cbind(e$dup, e$sig, e$alt), cor.min = 0.9)
  expect_equal(d$keep, c(1L, 3L))
  expect_equal(as.numeric(d$beta.external[, 1]), e$dup)
})

test_that("a set of three copies collapses to its first member", {
  e <- ext_panel()
  d <- dedupExternals(cbind(e$sig, e$sig * 2, e$sig * 3, e$alt), cor.min = 0.9)
  expect_equal(d$keep, c(1L, 4L))
  expect_equal(nrow(d$dropped), 2)
})

test_that("the threshold decides what counts as near", {
  e <- ext_panel()
  panel <- cbind(e$sig, e$near, e$alt)
  # e$near sits at a cosine of roughly 0.999 with e$sig.
  expect_equal(dedupExternals(panel, cor.min = 0.9)$keep, c(1L, 3L))
  expect_equal(dedupExternals(panel, cor.min = 0.99999)$keep, 1:3)
})

test_that("a near-copy with its sign flipped is still a duplicate", {
  e <- ext_panel()
  # Similarity is measured in absolute value. A sign-flipped model spans the
  # same direction, so it is just as unidentifiable in the stacking solve, and
  # an exact flip is not required for that to bite.
  d <- dedupExternals(cbind(a = e$sig, b = e$neg, c = e$alt), cor.min = 0.9)
  expect_equal(d$keep, c(1L, 3L))
  expect_equal(d$dropped$duplicate.of, "a")
})

test_that("similarity is measured only against the models still kept", {
  # B is dropped as a copy of A. C is close to B but not to A, so with B gone
  # there is nothing left for C to duplicate and it must survive: dropping it
  # would remove a model on the strength of one that is not in the fit.
  chain <- similarity_chain()
  d <- dedupExternals(chain, cor.min = 0.9)
  expect_equal(d$keep, c(1L, 3L))
  expect_equal(d$dropped$index, 2L)
})

test_that("a column the kept set spans exactly is caught by the backstop", {
  e <- ext_panel()
  # Uncorrelated with either of its parents, so pairwise similarity cannot see
  # it, but the stacking Gram is singular all the same.
  d <- dedupExternals(cbind(e$sig, e$alt, e$sig + e$alt), cor.min = 0.9)
  expect_equal(d$keep, c(1L, 2L))
  expect_true(grepl("linearly dependent", d$dropped$reason))
})

test_that("an all-zero external is dropped", {
  e <- ext_panel()
  d <- dedupExternals(cbind(e$sig, e$zero, e$alt), cor.min = 0.9)
  expect_equal(d$keep, c(1L, 3L))
  expect_true(grepl("all zero", d$dropped$reason))
})

test_that("an entirely zero panel is reported and left alone", {
  e <- ext_panel()
  d <- dedupExternals(cbind(e$zero, e$zero), cor.min = 0.9)
  # Returning a zero-column panel would break the fitter for no gain.
  expect_false(d$applied)
  expect_equal(ncol(d$beta.external), 2)
  expect_true(length(d$note) > 0)
})

test_that("the intercept row is held out of the comparison and carried along", {
  e <- ext_panel()
  panel <- rbind(c(3, 5, 7), cbind(e$sig, e$dup, e$alt))
  d <- dedupExternals(panel, cor.min = 0.9, intercept.row = TRUE)

  expect_equal(d$keep, c(1L, 3L))
  expect_equal(nrow(d$beta.external), e$p + 1)
  expect_equal(as.numeric(d$beta.external[1, ]), c(3, 7))
})

test_that("intercepts alone never make two models look different", {
  e <- ext_panel()
  # Same predictor coefficients, different intercepts. The intercept is not a
  # direction, so these are still the same model.
  panel <- rbind(c(0, 99), cbind(e$sig, e$sig))
  expect_equal(dedupExternals(panel, cor.min = 0.9, intercept.row = TRUE)$keep, 1L)
})

test_that("a single external and a disabled threshold are both no-ops", {
  e <- ext_panel()
  one <- dedupExternals(matrix(e$sig, ncol = 1), cor.min = 0.9)
  expect_false(one$applied)
  expect_equal(ncol(one$beta.external), 1)

  off <- dedupExternals(cbind(e$sig, e$dup), cor.min = NULL)
  expect_false(off$applied)
  expect_equal(off$keep, 1:2)
  expect_equal(dedupExternals(cbind(e$sig, e$dup), cor.min = NA)$keep, 1:2)
})

test_that("an out-of-range threshold is refused", {
  e <- ext_panel()
  panel <- cbind(e$sig, e$alt)
  expect_error(dedupExternals(panel, cor.min = 0), "in \\(0, 1\\]")
  expect_error(dedupExternals(panel, cor.min = 1.5), "in \\(0, 1\\]")
})


# -- reduceExternalsPCA ------------------------------------------------------

test_that("the fewest components reaching the variance target are kept", {
  B <- correlated_panel()
  r <- reduceExternalsPCA(B, pca.var = 0.8)

  expect_true(r$applied)
  expect_true(r$n.pcs < ncol(B))
  expect_true(r$cumulative[r$n.pcs] >= 0.8)
  # Minimal: one component fewer would fall short.
  expect_true(r$n.pcs == 1 || r$cumulative[r$n.pcs - 1] < 0.8)
  expect_equal(ncol(r$beta.external), r$n.pcs)
  expect_equal(colnames(r$beta.external), paste0("PC", seq_len(r$n.pcs)))
})

test_that("pca.var moves the number of components kept", {
  B <- correlated_panel()
  expect_true(reduceExternalsPCA(B, pca.var = 0.5)$n.pcs <=
                reduceExternalsPCA(B, pca.var = 0.95)$n.pcs)
  expect_equal(reduceExternalsPCA(B, pca.var = 1)$n.pcs, ncol(B))
})

test_that("the result is the projection of the models, on the coefficient scale", {
  B <- correlated_panel()
  r <- reduceExternalsPCA(B, pca.var = 0.8)
  # B %*% V, not an orthonormal basis: eta shrinks toward these coefficients, so
  # a rescaled basis would silently rescale what eta means.
  expect_equal(unname(r$beta.external), unname(B %*% r$rotation))
  expect_equal(nrow(r$beta.external), nrow(B))
})

test_that("the decomposition is uncentred", {
  B <- correlated_panel()
  r <- reduceExternalsPCA(B, pca.var = 0.8)
  expect_equal(unname(r$rotation), unname(svd(B)$v[, seq_len(r$n.pcs), drop = FALSE]))
  # A centred decomposition would subtract a mean model, which is not a
  # meaningful reference when each column is a model rather than a draw.
  centred <- prcomp(B, center = TRUE)$rotation[, seq_len(r$n.pcs), drop = FALSE]
  expect_false(isTRUE(all.equal(unname(r$rotation), unname(centred))))
})

test_that("the intercept row is reduced by the same combination as its column", {
  B <- correlated_panel()
  icpt <- c(1, 2, 3, 4)
  r <- reduceExternalsPCA(rbind(icpt, B), pca.var = 0.8, intercept.row = TRUE)

  expect_equal(nrow(r$beta.external), nrow(B) + 1)
  expect_equal(as.numeric(r$beta.external[1, ]), as.numeric(icpt %*% r$rotation))
  expect_equal(unname(r$beta.external[-1, , drop = FALSE]), unname(B %*% r$rotation))
})

test_that("a single external and an invalid target are handled", {
  B <- correlated_panel()
  one <- reduceExternalsPCA(B[, 1, drop = FALSE], pca.var = 0.8)
  expect_false(one$applied)
  expect_equal(ncol(one$beta.external), 1)

  expect_error(reduceExternalsPCA(B, pca.var = 0), "in \\(0, 1\\]")
  expect_error(reduceExternalsPCA(B, pca.var = 2), "in \\(0, 1\\]")
})


# -- The fitters -------------------------------------------------------------

test_that("BRIERi drops a duplicate and the eta grid follows the survivors", {
  d <- make_individual_data()
  e <- ext_panel(d$p)
  panel <- rbind(0, cbind(e$sig, e$dup, e$alt))

  fit <- suppressMessages(BRIERi(
    d$X, d$y, family = "gaussian",
    eta.list = list(c(0, 1), c(0, 2), c(0, 3)),
    beta.external = panel, multi.method = "ind",
    penalty = "LASSO", nlambda = 10, parallel = FALSE, ncores = 1
  ))

  expect_equal(fit$M, 2)
  expect_equal(ncol(fit$eta.grid), 2)
  # The surviving sources keep their own eta grids rather than the fit failing
  # on a length mismatch the caller had no way to anticipate.
  expect_equal(fit$eta.list, list(c(0, 1), c(0, 3)))
  expect_equal(ncol(fit$beta.external), 2)
  expect_equal(fit$external.dedup$dropped$index, 2L)
})

test_that("the de-duplication is announced rather than silent", {
  d <- make_individual_data()
  e <- ext_panel(d$p)
  panel <- rbind(0, cbind(e$sig, e$dup))

  expect_message(
    BRIERi(d$X, d$y, family = "gaussian", eta.list = list(c(0, 1), c(0, 1)),
           beta.external = panel, multi.method = "ind",
           penalty = "LASSO", nlambda = 10, parallel = FALSE, ncores = 1),
    "dedupExternals"
  )
})

test_that("dedup.cor = NULL restores the un-screened panel", {
  d <- make_individual_data()
  e <- ext_panel(d$p)
  panel <- rbind(0, cbind(e$sig, e$dup, e$alt))

  fit <- BRIERi(
    d$X, d$y, family = "gaussian",
    eta.list = list(c(0, 1), c(0, 2), c(0, 3)),
    beta.external = panel, multi.method = "ind", dedup.cor = NULL,
    penalty = "LASSO", nlambda = 10, parallel = FALSE, ncores = 1
  )

  expect_equal(fit$M, 3)
  expect_null(fit$external.dedup)
})

test_that("a panel with nothing redundant is fitted identically either way", {
  d <- make_individual_data()
  e <- ext_panel(d$p)
  panel <- rbind(0, cbind(e$sig, e$alt))
  args <- list(
    X = d$X, y = d$y, family = "gaussian", eta.list = list(c(0, 1), c(0, 1)),
    beta.external = panel, multi.method = "ind",
    penalty = "LASSO", nlambda = 10, parallel = FALSE, ncores = 1
  )

  on.  <- do.call(BRIERi, args)
  off. <- do.call(BRIERi, c(args, list(dedup.cor = NULL)))
  expect_equal(on.$res[[2]]$beta, off.$res[[2]]$beta)
})

test_that("de-duplication rescues a stacking solve that is otherwise singular", {
  d <- make_summary_data()
  e <- ext_panel(nrow(d$sumstats))
  panel <- cbind(e$sig, e$dup, e$alt)

  # t(B) XtX B is singular the moment two columns carry the same information.
  expect_error(
    BRIERs(d$sumstats, d$XtX, family = "gaussian", eta.list = c(0, 1),
           beta.external = panel, multi.method = "stacking", dedup.cor = NULL,
           nlambda = 10, parallel = FALSE, ncores = 1)
  )

  fit <- suppressWarnings(suppressMessages(BRIERs(
    d$sumstats, d$XtX, family = "gaussian", eta.list = c(0, 1),
    beta.external = panel, multi.method = "stacking",
    nlambda = 10, parallel = FALSE, ncores = 1
  )))
  expect_s3_class(fit, "BRIER")
  expect_equal(fit$M, 1)
})

test_that("BRIERi.cv drops a duplicate too", {
  d <- make_individual_data()
  e <- ext_panel(d$p)
  panel <- rbind(0, cbind(e$sig, e$dup, e$alt))

  fit <- suppressMessages(BRIERi.cv(
    d$X, d$y, family = "gaussian",
    eta.list = list(c(0, 1), c(0, 2), c(0, 3)),
    beta.external = panel, multi.method = "ind",
    penalty = "LASSO", nlambda = 10, nfolds = 3, seed = 1,
    parallel = FALSE, ncores = 1
  ))

  expect_equal(fit$M, 2)
  expect_equal(ncol(fit$eta.grid), 2)
})


# -- multi.method = "PCstacking" ------------------------------------------------

test_that("PCstacking reduces then stacks, leaving a single eta", {
  d <- make_individual_data()
  B <- correlated_panel(d$p)
  panel <- rbind(0, B)

  fit <- suppressMessages(BRIERi(
    d$X, d$y, family = "gaussian", eta.list = c(0, 1, 5),
    beta.external = panel, multi.method = "PCstacking",
    penalty = "LASSO", nlambda = 10, parallel = FALSE, ncores = 1
  ))

  expect_equal(fit$M, 1)
  expect_equal(ncol(fit$eta.grid), 1)
  expect_true(fit$external.pca$applied)
  expect_true(fit$external.pca$n.pcs < ncol(B))
  # The stored panel is the one the caller supplied, so plot.box still compares
  # against the real external models rather than against components.
  expect_equal(ncol(fit$beta.external), ncol(B))
})

test_that("PCstacking is not the same fit as plain stacking", {
  d <- make_individual_data()
  panel <- rbind(0, correlated_panel(d$p))
  args <- list(
    X = d$X, y = d$y, family = "gaussian", eta.list = c(0, 1),
    beta.external = panel, penalty = "LASSO", nlambda = 10,
    parallel = FALSE, ncores = 1
  )

  pcs <- suppressMessages(do.call(BRIERi, c(args, list(multi.method = "PCstacking"))))
  stk <- suppressMessages(do.call(BRIERi, c(args, list(multi.method = "stacking"))))
  expect_false(isTRUE(all.equal(pcs$y.external, stk$y.external)))
})

test_that("pca.var reaches the fitters rather than being fixed at its default", {
  d <- make_individual_data()
  panel <- rbind(0, correlated_panel(d$p))
  args <- list(
    X = d$X, y = d$y, family = "gaussian", eta.list = c(0, 1),
    beta.external = panel, multi.method = "PCstacking",
    penalty = "LASSO", nlambda = 10, parallel = FALSE, ncores = 1
  )

  loose  <- suppressMessages(do.call(BRIERi, c(args, list(pca.var = 0.5))))
  strict <- suppressMessages(do.call(BRIERi, c(args, list(pca.var = 0.999))))
  expect_true(strict$external.pca$n.pcs > loose$external.pca$n.pcs)
})

test_that("PCstacking works in the summary module and in cross-validation", {
  ds <- make_summary_data()
  Bs <- correlated_panel(nrow(ds$sumstats))
  s <- suppressWarnings(suppressMessages(BRIERs(
    ds$sumstats, ds$XtX, family = "gaussian", eta.list = c(0, 1),
    beta.external = Bs, multi.method = "PCstacking",
    nlambda = 10, parallel = FALSE, ncores = 1
  )))
  expect_equal(s$M, 1)
  expect_true(s$external.pca$applied)

  d <- make_individual_data()
  cv <- suppressWarnings(suppressMessages(BRIERi.cv(
    d$X, d$y, family = "gaussian", eta.list = c(0, 1),
    beta.external = rbind(0, correlated_panel(d$p)), multi.method = "PCstacking",
    penalty = "LASSO", nlambda = 10, nfolds = 3, seed = 1,
    parallel = FALSE, ncores = 1
  )))
  expect_equal(cv$M, 1)
  expect_true(cv$external.pca$applied)
})

test_that("a PCstacking fit is a plain BRIER object downstream", {
  d <- make_individual_data()
  fit <- suppressMessages(BRIERi(
    d$X, d$y, family = "gaussian", eta.list = c(0, 1, 5),
    beta.external = rbind(0, correlated_panel(d$p)), multi.method = "PCstacking",
    penalty = "LASSO", nlambda = 10, parallel = FALSE, ncores = 1
  ))

  sel <- BRIERi.selection(fit, criteria = "gaussian.mspe",
                          X.val = d$X.val, y.val = d$y.val)
  expect_length(coef(sel), d$p + 1)
  expect_length(predict(sel, X = d$X.val), nrow(d$X.val))
})


# -- The Bayesian-optimization entry points ----------------------------------

test_that("BRIERi.bopt subsets the search box to the surviving sources", {
  d <- make_individual_data()
  e <- ext_panel(d$p)
  panel <- rbind(0, cbind(e$sig, e$dup, e$alt))

  # Points of length 2, not 3: if the box were not subset with the panel the
  # objective would still expect three arguments and this would not run.
  local_mocked_bindings(
    .bopt_optimizer = function() fake_optimizer(list(c(1, 1), c(0, 2), c(3, 1)))
  )

  fit <- suppressMessages(BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = panel,
    multi.method = "ind", criteria = "BIC",
    bounds = list(eta_1 = c(0, 5), eta_2 = c(0, 5), eta_3 = c(0, 5)),
    init.points = 2, n.iter = 1, nlambda = 10, verbose = FALSE
  ))

  expect_equal(fit$M, 2)
  expect_equal(ncol(fit$eta.grid), 2)
  expect_named(fit$bopt$bounds, c("eta_1", "eta_2"))
  expect_equal(fit$external.dedup$dropped$index, 2L)
})

test_that("BRIERs.bopt accepts PCstacking and searches a single eta", {
  ds <- make_summary_data()
  Bs <- correlated_panel(nrow(ds$sumstats))

  local_mocked_bindings(
    .bopt_optimizer = function() fake_optimizer(list(0.5, 2, 4))
  )

  fit <- suppressWarnings(suppressMessages(BRIERs.bopt(
    ds$sumstats, ds$XtX, family = "gaussian", beta.external = Bs,
    multi.method = "PCstacking", criteria = "gaussian.mspe",
    X.val = ds$X.val, y.val = ds$y.val,
    init.points = 2, n.iter = 1, nlambda = 10, verbose = FALSE
  )))

  expect_equal(fit$M, 1)
  expect_true(fit$external.pca$applied)
  expect_s3_class(fit, "BRIER.bopt")
})
