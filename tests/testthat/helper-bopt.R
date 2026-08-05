# Helpers for the Bayesian optimization tests.
#
# The optimizer itself is a suggested dependency, so these tests do not use it.
# They substitute a deterministic stand-in that honours the same contract:
# call FUN with one named argument per bound, collect Score and Pred, and
# return History, Best_Par and Best_Value. What is under test is the BRIER
# code around the optimizer (the objective, the caching, the assembly of the
# result), not the search strategy, so a deterministic sequence of points is
# both sufficient and more reproducible than a real search.


# Deterministic stand-in for rBayesianOptimization::BayesianOptimization.
# `points` fixes the sequence of proposals so every test is reproducible. It
# stands in for the random initial draws and the sequential proposals together,
# so a caller must set init.points + n.iter to its length; the check below keeps
# the tests honest about the argument combination they claim to exercise. Rows
# of init_grid_dt are evaluated first, as the real optimizer does.
fake_optimizer <- function(points) {
  force(points)
  function(FUN, bounds, init_grid_dt = NULL, init_points = 0, n_iter = 1,
           acq = "ucb", kappa = 2.576, eps = 0, kernel = NULL, verbose = TRUE) {
    if (length(points) != init_points + n_iter) {
      stop(sprintf(
        "fake_optimizer: %d point(s) supplied but init.points + n.iter = %d.",
        length(points), init_points + n_iter
      ), call. = FALSE)
    }
    nm <- names(bounds)
    pts <- list()
    if (!is.null(init_grid_dt)) {
      for (i in seq_len(nrow(init_grid_dt))) {
        pts[[length(pts) + 1L]] <- as.numeric(init_grid_dt[i, ])
      }
    }
    for (i in seq_along(points)) {
      pts[[length(pts) + 1L]] <- as.numeric(points[[i]])
    }
    scores <- numeric(0)
    preds <- list()
    for (i in seq_along(pts)) {
      r <- do.call(FUN, stats::setNames(as.list(pts[[i]]), nm))
      scores <- c(scores, r$Score)
      preds[[i]] <- r$Pred
    }
    history <- cbind(
      data.frame(Round = seq_along(pts)),
      as.data.frame(do.call(rbind, pts)),
      data.frame(Value = scores)
    )
    colnames(history) <- c("Round", nm, "Value")
    best <- which.max(scores)
    list(
      Best_Par   = stats::setNames(as.numeric(pts[[best]]), nm),
      Best_Value = scores[best],
      History    = history,
      Pred       = preds
    )
  }
}


# Small individual-level problem: cheap to fit, but with real signal so the
# criterion actually varies with eta.
make_individual_data <- function(n = 80, p = 25, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(rep(1.5, 5), rep(0, p - 5))
  y <- as.numeric(X %*% beta + rnorm(n))
  Xv <- matrix(rnorm(n * p), n, p)
  yv <- as.numeric(Xv %*% beta + rnorm(n))
  Xt <- matrix(rnorm(n * p), n, p)
  yt <- as.numeric(Xt %*% beta + rnorm(n))
  list(
    X = X, y = y, X.val = Xv, y.val = yv, X.test = Xt, y.test = yt,
    beta = beta, n = n, p = p,
    # one informative external model, intercept row first
    beta.external = matrix(c(0, beta * 0.8), ncol = 1),
    # Three external models of decreasing usefulness. They must not be
    # proportional to one another: a rescaling is the same model, and the
    # fitters drop it before fitting (see dedupExternals), which would leave a
    # test that asked for three sources looking at two. The second model
    # therefore carries its own support on variables 6 to 10, so its cosine
    # similarity with the first is about 0.71 rather than 1.
    beta.external3 = cbind(
      c(0, beta * 0.8),
      c(0, beta * 0.4 + c(rep(0, 5), rep(0.6, 5), rep(0, p - 10))),
      c(0, rev(beta) * 0.2)
    )
  )
}


# Matching summary-level problem built from the same generator.
make_summary_data <- function(n = 80, p = 25, seed = 1) {
  d <- make_individual_data(n, p, seed)
  Xs <- standardize_X(d$X)[[1]]
  ys <- as.numeric(standardize_X(as.matrix(d$y))[[1]])
  ld <- calLD(d$X)
  sumstats <- data.frame(
    varnames = paste0("V", seq_len(ncol(Xs))),
    corr     = as.numeric(crossprod(Xs, ys)) / length(ys),
    n        = d$n
  )
  list(
    sumstats = sumstats,
    XtX = ld$XtX,
    nz = ld$nz,
    beta.external = matrix(d$beta * 0.8, ncol = 1),
    X.val = standardize_X(d$X.val)[[1]],
    y.val = as.numeric(standardize_X(as.matrix(d$y.val))[[1]]),
    n = d$n
  )
}


# Pooled problem for BRIERfull: target cohort 0 plus one external cohort.
make_pooled_data <- function(n = 60, p = 25, seed = 2) {
  d <- make_individual_data(n, p, seed)
  list(
    X = rbind(d$X, d$X.test),
    y = c(d$y, d$y.test),
    cohort = c(rep(0L, nrow(d$X)), rep(1L, nrow(d$X.test))),
    X.val = d$X.val,
    y.val = d$y.val
  )
}
