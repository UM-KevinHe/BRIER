# The three Bayesian optimization entry points, driven by a deterministic
# stand-in for the optimizer (see helper-bopt.R). These fit real models, so
# the problems are kept small.

# ---------------------------------------------------------------- BRIERi

test_that("BRIERi.bopt returns a usable selection object", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 0.5, 3, 8)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "gaussian.mspe", X.val = d$X.val, y.val = d$y.val,
    bounds = list(eta_1 = c(0, 10)), init.points = 2, n.iter = 2,
    nlambda = 20, verbose = FALSE
  )

  expect_s3_class(fit, "BRIER.bopt")
  expect_s3_class(fit, "BRIER.selection")
  expect_s3_class(fit, "BRIER")
  expect_equal(fit$criteria, "gaussian.mspe")
  expect_equal(fit$M, 1)
  expect_equal(length(fit$res), nrow(fit$eta.grid))
})

test_that("res stays aligned with eta.grid and the fits carry their own eta", {
  d <- make_individual_data()
  # deliberately unsorted, with an exact duplicate
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(3, 0, 8, 3, 0.5)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "BIC", bounds = list(eta_1 = c(0, 10)),
    init.points = 2, n.iter = 3, nlambda = 20, verbose = FALSE
  )

  expect_false(is.unsorted(fit$eta.grid[, 1]))
  expect_equal(nrow(fit$eta.grid), 4)                 # the duplicate is dropped
  expect_equal(as.numeric(fit$eta.grid[, 1]), c(0, 0.5, 3, 8))
  for (i in seq_along(fit$res)) {
    expect_lt(max(abs(fit$res[[i]]$eta - fit$eta.grid[i, ])), 1e-10)
  }
  expect_equal(fit$eta.lambda$eta.index, seq_len(nrow(fit$eta.grid)))
})

test_that("the reported optimum is the argmin of the criterion", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 1, 4, 9)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "gaussian.mspe", X.val = d$X.val, y.val = d$y.val,
    init.points = 2, n.iter = 2, nlambda = 20, verbose = FALSE
  )

  best <- which.min(fit$eta.lambda$measure.min)
  expect_equal(fit$eta.min.index, best)
  expect_equal(unname(fit$eta.min), unname(fit$eta.grid[best, ]))
  expect_equal(fit$lambda.min, fit$eta.lambda$lambda.min[best])
  expect_equal(fit$lambda.min.index, fit$eta.lambda$lambda.min.index[best])
})

test_that("coef and predict work on the result without extra arguments", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 2, 6)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "gaussian.mspe", X.val = d$X.val, y.val = d$y.val,
    init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
  )

  cf <- coef(fit)
  expect_length(cf, d$p + 1)                          # intercept plus predictors
  expect_true(all(is.finite(cf)))

  pr <- as.vector(predict(fit, X = d$X.test, type = "response"))
  expect_length(pr, nrow(d$X.test))
  expect_true(all(is.finite(pr)))

  # an eta taken from the grid resolves; one that is not on it does not
  expect_silent(predict(fit, X = d$X.test,
                        eta = matrix(fit$eta.grid[2, ], nrow = 1), type = "response"))
  expect_error(predict(fit, X = d$X.test, eta = matrix(999, nrow = 1)), "does not match")
})

test_that("information criteria need no held-out data", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 1, 5)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "BIC", init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
  )
  expect_equal(fit$criteria, "BIC")
  expect_true(is.finite(fit$lambda.min))
})

test_that("multiple external sources give one weight per source", {
  d <- make_individual_data()
  local_mocked_bindings(
    .bopt_optimizer = function() fake_optimizer(list(c(1, 1, 1), c(0, 2, 4), c(3, 0, 1)))
  )

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external3,
    multi.method = "ind", criteria = "BIC",
    bounds = list(eta_1 = c(0, 5), eta_2 = c(0, 5), eta_3 = c(0, 5)),
    init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
  )

  expect_equal(fit$M, 3)
  expect_equal(ncol(fit$eta.grid), 3)
  expect_length(fit$eta.min, 3)
  expect_true(all(c("eta_1", "eta_2", "eta_3") %in% names(fit$eta.lambda)))
  for (i in seq_along(fit$res)) {
    expect_lt(max(abs(fit$res[[i]]$eta - fit$eta.grid[i, ])), 1e-10)
  }
})

test_that("the search record is kept on the result", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 2, 7)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "BIC", bounds = list(eta_1 = c(0, 9)),
    init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
  )

  expect_equal(fit$bopt$bounds$eta_1, c(0, 9))
  expect_equal(nrow(fit$bopt$history), 3)
  expect_true(all(c("Round", "eta_1", "Value") %in% names(fit$bopt$history)))
  expect_equal(fit$bopt$n.failed, 0)
  expect_null(fit$bopt$failed.eta)
})

test_that("the score reported to the optimizer is the negated criterion", {
  # The optimizer maximizes while every BRIER criterion is minimise-better, so
  # the sign of what the objective returns is load bearing. The assembled result
  # is computed from eta.lambda and would look correct even with the sign wrong,
  # which is why this reads the score the optimizer actually saw.
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 1, 4)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "gaussian.mspe", X.val = d$X.val, y.val = d$y.val,
    init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
  )

  seen <- fit$bopt$history
  seen <- seen[order(seen$eta_1), ]
  known <- fit$eta.lambda[order(fit$eta.lambda$eta_1), ]
  expect_equal(seen$Value, -known$measure.min, tolerance = 0)
  expect_true(all(seen$Value < 0))          # MSPE is positive, so scores are negative
  expect_equal(unname(fit$bopt$best.par), unname(fit$eta.min))
})

test_that("a criterion whose natural orientation is reversed is still minimised", {
  # gaussian.rsq is negated inside the selection helper, so larger R^2 must come
  # back as a smaller measure.min and the selected eta must be the best R^2.
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 1, 4)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "gaussian.rsq", X.val = d$X.val, y.val = d$y.val,
    init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
  )

  expect_true(all(fit$eta.lambda$measure.min <= 0))     # negated R^2
  expect_equal(fit$eta.min.index, which.min(fit$eta.lambda$measure.min))
  expect_equal(fit$bopt$history$Value, -fit$eta.lambda$measure.min[
    order(order(fit$bopt$history$eta_1))], tolerance = 0)
})

test_that("an initial design is evaluated before the sequential proposals", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(6)))

  fit <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "BIC", init.grid = c(0, 1), init.points = 0, n.iter = 1,
    nlambda = 20, verbose = FALSE
  )
  expect_equal(as.numeric(fit$eta.grid[, 1]), c(0, 1, 6))
})


# ------------------------------------------------------- failure tolerance

test_that("a failing evaluation is recorded rather than aborting the run", {
  d <- make_individual_data()
  calls <- 0
  fake_with_failure <- function(FUN, bounds, init_grid_dt = NULL, init_points = 0,
                                n_iter = 1, acq = "ucb", kappa = 2.576, eps = 0,
                                kernel = NULL, verbose = TRUE) {
    scores <- c(
      FUN(eta_1 = 0)$Score,
      FUN(eta_1 = -99)$Score,        # invalid: BRIERi.eta rejects a negative eta
      FUN(eta_1 = 2)$Score
    )
    list(Best_Par = c(eta_1 = 2), Best_Value = max(scores),
         History = data.frame(Round = 1:3, eta_1 = c(0, -99, 2), Value = scores),
         Pred = list(1, NA, 2))
  }
  local_mocked_bindings(.bopt_optimizer = function() fake_with_failure)

  expect_warning(
    fit <- BRIERi.bopt(
      d$X, d$y, family = "gaussian", beta.external = d$beta.external,
      criteria = "BIC", init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
    ),
    "Fit failed at eta"
  )
  expect_equal(fit$bopt$n.failed, 1)
  expect_equal(as.numeric(fit$bopt$failed.eta), -99)
  expect_equal(nrow(fit$eta.grid), 2)                 # only the successes are kept
  expect_equal(length(fit$res), 2)
})

test_that("a run in which every evaluation fails is an error, not a bad model", {
  d <- make_individual_data()
  all_fail <- function(FUN, bounds, ...) {
    s <- c(FUN(eta_1 = -1)$Score, FUN(eta_1 = -2)$Score)
    list(Best_Par = c(eta_1 = -1), Best_Value = max(s),
         History = data.frame(Round = 1:2, eta_1 = c(-1, -2), Value = s),
         Pred = list(NA, NA))
  }
  local_mocked_bindings(.bopt_optimizer = function() all_fail)

  expect_error(
    suppressWarnings(BRIERi.bopt(
      d$X, d$y, family = "gaussian", beta.external = d$beta.external,
      criteria = "BIC", init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
    )),
    "Every eta evaluation failed"
  )
})


# --------------------------------------------- equivalence with a grid search

test_that("the objective matches the grid objective exactly", {
  d <- make_individual_data()
  etas <- c(0, 0.5, 2, 7)

  grid <- BRIERi(
    d$X, d$y, family = "gaussian", eta.list = list(etas),
    beta.external = d$beta.external, nlambda = 20, parallel = FALSE
  )
  grid_sel <- BRIERi.selection(
    grid, criteria = "gaussian.mspe", X.val = d$X.val, y.val = d$y.val
  )

  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(as.list(etas)))
  bopt <- BRIERi.bopt(
    d$X, d$y, family = "gaussian", beta.external = d$beta.external,
    criteria = "gaussian.mspe", X.val = d$X.val, y.val = d$y.val,
    init.points = 2, n.iter = 2, nlambda = 20, verbose = FALSE
  )

  g <- grid_sel$eta.lambda[order(grid_sel$eta.lambda$eta_1), ]
  b <- bopt$eta.lambda[order(bopt$eta.lambda$eta_1), ]

  expect_equal(b$measure.min, g$measure.min, tolerance = 0)
  expect_equal(b$lambda.min, g$lambda.min, tolerance = 0)
  expect_equal(b$lambda.min.index, g$lambda.min.index)
  expect_equal(unname(bopt$eta.min), unname(grid_sel$eta.min))
})


# ---------------------------------------------------------------- BRIERs

test_that("BRIERs.bopt preserves the summary-module conventions", {
  d <- make_summary_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 1, 5)))

  fit <- BRIERs.bopt(
    d$sumstats, d$XtX, family = "gaussian", beta.external = d$beta.external,
    criteria = "gaussian.mspe", X.val = d$X.val, y.val = d$y.val,
    init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
  )

  expect_s3_class(fit, "BRIER.bopt")
  expect_equal(fit$p, nrow(d$sumstats))
  expect_false(is.null(fit$XtY))
  expect_false(is.null(fit$varnames))
  expect_length(coef(fit), nrow(d$sumstats))          # no intercept row
  for (i in seq_along(fit$res)) {
    expect_lt(max(abs(fit$res[[i]]$eta - fit$eta.grid[i, ])), 1e-10)
  }
})

test_that("BRIERs.bopt supports summary information criteria", {
  d <- make_summary_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 2, 6)))

  fit <- BRIERs.bopt(
    d$sumstats, d$XtX, family = "gaussian", beta.external = d$beta.external,
    criteria = "GIC", TN = d$n, h2 = 0,
    init.points = 2, n.iter = 1, nlambda = 20, verbose = FALSE
  )
  expect_equal(fit$criteria, "GIC")
  expect_true(is.finite(fit$lambda.min))
})

test_that("BRIERs.bopt requires TN for the information criteria", {
  d <- make_summary_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0, 1)))
  expect_error(
    BRIERs.bopt(
      d$sumstats, d$XtX, family = "gaussian", beta.external = d$beta.external,
      criteria = "GIC", init.points = 2, n.iter = 1, verbose = FALSE
    ),
    "TN"
  )
})


# -------------------------------------------------------------- BRIERfull

test_that("BRIERfull.bopt weights the external cohorts and keeps the cohort vector", {
  d <- make_pooled_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(0.5, 1, 3)))

  fit <- BRIERfull.bopt(
    d$X, d$y, d$cohort, family = "gaussian", criteria = "gaussian.mspe",
    X.val = d$X.val, y.val = d$y.val, init.points = 2, n.iter = 1,
    nlambda = 20, verbose = FALSE
  )

  expect_s3_class(fit, "BRIER.bopt")
  expect_equal(fit$cohort, d$cohort)
  expect_equal(fit$M, 1)
  expect_equal(length(fit$res), nrow(fit$eta.grid))
  expect_length(as.vector(predict(fit, X = d$X.val, type = "response")), nrow(d$X.val))
})

test_that("BRIERfull.bopt refuses a pool with no external cohort", {
  d <- make_pooled_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(1)))
  expect_error(
    BRIERfull.bopt(
      d$X, d$y, rep(0L, length(d$cohort)), family = "gaussian",
      criteria = "gaussian.mspe", X.val = d$X.val, y.val = d$y.val,
      init.points = 0, n.iter = 1, verbose = FALSE
    ),
    "At least one external cohort"
  )
})


# ------------------------------------------------------------- guardrails

test_that("validation criteria require held-out data", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(1)))
  expect_error(
    BRIERi.bopt(
      d$X, d$y, family = "gaussian", beta.external = d$beta.external,
      criteria = "gaussian.mspe", init.points = 0, n.iter = 1, verbose = FALSE
    ),
    "X.val and y.val are required"
  )
})

test_that("BRIERfull.bopt requires validation data, having no information criteria", {
  d <- make_pooled_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(1)))
  expect_error(
    BRIERfull.bopt(
      d$X, d$y, d$cohort, family = "gaussian", criteria = "gaussian.mspe",
      init.points = 0, n.iter = 1, verbose = FALSE
    ),
    "validation criteria only"
  )
})

test_that("an initialization too small to seed the surrogate is refused", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(1)))
  expect_error(
    BRIERi.bopt(
      d$X, d$y, family = "gaussian", beta.external = d$beta.external,
      criteria = "BIC", init.points = 1, n.iter = 2, verbose = FALSE
    ),
    "init.points >= 2"
  )
  expect_error(
    BRIERi.bopt(
      d$X, d$y, family = "gaussian", beta.external = d$beta.external,
      criteria = "BIC", init.points = 4, n.iter = 0, verbose = FALSE
    ),
    "n.iter must be at least 1"
  )
})

test_that("a mis-shaped external coefficient matrix is refused", {
  d <- make_individual_data()
  local_mocked_bindings(.bopt_optimizer = function() fake_optimizer(list(1)))
  expect_error(
    BRIERi.bopt(
      d$X, d$y, family = "gaussian",
      beta.external = matrix(0, nrow = d$p, ncol = 1),   # missing the intercept row
      criteria = "BIC", init.points = 0, n.iter = 1, verbose = FALSE
    ),
    "Please include intercept"
  )
})
