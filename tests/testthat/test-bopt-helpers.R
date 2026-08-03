# Argument handling for the Bayesian optimization entry points. These need no
# model fit and no optimizer, so they run everywhere.

test_that("bounds default to one [0, 10] box per external source", {
  b <- BRIER:::.bopt_bounds(NULL, 3, "external model(s)")
  expect_length(b, 3)
  expect_named(b, c("eta_1", "eta_2", "eta_3"))
  expect_true(all(vapply(b, function(x) identical(x, c(0, 10)), logical(1))))
})

test_that("user supplied bound names are preserved", {
  b <- BRIER:::.bopt_bounds(list(a = c(0, 5), b = c(0.1, 2)), 2, "x")
  expect_named(b, c("a", "b"))
  expect_equal(b$b, c(0.1, 2))
})

test_that("unnamed bounds are given canonical names", {
  b <- BRIER:::.bopt_bounds(list(c(0, 5), c(0, 2)), 2, "x")
  expect_named(b, c("eta_1", "eta_2"))
})

test_that("invalid bounds are rejected", {
  expect_error(BRIER:::.bopt_bounds(list(c(-1, 5)), 1, "x"), "cannot be negative")
  expect_error(BRIER:::.bopt_bounds(list(c(5, 5)), 1, "x"), "increasing")
  expect_error(BRIER:::.bopt_bounds(list(c(5, 1)), 1, "x"), "increasing")
  expect_error(BRIER:::.bopt_bounds(list(c(0, Inf)), 1, "x"), "length 2")
  expect_error(BRIER:::.bopt_bounds(list(c(0, 1, 2)), 1, "x"), "length 2")
  expect_error(BRIER:::.bopt_bounds(list(c(0, 5)), 2, "x"), "element")
  expect_error(BRIER:::.bopt_bounds("nonsense", 1, "x"), "named list")
})

test_that("a plain numeric init.grid is accepted only for a single source", {
  b1 <- BRIER:::.bopt_bounds(NULL, 1, "x")
  g <- BRIER:::.bopt_init_grid(c(0, 1, 2), b1, 1)
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 3)
  expect_named(g, "eta_1")

  b2 <- BRIER:::.bopt_bounds(NULL, 2, "x")
  expect_error(BRIER:::.bopt_init_grid(c(0, 1), b2, 2), "single external source")
})

test_that("init.grid must fit inside the bounds and match the source count", {
  b1 <- BRIER:::.bopt_bounds(NULL, 1, "x")
  expect_error(BRIER:::.bopt_init_grid(c(0, 99), b1, 1), "outside bounds")
  expect_error(BRIER:::.bopt_init_grid(c(0, NA), b1, 1), "finite")
  b2 <- BRIER:::.bopt_bounds(NULL, 2, "x")
  expect_error(
    BRIER:::.bopt_init_grid(data.frame(a = 1, b = 2, c = 3), b2, 2),
    "must have 2 column"
  )
  expect_null(BRIER:::.bopt_init_grid(NULL, b1, 1))
})

test_that("init.grid columns are renamed to the bound names", {
  b2 <- BRIER:::.bopt_bounds(list(lo = c(0, 5), hi = c(0, 5)), 2, "x")
  g <- BRIER:::.bopt_init_grid(data.frame(z = c(1, 2), q = c(3, 4)), b2, 2)
  expect_named(g, c("lo", "hi"))
})

test_that("the validation criteria list excludes information criteria", {
  v <- BRIER:::.bopt_val_criteria()
  expect_true(all(c("gaussian.mspe", "gaussian.rsq", "binomial.auc", "poisson.dev") %in% v))
  expect_false(any(c("BIC", "AIC", "gcv", "Cp", "GIC", "pseu.val") %in% v))
})

test_that("the optimizer accessor reports a clear error when the package is absent", {
  skip_if(requireNamespace("rBayesianOptimization", quietly = TRUE))
  expect_error(BRIER:::.bopt_optimizer(), "rBayesianOptimization")
})
