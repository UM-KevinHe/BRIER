# .stop_on_worker_errors (R/utils-parallel.R): a failed parallel::mclapply worker must stop the fit
# instead of being returned as a "result" (2026-11-13: an invalid MCP gamma failed in every BRIERfull
# worker, the fit returned, and each selection then died on an unrelated error).

test_that("good worker results pass through unchanged", {
  res = list(list(a = 1), list(a = 2))
  expect_identical(BRIER:::.stop_on_worker_errors(res, list(0, 1)), res)
})

test_that("a try-error from a worker stops with its message and the eta", {
  err = structure("Error in f() : gamma must be greater than 1 for MCP\n", class = "try-error")
  expect_error(
    BRIER:::.stop_on_worker_errors(list(list(a = 1), err), list(0, 0.5)),
    "eta = 0.5.*gamma must be greater than 1 for MCP"
  )
})

test_that("a killed worker (NULL result) stops", {
  expect_error(
    BRIER:::.stop_on_worker_errors(list(NULL, list(a = 1)), list(2, 3)),
    "returned no result"
  )
})

test_that("a parallel fit with an invalid penalty now errors instead of returning", {
  skip_on_os("windows")
  set.seed(1)
  n = 60; p = 8
  X = matrix(rnorm(n * p), n, p)
  y = rnorm(n)
  b_ext = matrix(c(0, rnorm(p)), ncol = 1)
  expect_error(
    BRIERi(X, y, beta.external = b_ext, eta.list = c(0, 1), penalty = "MCP", gamma = 0.5,
           ncores = 2, parallel = TRUE),
    "fitting failed in the parallel workers|gamma"
  )
})
