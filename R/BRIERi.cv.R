#' Cross-validation for BRIER.I
#'
#' Performs k-fold cross-validation for \code{\link{BRIERi}}, fitting the model
#' across the eta grid and selecting the optimal eta and lambda by minimum
#' cross-validated error (cve).
#'
#' @param X A numeric matrix of predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param family A string specifying the response distribution: "gaussian",
#'   "binomial", or "poisson".
#' @param eta.list A list of numeric vectors (one per external model) specifying
#'   the eta grid for that model. A single numeric vector is replicated for
#'   all external models with a warning.
#' @param beta.external A numeric matrix of external model coefficients
#'   ((p+1) x M). Must include the intercept as the first row.
#' @param multi.method A string specifying how to combine multiple external
#'   models: "ind", "PCA", or "stacking". See \code{\link{BRIERi}}.
#' @param optim.args A list of arguments passed to \code{optim()} when using
#'   stacking with binomial or Poisson families.
#' @param ... Additional arguments passed to \code{BRIERi.eta}.
#' @param nfolds Integer. Number of CV folds.
#' @param fold Optional integer vector of fold assignments (length n). If
#'   provided, overrides \code{nfolds} and \code{seed}.
#' @param seed Optional integer seed for reproducible fold assignment.
#' @param returnY Logical. If TRUE, the per-fold predictions \code{Y} and
#'   per-fold deviances \code{E} are stored in each \code{BRIER.eta} object.
#' @param trace Logical. If TRUE, prints progress messages during fitting.
#' @param ncores Integer. Number of cores for parallel fitting.
#' @param parallel Logical. If TRUE and on a non-Windows platform, fits the eta
#'   grid in parallel using \code{parallel::mclapply}.
#'
#' @return An object of class \code{c("BRIER.cv", "BRIER")} containing:
#' \describe{
#'   \item{y}{The response vector.}
#'   \item{y.external}{A matrix of external predictions (n x M).}
#'   \item{eta.list}{The list of per-model eta grids.}
#'   \item{eta.grid}{The full combinatorial eta grid.}
#'   \item{res}{A list of \code{BRIER.eta} objects with CV fields attached
#'     (\code{cve}, \code{cvse}, \code{lambda.min}, \code{lambda.min.index}).}
#'   \item{null.dev}{The null deviance.}
#'   \item{criteria}{Always \code{"cve"} for CV.}
#'   \item{nfolds}{The number of folds used.}
#'   \item{fold}{The fold assignment vector.}
#'   \item{eta.min}{The eta combination minimising cve.}
#'   \item{eta.min.index}{The corresponding row index in \code{eta.grid}.}
#'   \item{lambda.min}{The lambda value minimising cve at \code{eta.min}.}
#'   \item{lambda.min.index}{The corresponding lambda index.}
#'   \item{eta.lambda}{A data.frame summarising the optimal lambda per eta.}
#' }
#'
#' @seealso \code{\link{BRIERi}}, \code{\link{BRIERi.selection}},
#'   \code{\link{coef.BRIER}}, \code{\link{predict.BRIER}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200
#' p <- 50
#' X <- matrix(rnorm(n * p), ncol = p)
#' beta_true <- c(rep(1, 5), rep(0, p - 5))
#' y <- X %*% beta_true + rnorm(n)
#' beta.external <- matrix(c(0, beta_true * 0.8), ncol = 1)
#'
#' cv.fit <- BRIERi.cv(
#'   X, y, family = "gaussian",
#'   eta.list = list(c(0, 0.5, 1, 2)),
#'   beta.external = beta.external,
#'   nfolds = 5, seed = 42
#' )
#' coef(cv.fit)
#' }
#'
#' @export
BRIERi.cv <- function(
  X, y, family = c("gaussian", "binomial", "poisson"),
  eta.list = c(0, exp(seq(log(0.1), log(10), length.out = 20))), 
  beta.external = rep(0, ncol(X) + 1),
  multi.method = c("ind", "PCA", "stacking"), optim.args = list(),
  ...,
  nfolds = 5, fold = NULL, seed = NULL,
  returnY = FALSE, trace = FALSE,
  ncores = max(1L, parallel::detectCores() - 1L),
  parallel = (ncores > 1L)
) {

  family <- match.arg(family)
  multi.method <- match.arg(multi.method)

  # -- Validation of y --
  if (is.data.frame(y)) { y <- as.matrix(y) }
  if (is.vector(y)) { y <- matrix(y) }
  if (ncol(y) > 1) { stop("y must be a vector or a single column matrix.", call. = FALSE) }
  if (any(is.na(y))) { stop("Missing data (NA's) detected in y.", call. = FALSE) }

  if (family == "binomial" && typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2) { stop("Attempting to use family = 'binomial' with non-binary data.", call. = FALSE) }
    if (!identical(names(tab), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(y = ", names(tab)[2], ")"))
      y <- as.double(as.character(y) == names(tab)[2])
    }
  }

  if (typeof(y) != "double") {
    tryCatch(
      storage.mode(y) <- "double",
      warning = function(w) {
        stop("y must be numeric or able to be coerced to numeric.", call. = FALSE)
      }
    )
  }

  # -- Validation of X --
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~ 0 + ., data = X), silent = TRUE)
    if (inherits(tmp, "try-error")) { stop("X must be a matrix or able to be coerced to a matrix.", call. = FALSE) }
  }
  if (storage.mode(X) == "integer") { storage.mode(X) <- "double" }
  if (any(is.na(X))) { stop("Missing data (NA's) detected in X.", call. = FALSE) }
  if (nrow(X) != nrow(y)) { stop("X and y do not have the same number of observations.", call. = FALSE) }

  # -- Compute external information --
  beta.external <- as.matrix(beta.external)
  if (any(is.na(beta.external))) { stop("Missing data (NA's) detected in beta.external.", call. = FALSE) }
  if (nrow(beta.external) != ncol(X) + 1) {
    stop("beta.external must have p+1 rows (intercept + p predictors).", call. = FALSE)
  }
  ext <- calcExtY(X, y, beta.external, family, multi.method, optim.args)
  y.external <- ext$y.external
  M <- ncol(y.external)

  # -- Validate eta.list --
  if (!is.list(eta.list)) {
    if (!is.numeric(eta.list)) { stop("eta.list must be a numeric vector or a list of numeric vectors.", call. = FALSE) }
    warning("eta.list is not a list. Replicating the same eta grid for all ", M, " external model(s).", call. = FALSE)
    eta.list <- rep(list(eta.list), M)
  }
  if (length(eta.list) != M) {
    stop(paste0("eta.list has ", length(eta.list), " element(s) but there are ", M, " external model(s)."), call. = FALSE)
  }
  if (M >= 5) {
    warning(
      "5 or more external models provided. Please consider using ensemble ",
      "weighting to aggregate multiple external models.", call. = FALSE
    )
  }
  for (k in seq_len(M)) {
    if (!is.numeric(eta.list[[k]])) { stop(paste0("eta.list[[", k, "]] must be numeric."), call. = FALSE) }
  }

  # -- Build eta grid --
  eta.grid <- as.matrix(expand.grid(eta.list))
  colnames(eta.grid) <- paste0("eta_", seq_len(M))
  if (nrow(eta.grid) > 1000) {
    warning("Large eta grid (", nrow(eta.grid), " combinations). Consider coarser eta.list grids.", call. = FALSE)
  }
  eta.grid.list <- lapply(seq_len(nrow(eta.grid)), function(i) eta.grid[i, ])
  if (trace) { cat("Total eta combinations:", nrow(eta.grid), "\n") }

  # -- Set up folds --
  n <- nrow(y)
  if (!is.null(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  if (is.null(fold)) {
    fold <- sample(seq_len(n) %% nfolds)
    fold[fold == 0] <- nfolds
  } else {
    nfolds <- max(fold)
  }

  # -- Prepare fit args --
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$family <- family
  fit.args$y.external <- y.external

  # -- Null dev --

  null.dev <- calcNullDev(
    X, y, rep(1, nrow(X)) / nrow(X),
    if (!is.null(fit.args$penalty.factor)) fit.args$penalty.factor else rep(1, ncol(X)),
    family
  )

  # -- CV for one eta combination --
  cv_one_eta <- function(eta, fit.args, fold, nfolds, n, trace, returnY) {

    if (trace) { cat("Fitting at eta = (", paste(round(eta, 3), collapse = ", "), ")\n") }

    # full data fit — returns a BRIER.eta object
    fit.args$eta <- as.numeric(eta)
    fit_full <- do.call(BRIERi.eta, fit.args)
    lambda <- fit_full$lambda

    # CV folds
    if (trace) { cat("  Starting ", nfolds, "-fold CV\n") }
    E <- Y <- matrix(NA, nrow = n, ncol = length(lambda))

    cv.args <- fit.args
    cv.args$lambda <- lambda

    for (i in seq_len(nfolds)) {
      if (trace) { cat("  Fold ", i, "\n") }
      res_fold <- cv_fold(i, cv.args, fold)
      nl <- res_fold$nl
      Y[fold == i, seq_len(nl)] <- res_fold$yhat
      E[fold == i, seq_len(nl)] <- res_fold$loss
    }

    # drop lambda values with non-finite results
    ind <- which(apply(is.finite(E), 2, all))
    E <- E[, ind, drop = FALSE]
    Y <- Y[, ind, drop = FALSE]
    lambda.cv <- lambda[ind]

    cve <- colMeans(E)
    cvse <- apply(E, 2, sd) / sqrt(n)
    min.idx <- which.min(cve)

    # attach CV results to the BRIER.eta object
    fit_full$cve            <- cve
    fit_full$cvse           <- cvse
    fit_full$lambda.cv      <- lambda.cv
    fit_full$cve.min        <- cve[min.idx]
    fit_full$lambda.min     <- lambda.cv[min.idx]
    fit_full$lambda.min.index <- min.idx
    if (returnY) {
      fit_full$Y <- Y
      fit_full$E <- E
    }

    fit_full
  }

  # -- Parallel / serial loop --
  use_parallel <- isTRUE(parallel) && ncores > 1L && .Platform$OS.type != "windows"
  if (.Platform$OS.type == "windows" && ncores > 1L) {
    message("Parallel execution not supported on Windows; using serial.")
  }

  if (!use_parallel) {
    res <- lapply(eta.grid.list, function(eta) {
      cv_one_eta(eta, fit.args, fold, nfolds, n, trace, returnY)
    })
  } else {
    res <- parallel::mclapply(
      eta.grid.list,
      function(eta) cv_one_eta(eta, fit.args, fold, nfolds, n, trace, returnY),
      mc.cores = ncores
    )
  }

  # -- Build eta.lambda summary table --
  eta.lambda <- do.call(rbind, lapply(seq_len(length(res)), function(i) {
    r <- res[[i]]
    eta_vals <- as.data.frame(t(r$eta))
    colnames(eta_vals) <- paste0("eta_", seq_along(r$eta))
    cbind(
      data.frame(eta.index = i),
      eta_vals,
      data.frame(
        criteria         = "cve",
        measure.min      = r$cve.min,
        lambda.min.index = r$lambda.min.index,
        lambda.min       = r$lambda.min
      )
    )
  }))

  # -- Find optimal --
  eta.min.index    <- which.min(eta.lambda$measure.min)
  eta.min          <- eta.grid[eta.min.index, ]
  lambda.min.index <- eta.lambda$lambda.min.index[eta.min.index]
  lambda.min       <- eta.lambda$lambda.min[eta.min.index]

  if (M == 1) {
    cat("Best eta:", round(eta.min, 3), "\n")
  } else {
    cat("Best eta: (", paste(round(eta.min, 3), collapse = ", "), ")\n")
  }
  cat("Best lambda:", round(lambda.min, 3), "\n")

  # -- Output: same structure as BRIERi, with CV fields added --
  out <- list(
    y              = y,
    y.external     = y.external,
    beta.external  = beta.external,
    family         = family,
    eta.list       = eta.list,
    eta.grid       = eta.grid,
    res            = res,           # list of BRIER.eta objects (with CV fields attached)
    null.dev       = null.dev,
    n              = nrow(X),
    p              = ncol(X),
    M              = M,

    # CV-specific fields
    criteria         = "cve",
    nfolds           = nfolds,
    fold             = fold,
    eta.min          = eta.min,
    eta.min.index    = eta.min.index,
    lambda.min       = lambda.min,
    lambda.min.index = lambda.min.index,
    eta.lambda       = eta.lambda
  )
  class(out) <- c("BRIER.cv", "BRIER")
  out
}


# -- Internal helper: single CV fold --

cv_fold <- function(i, cv.args, fold) {

  train_idx <- fold != i
  test_idx  <- fold == i

  # subset training data
  fold.args <- cv.args
  fold.args$X <- cv.args$X[train_idx, , drop = FALSE]
  fold.args$y <- cv.args$y[train_idx, , drop = FALSE]
  if (!is.null(cv.args$y.external)) {
    fold.args$y.external <- cv.args$y.external[train_idx, , drop = FALSE]
  }
  if (!is.null(cv.args$weights)) {
    fold.args$weights <- cv.args$weights[train_idx]
  }

  # test data
  X.test <- cv.args$X[test_idx, , drop = FALSE]
  y.test <- cv.args$y[test_idx, , drop = FALSE]

  # fit on training fold
  fit.i <- do.call(BRIERi.eta, fold.args)

  # predict on test fold
  yhat <- predict.BRIER.eta(fit.i, X = X.test, type = "response")
  if (is.vector(yhat)) { yhat <- matrix(yhat, ncol = 1) }

  # compute loss
  wt <- rep(1, nrow(y.test)) / nrow(y.test)
  loss <- calcDev(y.test, yhat, wt, fit.i$family)

  list(
    loss = loss,
    yhat = yhat,
    nl   = length(fit.i$lambda)
  )
}
