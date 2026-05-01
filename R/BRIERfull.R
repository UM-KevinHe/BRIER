#' BRIER.FULL: integration of individual-level external data
#'
#' Fits a penalized regression model that integrates individual-level data from
#' one or more external cohorts with a target cohort. External observations are
#' assigned per-cohort weights (eta) which are tuned over a user-specified grid.
#'
#' @param X A numeric matrix of predictors (n x p) combining target and external
#'   cohort observations.
#' @param y A numeric response vector of length n.
#' @param cohort An integer vector of length n indicating cohort membership.
#'   Use 0 for the target cohort and 1, 2, ... for external cohorts.
#' @param family A string specifying the response distribution: "gaussian",
#'   "binomial", or "poisson".
#' @param eta.list A list of numeric vectors (one per external cohort) specifying
#'   the eta grid for that cohort. A single numeric vector is replicated for
#'   all external cohorts with a warning.
#' @param ... Additional arguments passed to \code{BRIERi_fit} (e.g. \code{penalty},
#'   \code{alpha}, \code{nlambda}, \code{penalty.factor}).
#' @param trace Logical. If TRUE, prints progress messages during fitting.
#' @param ncores Integer. Number of cores for parallel fitting.
#' @param parallel Logical. If TRUE and on a non-Windows platform, fits the eta
#'   grid in parallel using \code{parallel::mclapply}.
#'
#' @return An object of class \code{"BRIER"} containing:
#' \describe{
#'   \item{y}{The response vector.}
#'   \item{y.external}{NULL (placeholder for compatibility with \code{BRIERi}).}
#'   \item{family}{The response family.}
#'   \item{cohort}{The cohort indicator vector.}
#'   \item{eta.list}{The list of per-cohort eta grids.}
#'   \item{eta.grid}{The full combinatorial eta grid (matrix, rows are combinations).}
#'   \item{res}{A list of \code{BRIER.fit} objects, one per eta combination.}
#'   \item{null.dev}{The null deviance computed on the target cohort only.}
#' }
#'
#' @seealso \code{\link{BRIERi}}, \code{\link{BRIERi_fit}},
#'   \code{\link{BRIERi.selection}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n_target <- 200
#' n_ext <- 500
#' p <- 50
#' X <- matrix(rnorm((n_target + n_ext) * p), ncol = p)
#' beta <- c(rep(1, 5), rep(0, p - 5))
#' y <- X %*% beta + rnorm(n_target + n_ext)
#' cohort <- c(rep(0, n_target), rep(1, n_ext))
#'
#' fit <- BRIERfull(
#'   X, y, cohort, family = "gaussian",
#'   eta.list = list(c(0, 0.5, 1, 2)),
#'   penalty = "LASSO"
#' )
#' }
#'
#' @export
BRIERfull <- function(
  X, y, cohort, family = c("gaussian", "binomial", "poisson"),
  eta.list,
  ...,
  trace = FALSE,
  ncores = max(1L, parallel::detectCores() - 1L),
  parallel = (ncores > 1L)
) {

  family <- match.arg(family)

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

  # -- Validation of cohort --
  cohort <- as.integer(cohort)
  if (length(cohort) != nrow(X)) { stop("Length of cohort must equal nrow(X).", call. = FALSE) }
  if (!any(cohort == 0)) { stop("cohort must contain target cohort coded as 0.", call. = FALSE) }
  if (any(cohort < 0)) { stop("cohort must be >= 0 (0 = target, 1..M = external).", call. = FALSE) }
  ext_cohorts <- sort(unique(cohort[cohort != 0]))
  M <- length(ext_cohorts)
  if (M == 0) { stop("At least one external cohort (cohort > 0) is required.", call. = FALSE) }

  # -- Validate eta.list --
  if (!is.list(eta.list)) {
    if (!is.numeric(eta.list)) { stop("eta.list must be a numeric vector or a list of numeric vectors.", call. = FALSE) }
    warning("eta.list is not a list. Replicating the same eta grid for all ",
            M, " external cohort(s).", call. = FALSE)
    eta.list <- rep(list(eta.list), M)
  }
  if (length(eta.list) != M) {
    stop(paste0("eta.list has ", length(eta.list), " element(s) but there are ",
                M, " external cohort(s)."), call. = FALSE)
  }
  for (k in seq_len(M)) {
    if (!is.numeric(eta.list[[k]])) { stop(paste0("eta.list[[", k, "]] must be numeric."), call. = FALSE) }
  }

  # -- Build eta grid --
  eta.grid <- as.matrix(expand.grid(eta.list))
  colnames(eta.grid) <- paste0("eta_", seq_len(M))
  if (nrow(eta.grid) > 1000) {
    warning("Large eta grid (", nrow(eta.grid), " combinations). ",
            "Consider coarser eta.list grids.", call. = FALSE)
  }
  eta.grid.list <- lapply(seq_len(nrow(eta.grid)), function(i) eta.grid[i, ])
  if (trace) { cat("Total eta combinations:", nrow(eta.grid), "\n") }

  # -- Prepare fit args --
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$family <- family
  fit.args$eta <- 0
  fit.args$y.external <- NULL

  ## Calculate null deviance (target cohort only)
  penalty.factor <- if (!is.null(fit.args$penalty.factor)) {
    fit.args$penalty.factor
  } else {
    rep(1, ncol(X))
  }
  target_idx <- cohort == 0
  n_target <- sum(target_idx)
  null.dev <- calcNullDev(
    X[target_idx, , drop = FALSE], y[target_idx, , drop = FALSE],
    rep(1, n_target) / n_target, penalty.factor, family
  )

  # -- Single eta fit --
  single_eta_fit <- function(eta_row, fit.args, cohort, ext_cohorts, trace) {
    if (trace) { cat("Fitting at eta = (", paste(round(eta_row, 3), collapse = ", "), ")\n") }

    # build per-observation weights: target = 1, external cohort k = eta_k
    w <- rep(1.0, length(cohort))
    for (k in seq_along(ext_cohorts)) {
      w[cohort == ext_cohorts[k]] <- eta_row[k]
    }

    fit.args$weights <- w
    fit_full <- do.call(BRIERi_fit, fit.args)

    # store eta on the fit object for consistency with BRIERi output
    fit_full$eta <- as.numeric(eta_row)
    fit_full
  }

  # -- Parallel / serial loop --
  use_parallel <- isTRUE(parallel) && ncores > 1L && .Platform$OS.type != "windows"
  if (.Platform$OS.type == "windows" && ncores > 1L) {
    message("Parallel execution not supported on Windows; using serial.")
  }

  if (!use_parallel) {
    res <- lapply(eta.grid.list, FUN = single_eta_fit,
                  fit.args = fit.args, cohort = cohort,
                  ext_cohorts = ext_cohorts, trace = trace)
  } else {
    res <- parallel::mclapply(
      eta.grid.list, FUN = single_eta_fit,
      fit.args = fit.args, cohort = cohort,
      ext_cohorts = ext_cohorts, trace = trace,
      mc.cores = ncores
    )
  }

  # -- Output: same structure as BRIERi --
  out <- list(
    y          = y,
    y.external = NULL,
    family     = family,
    cohort     = cohort,
    eta.list   = eta.list,
    eta.grid   = eta.grid,
    res        = res,
    null.dev   = null.dev
  )
  class(out) <- "BRIER"
  out
}