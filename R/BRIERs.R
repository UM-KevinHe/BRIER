#' BRIER.S: integration of pretrained external models with GWAS summary statistics
#'
#' Fits a penalized regression model that integrates one or more pretrained
#' external prediction models with target-cohort GWAS summary statistics and an
#' LD reference panel. External information is blended into the marginal
#' correlation vector via a tunable eta parameter (one per external model), and
#' lambda is selected over a regularisation path.
#'
#' @param sumstats A data.frame of GWAS summary statistics with at minimum a
#'   \code{corr} column (marginal correlation between genotype and outcome).
#'   Optionally \code{CHR} and \code{BP} columns for variant naming, or a
#'   \code{varnames} column for explicit names.
#' @param XtX A p x p sparse matrix representing the LD reference panel.
#'   Will be coerced to \code{Matrix::sparseMatrix} if not already sparse.
#' @param family A string specifying the response distribution: "gaussian",
#'   "binomial", or "poisson".
#' @param eta.list A list of numeric vectors (one per external model) specifying
#'   the eta grid for that model. A single numeric vector is replicated for
#'   all external models with a warning.
#' @param beta.external A numeric matrix of external model coefficients
#'   (p x M, no intercept). Defaults to a vector of zeros.
#' @param multi.method A string specifying how to combine multiple external
#'   models: "ind", "PCA", or "stacking". See \code{\link{BRIERi}}.
#' @param optim.args A list of arguments passed to \code{optim()} when using
#'   stacking.
#' @param ... Additional arguments passed to \code{BRIERs_fit}.
#' @param trace Logical. If TRUE, prints progress messages during fitting.
#' @param ncores Integer. Number of cores for parallel fitting.
#' @param parallel Logical. If TRUE and on a non-Windows platform, fits the eta
#'   grid in parallel using \code{parallel::mclapply}.
#'
#' @return An object of class \code{"BRIER"} containing:
#' \describe{
#'   \item{XtY}{The marginal correlation vector.}
#'   \item{XtY.external}{Computed external summary predictions.}
#'   \item{beta.external}{The external coefficients (possibly aggregated).}
#'   \item{family}{The response family.}
#'   \item{eta.list}{The list of per-model eta grids.}
#'   \item{eta.grid}{The full combinatorial eta grid.}
#'   \item{res}{A list of \code{BRIER.fit} objects, one per eta combination.}
#'   \item{null.dev}{Placeholder set to 0 (null deviance not defined for
#'     summary-stat models).}
#'   \item{varnames}{The variant names.}
#' }
#'
#' @seealso \code{\link{BRIERs_fit}}, \code{\link{BRIERs.selection}},
#'   \code{\link{BRIERi}}, \code{\link{BRIERfull}}
#'
#' @examples
#' \dontrun{
#' # sumstats: data.frame with 'corr', 'CHR', 'BP' columns
#' # XtX: sparse p x p LD matrix
#' # beta.external: p x M matrix of external coefficients
#'
#' fit <- BRIERs(
#'   sumstats, XtX, family = "gaussian",
#'   eta.list = list(c(0, 0.5, 1, 2)),
#'   beta.external = beta.external,
#'   penalty = "LASSO"
#' )
#' }
#'
#' @export
BRIERs <- function(
  sumstats, XtX, family = c("gaussian", "binomial", "poisson"),
  eta.list, beta.external = rep(0, nrow(sumstats)),
  multi.method = c("ind", "PCA", "stacking"), optim.args = list(),
  ...,
  trace = FALSE,
  ncores = max(1L, parallel::detectCores() - 1L),
  parallel = (ncores > 1L)
) {

  family <- match.arg(family)
  multi.method <- match.arg(multi.method)

  # -- Validate GWAS summary --
  if (is.matrix(sumstats)) { sumstats <- as.data.frame(sumstats) }
  if (!is.data.frame(sumstats)) { stop("sumstats must be a data.frame.", call. = FALSE) }
  if (is.null(sumstats$corr)) { stop("sumstats must contain a 'corr' column.", call. = FALSE) }
  XtY <- as.numeric(as.vector(sumstats$corr))
  if (any(!is.finite(XtY))) { stop("sumstats$corr must be finite (no NA/Inf).", call. = FALSE) }
  p <- length(XtY)

  ## Build variable names
  chr <- sumstats$CHR
  pos <- sumstats$BP
  varnames <- if (!is.null(sumstats$varnames)) {
    sumstats$varnames
  } else if (!is.null(chr) && !is.null(pos)) {
    paste(chr, pos, sep = ":")
  } else {
    paste0("X", seq_len(p))
  }

  # -- Validate XtX --
  if (!inherits(XtX, "Matrix")) { XtX <- Matrix::Matrix(XtX, sparse = TRUE) }
  if (nrow(XtX) != p || ncol(XtX) != p) {
    stop(paste0("XtX dimensions (", nrow(XtX), " x ", ncol(XtX),
                ") do not match sumstats (", p, " variants)."), call. = FALSE)
  }

  # -- Validate beta.external --
  beta.external <- as.matrix(beta.external)
  if (any(is.na(beta.external))) { stop("Missing data (NA's) detected in beta.external.", call. = FALSE) }
  if (nrow(beta.external) != p) {
    stop(paste0(
      "beta.external must have ", p, " rows (no intercept). Got ", nrow(beta.external), "."
    ), call. = FALSE)
  }

  # -- Compute external summary predictions --
  ext <- calcExtXtY(XtX, XtY, beta.external, multi.method)
  XtY.external <- ext$XtY.external
  M <- ncol(XtY.external)

  # -- Validate eta.list --
  if (!is.list(eta.list)) {
    if (!is.numeric(eta.list)) { stop("eta.list must be a numeric vector or a list of numeric vectors.", call. = FALSE) }
    warning("eta.list is not a list. Replicating the same eta grid for all ",
            M, " external model(s).", call. = FALSE)
    eta.list <- rep(list(eta.list), M)
  }
  if (length(eta.list) != M) {
    stop(paste0("eta.list has ", length(eta.list), " element(s) but there are ",
                M, " external model(s)."), call. = FALSE)
  }
  for (k in seq_len(M)) {
    if (!is.numeric(eta.list[[k]])) { stop(paste0("eta.list[[", k, "]] must be numeric."), call. = FALSE) }
  }
  if (M >= 5) {
    warning("5 or more external models provided. Please consider using ensemble ",
            "weighting to aggregate multiple external models.", call. = FALSE)
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
  fit.args$XtX <- XtX
  fit.args$XtY <- XtY
  fit.args$XtY.external <- XtY.external
  fit.args$family <- family
  fit.args$varnames <- varnames

  # -- Single eta fit --
  single_eta_fit <- function(eta_row, fit.args, trace) {
    if (trace) { cat("Fitting at eta = (", paste(round(eta_row, 3), collapse = ", "), ")\n") }
    fit.args$eta <- as.numeric(eta_row)
    do.call(BRIERs_fit, fit.args)
  }

  # -- Parallel / serial loop --
  use_parallel <- isTRUE(parallel) && ncores > 1L && .Platform$OS.type != "windows"
  if (.Platform$OS.type == "windows" && ncores > 1L) {
    message("Parallel execution not supported on Windows; using serial.")
  }

  if (!use_parallel) {
    res <- lapply(eta.grid.list, FUN = single_eta_fit, fit.args = fit.args, trace = trace)
  } else {
    res <- parallel::mclapply(
      eta.grid.list, FUN = single_eta_fit,
      fit.args = fit.args, trace = trace,
      mc.cores = ncores
    )
  }

  # -- Output: same structure as BRIERi --
  out <- list(
    XtY            = XtY,
    XtY.external   = XtY.external,
    beta.external  = beta.external,
    family         = family,
    eta.list       = eta.list,
    eta.grid       = eta.grid,
    res            = res,
    null.dev = 0,
    varnames       = varnames
  )
  class(out) <- "BRIER"
  out
}


#' Core fitting function for BRIER.S (summary statistics)
#'
#' Internal workhorse that fits a penalized regression model for a single eta
#' vector using GWAS summary statistics and an LD reference panel. Called by
#' \code{\link{BRIERs}}. No intercept is included since summary statistics
#' have no observation-level structure.
#'
#' @param XtX A p x p sparse LD matrix.
#' @param XtY A numeric vector of marginal correlations (length p).
#' @param eta A numeric vector of external-model weights. Use \code{eta = 0}
#'   to ignore external information.
#' @param XtY.external Optional p x M matrix of external summary predictions.
#' @param family A string: "gaussian", "binomial", or "poisson".
#' @param penalty A string: "LASSO", "MCP", or "SCAD".
#' @param penalty.factor A numeric vector of per-variable penalty multipliers
#'   (length p).
#' @param alpha A numeric scalar for the elastic net mixing parameter.
#' @param gamma A numeric scalar for MCP and SCAD penalties.
#' @param lambda Optional user-specified lambda sequence.
#' @param nlambda An integer specifying the number of lambda values.
#' @param lambda.min A numeric scalar for the minimum lambda ratio.
#' @param log.lambda Logical. If TRUE, the lambda sequence is log-spaced.
#' @param max.iter An integer specifying the maximum total iterations.
#' @param eps A numeric scalar for the convergence tolerance.
#' @param nvar.max An integer specifying the maximum number of selected variables.
#' @param varnames Optional character vector of variant names (length p).
#'
#' @return An object of class \code{"BRIER.fit"} containing model specification,
#'   coefficients, fit statistics, and data. The \code{summary = TRUE} flag
#'   indicates a summary-stat model (no intercept).
#'
#' @seealso \code{\link{BRIERs}}, \code{\link{coef.BRIER.fit}},
#'   \code{\link{predict.BRIER.fit}}
#'
#' @keywords internal
#' @export
BRIERs_fit <- function(
  XtX, XtY, eta = 0, XtY.external = NULL,
  family = c("gaussian", "binomial", "poisson"),
  penalty = c("LASSO", "MCP", "SCAD"),
  penalty.factor = rep(1, length(XtY)),
  alpha = 1, gamma = ifelse(penalty == "SCAD", 3.7, 3),
  lambda, nlambda = 100, lambda.min = 0.001, log.lambda = TRUE,
  max.iter = 1e6, eps = 1e-4, nvar.max = NULL,
  varnames = NULL
) {

  penalty <- match.arg(penalty)
  family <- match.arg(family)

  XtY <- as.numeric(as.vector(XtY))
  p <- length(XtY)

  # -- Validate eta --
  eta <- as.numeric(as.vector(eta))
  if (any(!is.finite(eta))) { stop("eta must be finite (no NA/Inf).", call. = FALSE) }
  if (any(eta < 0)) { stop("eta must be >= 0.", call. = FALSE) }

  # -- Validate XtY.external --
  if (!is.null(XtY.external)) {
    XtY.external <- as.matrix(XtY.external)
    if (nrow(XtY.external) != p) {
      stop(paste0("XtY.external must have ", p, " rows. Got ", nrow(XtY.external), "."), call. = FALSE)
    }
    if (length(eta) != ncol(XtY.external)) {
      stop("Length of eta must match number of columns in XtY.external.", call. = FALSE)
    }
  } else {
    if (any(eta != 0)) { stop("XtY.external must be provided when any eta != 0.", call. = FALSE) }
  }

  # -- Validate dimensions --
  if (!inherits(XtX, "Matrix")) { XtX <- Matrix::Matrix(XtX, sparse = TRUE) }
  if (nrow(XtX) != p || ncol(XtX) != p) {
    stop(paste0("XtX must be ", p, " x ", p, "."), call. = FALSE)
  }

  # -- Validate penalty.factor --
  if (!is.double(penalty.factor)) { penalty.factor <- as.double(penalty.factor) }
  if (any(is.na(penalty.factor))) { stop("Missing data (NA's) detected in penalty.factor.", call. = FALSE) }
  if (length(penalty.factor) != p) { stop("Length of penalty.factor must equal length of XtY.", call. = FALSE) }
  if (any(penalty.factor < 0)) { stop("penalty.factor must be non-negative.", call. = FALSE) }

  # -- Checking sparse parameters --
  if (nlambda < 2) {
    stop("nlambda must be at least 2", call. = FALSE)
  } else if (nlambda != round(nlambda)){
    stop("nlambda must be a positive integer", call. = FALSE)
  }
  if (alpha <= 0 || alpha > 1)  { stop("alpha must be between (0, 1]; choose a small positive number instead", call. = FALSE) }
  if (gamma <= 1 && penalty == "MCP") { stop("gamma must be greater than 1 for MCP.", call. = FALSE) }
  if (gamma <= 2 && penalty == "SCAD") { stop("gamma must be greater than 2 for SCAD.", call. = FALSE) }
  if (is.null(nvar.max)) { nvar.max <- p }

  # -- Variable names --
  if (is.null(varnames)) { varnames <- paste0("X", seq_len(p)) }

  # -- Compute summary XtY.tilde --
  if (all(eta == 0) || is.null(XtY.external)) {
    XtY.tilde <- XtY
  } else {
    XtY.tilde <- (XtY + XtY.external %*% eta) / (1 + sum(eta))
  }

  # -- Lambda sequence --
  if (missing(lambda)) {

    lambda.max <- maxlambda_summary(
      XtY.tilde, XtX,
      penalty.factor, alpha, eps, max.iter
    )

    if (log.lambda) {
      if (lambda.min == 0) {
        lambda.seq <- c(exp(seq(log(lambda.max), log(0.001 * lambda.max), length = nlambda - 1)), 0)
      } else {
        lambda.seq <- exp(seq(log(lambda.max), log(lambda.min * lambda.max), length = nlambda))
      }
    } else {
      if (lambda.min == 0) {
        lambda.seq <- c(seq(lambda.max, 0.001 * lambda.max, length = nlambda - 1), 0)
      } else {
        lambda.seq <- seq(lambda.max, lambda.min * lambda.max, length = nlambda)
      }
    }
    user <- FALSE
  } else {
    nlambda <- length(lambda)
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
    user <- TRUE
  }

  # -- Fit --
  res <- cd_summary_fit(
    XtY.tilde, XtX, 
    penalty.factor, penalty, lambda.seq,
    alpha, gamma, eps, max.iter, nvar.max, user
  )

  b      <- res$beta
  iter   <- res$iter
  df     <- res$df

  # -- Eliminate saturated lambdas --
  ind <- !is.na(iter) & colSums(!is.finite(b)) == 0
  b      <- b[, ind, drop = FALSE]
  lambda <- lambda.seq[ind]
  df     <- df[ind]
  iter   <- iter[ind]

  if (length(iter) == 0 || iter[1] == max.iter) {
    stop("Algorithm failed to converge for any values of lambda.", call. = FALSE)
  }
  if (any(iter == max.iter)) {
    warning("Algorithm failed to converge for some values of lambda.", call. = FALSE)
  }

  # -- Build full beta matrix --
  beta <- b
  dimnames(beta) <- list(varnames, round(lambda, 4))
  # deviance <- colSums(beta * (XtX %*% beta)) - 2 * crossprod(beta, XtY)[, 1]
  XtXbeta <- as.matrix(XtX %*% beta)
  beta_mat <- as.matrix(beta)
  deviance <- colSums(beta_mat * XtXbeta) - 2 * as.numeric(crossprod(beta_mat, XtY))

  # -- Degree of freedom --
  k <- colSums(abs(beta) >= eps)

  out <- list(
    # Model specification
    family         = family,
    penalty        = penalty,
    penalty.factor = penalty.factor,
    eta            = eta,
    lambda         = lambda,
    alpha          = alpha,
    gamma          = gamma,

    # Coefficients
    beta           = beta,
    df             = df,
    k              = k,

    # Fit statistics
    deviance       = deviance,
    iter           = iter,

    # Data
    XtY = XtY, 
    XtY.tilde = XtY.tilde,
    XtY.external = XtY.external,
    summary        = TRUE
  )
  class(out) <- "BRIER.fit"
  out
}

# -- Internal helper: aggregate external summary predictions --

calcExtXtY <- function(XtX, XtY, beta.external, multi.method) {
  if (multi.method == "PCA") {
    bb <- apply(beta.external, 2, function(x) x / sqrt(sum(x^2)))
    w <- prcomp(t(bb) %*% bb)$rotation[, 1]
    beta.external <- beta.external %*% abs(w)
  } else if (multi.method == "stacking") {
    w <- solve(crossprod(beta.external, XtX %*% beta.external), crossprod(beta.external, XtY))
    beta.external <- beta.external %*% w
  }
  
  XtY.external <- as.matrix(XtX %*% beta.external)
  list(beta.external = beta.external, XtY.external = XtY.external)
}




