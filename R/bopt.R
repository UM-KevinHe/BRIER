# =============================================================================
# Bayesian optimization of the integration weight eta for the three BRIER
# modules. The grid fitters (BRIERfull, BRIERi, BRIERs) evaluate every eta in
# expand.grid(eta.list), which grows multiplicatively with the number of
# external sources M. These functions replace that exhaustive grid with a
# sequential model-based search over a continuous eta box, keeping the inner
# lambda path and the selection criterion exactly as the grid fitters use them.
#
# Each evaluation refits at a single eta and scores it with the SAME internal
# helper the corresponding *.selection() function uses, so a bopt result and a
# grid result are directly comparable. The returned object is a standard
# BRIER.selection object, so coef(), predict(), plot.eta(), plot.box() and
# plot.importance() all work on it unchanged.
# =============================================================================


# -- Internal: validate the eta box and the initial design --

.bopt_bounds <- function(bounds, M, label) {

  if (is.null(bounds)) {
    bounds <- rep(list(c(0, 10)), M)
    names(bounds) <- paste0("eta_", seq_len(M))
  }
  if (!is.list(bounds)) { stop("bounds must be a named list of length-2 numeric vectors.", call. = FALSE) }
  if (length(bounds) != M) {
    stop(paste0(
      "bounds has ", length(bounds), " element(s) but there are ", M, " ", label, "."
    ), call. = FALSE)
  }
  if (is.null(names(bounds)) || any(names(bounds) == "")) {
    names(bounds) <- paste0("eta_", seq_len(M))
  }
  for (k in seq_len(M)) {
    b <- bounds[[k]]
    if (!is.numeric(b) || length(b) != 2L || any(!is.finite(b))) {
      stop(paste0("bounds[[", k, "]] must be a finite numeric vector of length 2."), call. = FALSE)
    }
    if (b[1] < 0) { stop("eta must be >= 0, so the lower bound cannot be negative.", call. = FALSE) }
    if (b[1] >= b[2]) { stop(paste0("bounds[[", k, "]] must be increasing (lower < upper)."), call. = FALSE) }
    bounds[[k]] <- as.numeric(b)
  }
  bounds
}


.bopt_init_grid <- function(init.grid, bounds, M) {

  if (is.null(init.grid)) { return(NULL) }

  if (is.numeric(init.grid) && is.null(dim(init.grid))) {
    if (M != 1L) {
      stop("init.grid can only be a plain numeric vector when there is a single external source.", call. = FALSE)
    }
    init.grid <- data.frame(eta = init.grid)
  }
  init.grid <- as.data.frame(init.grid)
  if (ncol(init.grid) != M) {
    stop(paste0("init.grid must have ", M, " column(s), one per external source. Got ", ncol(init.grid), "."), call. = FALSE)
  }
  colnames(init.grid) <- names(bounds)

  for (k in seq_len(M)) {
    col <- init.grid[[k]]
    if (!is.numeric(col) || any(!is.finite(col))) {
      stop(paste0("init.grid column '", names(bounds)[k], "' must be finite and numeric."), call. = FALSE)
    }
    if (any(col < bounds[[k]][1]) || any(col > bounds[[k]][2])) {
      stop(paste0(
        "init.grid column '", names(bounds)[k], "' falls outside bounds [",
        bounds[[k]][1], ", ", bounds[[k]][2], "]."
      ), call. = FALSE)
    }
  }
  init.grid
}


# -- Internal: resolve the optimizer --
#
# The dependency is reached through this one accessor rather than called
# directly, so the surrounding loop can be tested without rBayesianOptimization
# installed. The availability check lives here for the same reason.

.bopt_optimizer <- function() {
  if (!requireNamespace("rBayesianOptimization", quietly = TRUE)) {
    stop(
      "Bayesian optimization requires the 'rBayesianOptimization' package. ",
      "Install it with install.packages('rBayesianOptimization').",
      call. = FALSE
    )
  }
  rBayesianOptimization::BayesianOptimization
}


# -- Internal: the optimizer loop shared by all three modules --
#
# eval.eta is a closure taking a numeric eta vector of length M and returning
# list(fit = <BRIER.eta>, row = <one-row eta.lambda data.frame>). Every fit is
# cached here rather than routed through the optimizer's own bookkeeping, so
# the heavy fitted objects are never copied into the acquisition machinery.

.bopt_engine <- function(
  eval.eta, M, bounds,
  init.grid, init.points, n.iter,
  acq, kappa, acq.eps, kernel, verbose
) {

  optimize <- .bopt_optimizer()

  if (is.null(init.grid) && init.points < 2L) {
    stop(
      "Supply either an init.grid or init.points >= 2: the surrogate model ",
      "cannot be fitted from fewer than two evaluated points.",
      call. = FALSE
    )
  }
  if (n.iter < 1L) { stop("n.iter must be at least 1.", call. = FALSE) }

  arg.names <- names(bounds)

  cache <- new.env(parent = emptyenv())
  cache$fits <- list()
  cache$rows <- list()
  cache$etas <- list()
  cache$worst <- NA_real_
  cache$failed <- list()

  objective <- function() {

    eta <- as.numeric(unlist(mget(arg.names, envir = environment(), inherits = FALSE)))

    ev <- try(eval.eta(eta), silent = TRUE)

    if (inherits(ev, "try-error")) {
      cache$failed[[length(cache$failed) + 1L]] <- eta
      warning(
        "Fit failed at eta = (", paste(signif(eta, 4), collapse = ", "),
        "); scoring it as the worst point seen so far. Message: ",
        conditionMessage(attr(ev, "condition")),
        call. = FALSE
      )
      penalised <- if (is.na(cache$worst)) -1e6 else cache$worst - max(1, abs(cache$worst) * 0.1)
      return(list(Score = penalised, Pred = NA_real_))
    }

    score <- -as.numeric(ev$row$measure.min)
    if (!is.finite(score)) {
      cache$failed[[length(cache$failed) + 1L]] <- eta
      warning(
        "Non-finite criterion at eta = (", paste(signif(eta, 4), collapse = ", "), ").",
        call. = FALSE
      )
      penalised <- if (is.na(cache$worst)) -1e6 else cache$worst - max(1, abs(cache$worst) * 0.1)
      return(list(Score = penalised, Pred = NA_real_))
    }

    k <- length(cache$fits) + 1L
    cache$fits[[k]] <- ev$fit
    cache$rows[[k]] <- ev$row
    cache$etas[[k]] <- eta
    cache$worst <- if (is.na(cache$worst)) score else min(cache$worst, score)

    list(Score = score, Pred = k)
  }
  formals(objective) <- stats::setNames(rep(list(quote(expr = )), M), arg.names)

  opt <- optimize(
    FUN          = objective,
    bounds       = bounds,
    init_grid_dt = init.grid,
    init_points  = init.points,
    n_iter       = n.iter,
    acq          = acq,
    kappa        = kappa,
    eps          = acq.eps,
    kernel       = kernel,
    verbose      = verbose
  )

  n.ok <- length(cache$fits)
  if (n.ok == 0L) {
    stop("Every eta evaluation failed; there is no model to return.", call. = FALSE)
  }

  # -- Assemble the evaluated points into a standard grid-shaped result --

  eta.grid <- do.call(rbind, cache$etas)
  colnames(eta.grid) <- paste0("eta_", seq_len(M))

  keep <- !duplicated(apply(round(eta.grid, 10), 1, paste, collapse = "|"))
  eta.grid <- eta.grid[keep, , drop = FALSE]
  fits <- cache$fits[keep]
  rows <- cache$rows[keep]

  ord <- do.call(order, as.data.frame(eta.grid))
  eta.grid <- eta.grid[ord, , drop = FALSE]
  fits <- fits[ord]
  rows <- rows[ord]

  eta.lambda <- do.call(rbind, rows)
  eta.lambda$eta.index <- seq_len(nrow(eta.lambda))
  rownames(eta.lambda) <- NULL

  eta.min.index <- which.min(eta.lambda$measure.min)
  eta.min <- eta.grid[eta.min.index, ]
  lambda.min.index <- eta.lambda$lambda.min.index[eta.min.index]
  lambda.min <- eta.lambda$lambda.min[eta.min.index]

  if (M == 1L) {
    cat("Best eta:", round(eta.min, 3), "\n")
  } else {
    cat("Best eta: (", paste(round(eta.min, 3), collapse = ", "), ")\n")
  }
  cat("Best lambda:", round(lambda.min, 3), "\n")
  cat("Evaluations:", nrow(eta.grid), "unique of", length(cache$fits), "fitted\n")

  list(
    res              = fits,
    eta.grid         = eta.grid,
    eta.lambda       = eta.lambda,
    eta.min          = eta.min,
    eta.min.index    = eta.min.index,
    lambda.min       = lambda.min,
    lambda.min.index = lambda.min.index,
    bopt             = list(
      bounds      = bounds,
      init.points = init.points,
      n.iter      = n.iter,
      acq         = acq,
      kappa       = kappa,
      eps         = acq.eps,
      kernel      = kernel,
      history     = opt$History,
      best.par    = opt$Best_Par,
      best.value  = opt$Best_Value,
      n.failed    = length(cache$failed),
      failed.eta  = if (length(cache$failed)) do.call(rbind, cache$failed) else NULL
    )
  )
}


# -- Internal: score one fitted eta with the module's own selection helper --

.bopt_score_i <- function(fit, eta, criteria, X.val, y.val, n, var.y, y) {
  object <- list(res = list(fit), eta.grid = matrix(eta, nrow = 1L), y = y)
  if (criteria %in% .bopt_val_criteria()) {
    validation(1L, object, X.val, y.val, criteria)
  } else {
    ic_selection(1L, object, n, criteria, dispersion = NULL, var_y = var.y)
  }
}

.bopt_score_s <- function(fit, eta, criteria, X.val, y.val, XtX, sumstats, TN, h2) {
  object <- list(res = list(fit), eta.grid = matrix(eta, nrow = 1L))
  if (criteria %in% .bopt_val_criteria()) {
    validation(1L, object, X.val, y.val, criteria)
  } else {
    ic_selection_S(1L, object, XtX, sumstats, TN, h2, criteria)
  }
}

.bopt_val_criteria <- function() {
  c(
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  )
}


#' Bayesian optimization of eta for BRIER.I
#'
#' Tunes the integration weight \code{eta} for \code{\link{BRIERi}} with
#' sequential model-based optimization instead of an exhaustive grid. Each
#' proposed \code{eta} triggers one penalized fit over the full lambda path,
#' scored with the same criterion \code{\link{BRIERi.selection}} would apply,
#' so the result is directly comparable to a grid search.
#'
#' This is most useful when several external models are integrated with
#' \code{multi.method = "ind"}, where the grid grows multiplicatively in the
#' number of sources and quickly becomes infeasible.
#'
#' @param X A numeric matrix of predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param family A string: "gaussian", "binomial", or "poisson".
#' @param beta.external A \code{(p + 1) x M} matrix of external coefficients
#'   with the intercept in the first row.
#' @param multi.method External aggregation strategy: "ind", "PCA", "stacking",
#'   or "PCstacking". Passed to \code{\link{BRIERi}}'s external reduction, which is
#'   computed once and reused across every evaluation.
#' @param optim.args A list of control arguments for the stacking optimizer.
#' @param criteria Selection criterion, as in \code{\link{BRIERi.selection}}.
#'   Information criteria ("gcv", "AIC", "BIC", "Cp") need no held-out data;
#'   validation criteria require \code{X.val} and \code{y.val}.
#' @param X.val,y.val Held-out validation data, required for validation criteria.
#' @param n,var.y Optional overrides for the sample size and response variance
#'   used by the information criteria.
#' @param bounds A named list of length-2 numeric vectors giving the search box
#'   for each eta, for example \code{list(eta_1 = c(0, 10))}. Defaults to
#'   \code{c(0, 10)} per external model. The lower bound cannot be negative.
#' @param init.grid Optional initial design: a data.frame with one column per
#'   external model, or a plain numeric vector when M is 1. Evaluated before
#'   the optimizer starts.
#' @param init.points Number of random initial points drawn when
#'   \code{init.grid} is not supplied. At least two points are needed to fit
#'   the surrogate.
#' @param n.iter Number of sequential optimization steps after initialization.
#' @param acq Acquisition function: "ucb", "ei", or "poi".
#' @param kappa Exploration parameter for "ucb".
#' @param acq.eps Exploration parameter for "ei" and "poi". Named to avoid a
#'   collision with the coordinate-descent tolerance \code{eps}, which is
#'   passed through \code{...} to the fitter.
#' @param kernel Gaussian-process kernel specification.
#' @param verbose Logical. Print optimizer progress.
#' @param pca.var Fraction of the variance retained by the principal-component
#'   reduction under \code{multi.method = "PCstacking"}. Ignored otherwise.
#' @param dedup.cor Similarity threshold for dropping redundant external models
#'   before the search starts, keeping the first occurrence of each
#'   near-identical set. When a column is dropped, \code{bounds} and
#'   \code{init.grid} sized to the supplied panel are subset to match.
#'   \code{NULL} or \code{NA} disables the step. See \code{\link{dedupExternals}}.
#' @param ... Further arguments passed to \code{\link{BRIERi.eta}}, for example
#'   \code{penalty}, \code{alpha}, \code{gamma}, \code{penalty.factor},
#'   \code{nlambda}, or \code{eps}.
#'
#' @return An object of class \code{c("BRIER.bopt", "BRIER.selection", "BRIER")}.
#'   It carries the usual selection fields (\code{eta.min}, \code{lambda.min},
#'   \code{eta.lambda}, \code{res}, \code{eta.grid}) built from the evaluated
#'   points sorted by eta, so \code{coef}, \code{predict} and the \code{plot.*}
#'   functions work on it directly. The extra \code{bopt} element records the
#'   search: bounds, settings, optimizer history, and any failed evaluations.
#'
#' @seealso \code{\link{BRIERi}}, \code{\link{BRIERi.selection}},
#'   \code{\link{BRIERs.bopt}}, \code{\link{BRIERfull.bopt}}
#'
#' @examples
#' \dontrun{
#' sel <- BRIERi.bopt(
#'   X, y,
#'   family        = "gaussian",
#'   beta.external = beta.external,
#'   criteria      = "gaussian.mspe",
#'   X.val         = X.val,
#'   y.val         = y.val,
#'   bounds        = list(eta_1 = c(0, 10)),
#'   init.points   = 10,
#'   n.iter        = 15,
#'   penalty       = "LASSO"
#' )
#' coef(sel)
#' }
#'
#' @export
BRIERi.bopt <- function(
  X, y, family = c("gaussian", "binomial", "poisson"),
  beta.external = rep(0, ncol(X) + 1),
  multi.method = c("ind", "PCA", "stacking", "PCstacking"), optim.args = list(),
  criteria = c(
    "gcv", "AIC", "BIC", "Cp",
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  X.val = NULL, y.val = NULL, n = NULL, var.y = NULL,
  bounds = NULL, init.grid = NULL, init.points = 10L, n.iter = 15L,
  acq = "ucb", kappa = 2.576, acq.eps = 0,
  kernel = list(type = "exponential", power = 2),
  verbose = TRUE,
  pca.var = 0.8,
  dedup.cor = 0.9,
  ...
) {

  family <- match.arg(family)
  multi.method <- match.arg(multi.method)
  criteria <- match.arg(criteria)

  # -- Validation of y (mirrors BRIERi) --
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

  # -- Validation of external info --
  beta.external <- as.matrix(beta.external)
  if (any(is.na(beta.external))) { stop("Missing data (NA's) detected in beta.external.", call. = FALSE) }
  if (nrow(beta.external) != ncol(X) + 1) {
    stop("The dimension of external beta and X does not match. Please include intercept.", call. = FALSE)
  }

  if (criteria %in% .bopt_val_criteria()) {
    if (is.null(X.val) || is.null(y.val)) {
      stop("X.val and y.val are required for validation-based criteria.", call. = FALSE)
    }
    X.val <- as.matrix(X.val)
    y.val <- as.numeric(as.vector(y.val))
    if (nrow(X.val) != length(y.val)) {
      stop("X.val and y.val do not have the same number of observations.", call. = FALSE)
    }
  }

  ## Redundant sources go first: with "ind" each duplicate would otherwise add a
  ## search dimension the optimizer has to spend evaluations on.
  dedup <- .dedup_and_report(beta.external, dedup.cor, intercept.row = TRUE)
  beta.external <- dedup$beta.external
  if (dedup$applied && !is.null(bounds) && length(bounds) == dedup$M.in) {
    ## Names are dropped so the survivors are renumbered eta_1 .. eta_M rather
    ## than keeping the gaps of the panel they came from, which would disagree
    ## with the column names of eta.grid.
    bounds <- unname(bounds[dedup$keep])
  }
  if (!is.null(init.grid) && !is.null(dim(init.grid)) &&
      ncol(init.grid) == dedup$M.in && dedup$applied) {
    init.grid <- init.grid[, dedup$keep, drop = FALSE]
  }

  # -- External reduction and null deviance: computed once, reused every step --
  ext <- calcExtY(X, y, beta.external, family, multi.method, optim.args, pca.var)
  y.external <- ext$y.external
  M <- ncol(y.external)

  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$family <- family
  fit.args$y.external <- y.external

  penalty.factor <- if (!is.null(fit.args$penalty.factor)) fit.args$penalty.factor else rep(1, ncol(X))
  null.dev <- calcNullDev(X, y, rep(1, nrow(X)) / nrow(X), penalty.factor, family)

  n <- if (!is.null(n)) n else length(y)
  var.y <- if (!is.null(var.y)) var.y else var(y)

  bounds <- .bopt_bounds(bounds, M, "external model(s)")
  init.grid <- .bopt_init_grid(init.grid, bounds, M)

  eval.eta <- function(eta) {
    args <- fit.args
    args$eta <- as.numeric(eta)
    fit <- do.call(BRIERi.eta, args)
    row <- .bopt_score_i(fit, eta, criteria, X.val, y.val, n, var.y, y)
    list(fit = fit, row = row)
  }

  opt <- .bopt_engine(
    eval.eta, M, bounds, init.grid, init.points, n.iter,
    acq, kappa, acq.eps, kernel, verbose
  )

  out <- list(
    y                = y,
    y.external       = y.external,
    beta.external    = beta.external,
    family           = family,
    eta.list         = lapply(seq_len(M), function(k) sort(unique(opt$eta.grid[, k]))),
    eta.grid         = opt$eta.grid,
    res              = opt$res,
    null.dev         = null.dev,
    n                = nrow(X),
    p                = ncol(X),
    M                = M,
    external.dedup   = if (dedup$applied) dedup else NULL,
    external.pca     = ext$external.pca,
    criteria         = criteria,
    eta.min          = opt$eta.min,
    eta.min.index    = opt$eta.min.index,
    lambda.min       = opt$lambda.min,
    lambda.min.index = opt$lambda.min.index,
    eta.lambda       = opt$eta.lambda,
    bopt             = opt$bopt
  )
  class(out) <- c("BRIER.bopt", "BRIER.selection", "BRIER")
  out
}


#' Bayesian optimization of eta for BRIER.S
#'
#' Tunes the integration weight \code{eta} for \code{\link{BRIERs}} with
#' sequential model-based optimization instead of an exhaustive grid. Each
#' proposed \code{eta} triggers one summary-level fit over the full lambda
#' path, scored with the same criterion \code{\link{BRIERs.selection}} would
#' apply.
#'
#' @param sumstats A data.frame of target summary statistics containing a
#'   \code{corr} column.
#' @param XtX A p x p sparse LD matrix aligned to \code{sumstats}.
#' @param family A string: "gaussian", "binomial", or "poisson".
#' @param beta.external A \code{p x M} matrix of external coefficients, with no
#'   intercept row.
#' @param multi.method External aggregation strategy: "ind", "PCA", "stacking",
#'   or "PCstacking".
#' @param optim.args A list of control arguments for the stacking optimizer.
#' @param criteria Selection criterion, as in \code{\link{BRIERs.selection}}.
#'   Summary criteria are "Cp", "GIC" (need \code{TN}) and "pseu.val";
#'   validation criteria require \code{X.val} and \code{y.val}.
#' @param X.val,y.val Held-out validation data, required for validation
#'   criteria. \code{X.val} must already be standardized, and \code{y.val} must
#'   also be standardized when \code{family = "gaussian"}, because BRIER.S
#'   returns standardized coefficients and nothing is standardized internally.
#' @param XtX.val,sumstats.val Optional held-out LD matrix and summary
#'   statistics used by \code{"pseu.val"}. Default to the fitting \code{XtX}
#'   and \code{sumstats}.
#' @param TN Training cohort sample size, required for "Cp" and "GIC".
#' @param h2 Optional heritability estimate for "Cp" and "GIC". Defaults to 0.
#' @param bounds A named list of length-2 numeric vectors giving the search box
#'   for each eta. Defaults to \code{c(0, 10)} per external model.
#' @param init.grid Optional initial design, as in \code{\link{BRIERi.bopt}}.
#' @param init.points Number of random initial points when \code{init.grid} is
#'   not supplied.
#' @param n.iter Number of sequential optimization steps after initialization.
#' @param acq Acquisition function: "ucb", "ei", or "poi".
#' @param kappa Exploration parameter for "ucb".
#' @param acq.eps Exploration parameter for "ei" and "poi".
#' @param kernel Gaussian-process kernel specification.
#' @param verbose Logical. Print optimizer progress.
#' @param pca.var Fraction of the variance retained by the principal-component
#'   reduction under \code{multi.method = "PCstacking"}. Ignored otherwise.
#' @param dedup.cor Similarity threshold for dropping redundant external models
#'   before the search starts, keeping the first occurrence of each
#'   near-identical set. \code{NULL} or \code{NA} disables the step. See
#'   \code{\link{dedupExternals}}.
#' @param ... Further arguments passed to \code{\link{BRIERs.eta}}.
#'
#' @return An object of class \code{c("BRIER.bopt", "BRIER.selection", "BRIER")},
#'   shaped exactly like a \code{\link{BRIERs.selection}} result with an added
#'   \code{bopt} element recording the search.
#'
#' @seealso \code{\link{BRIERs}}, \code{\link{BRIERs.selection}},
#'   \code{\link{BRIERi.bopt}}, \code{\link{BRIERfull.bopt}}
#'
#' @examples
#' \dontrun{
#' sel <- BRIERs.bopt(
#'   sumstats, ld$XtX,
#'   family        = "gaussian",
#'   beta.external = beta.external,
#'   criteria      = "gaussian.mspe",
#'   X.val         = standardize_X(X.val)[[1]],
#'   y.val         = standardize_X(as.matrix(y.val))[[1]],
#'   init.points   = 10,
#'   n.iter        = 15
#' )
#' }
#'
#' @export
BRIERs.bopt <- function(
  sumstats, XtX, family = c("gaussian", "binomial", "poisson"),
  beta.external = rep(0, nrow(sumstats)),
  multi.method = c("ind", "PCA", "stacking", "PCstacking"), optim.args = list(),
  criteria = c(
    "Cp", "GIC", "pseu.val",
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  X.val = NULL, y.val = NULL,
  XtX.val = NULL, sumstats.val = NULL, TN = NULL, h2 = NULL,
  bounds = NULL, init.grid = NULL, init.points = 10L, n.iter = 15L,
  acq = "ucb", kappa = 2.576, acq.eps = 0,
  kernel = list(type = "exponential", power = 2),
  verbose = TRUE,
  pca.var = 0.8,
  dedup.cor = 0.9,
  ...
) {

  family <- match.arg(family)
  multi.method <- match.arg(multi.method)
  criteria <- match.arg(criteria)

  # -- Validate summary statistics (mirrors BRIERs) --
  if (is.matrix(sumstats)) { sumstats <- as.data.frame(sumstats) }
  if (!is.data.frame(sumstats)) { stop("sumstats must be a data.frame.", call. = FALSE) }
  if (is.null(sumstats$corr)) { stop("sumstats must contain a 'corr' column.", call. = FALSE) }
  XtY <- as.numeric(as.vector(sumstats$corr))
  if (any(!is.finite(XtY))) { stop("sumstats$corr must be finite (no NA/Inf).", call. = FALSE) }
  p <- length(XtY)

  chr <- sumstats$CHR
  pos <- sumstats$BP
  ref <- sumstats$REF
  alt <- sumstats$ALT
  varnames <- if (!is.null(sumstats$varnames)) {
    sumstats$varnames
  } else if (!is.null(chr) && !is.null(pos) && !is.null(ref) && !is.null(alt)) {
    paste(chr, pos, ref, alt, sep = ":")
  } else {
    paste0("X", seq_len(p))
  }

  if (!inherits(XtX, "Matrix")) { XtX <- Matrix::Matrix(XtX, sparse = TRUE) }
  if (nrow(XtX) != p || ncol(XtX) != p) {
    stop(paste0(
      "XtX dimensions (", nrow(XtX), " x ", ncol(XtX),
      ") do not match sumstats (", p, " variants)."
    ), call. = FALSE)
  }

  beta.external <- as.matrix(beta.external)
  if (any(is.na(beta.external))) { stop("Missing data (NA's) detected in beta.external.", call. = FALSE) }
  if (nrow(beta.external) != p) {
    stop(paste0("beta.external must have ", p, " rows (no intercept). Got ", nrow(beta.external), "."), call. = FALSE)
  }

  if (criteria %in% .bopt_val_criteria()) {
    if (is.null(X.val) || is.null(y.val)) {
      stop("X.val and y.val are required for validation-based criteria.", call. = FALSE)
    }
    X.val <- as.matrix(X.val)
    y.val <- as.numeric(as.vector(y.val))
    if (nrow(X.val) != length(y.val)) {
      stop("X.val and y.val do not have the same number of observations.", call. = FALSE)
    }
  } else {
    if (criteria %in% c("Cp", "GIC") && is.null(TN)) {
      stop("Training sample size TN must be provided for ", criteria, ".", call. = FALSE)
    }
    if (is.null(XtX.val)) { XtX.val <- XtX }
    if (is.null(sumstats.val)) { sumstats.val <- sumstats }
  }

  ## Redundant sources go first: the stacking solve inverts t(B) XtX B, which is
  ## singular the moment two columns carry the same information.
  dedup <- .dedup_and_report(beta.external, dedup.cor, intercept.row = FALSE)
  beta.external <- dedup$beta.external
  if (dedup$applied && !is.null(bounds) && length(bounds) == dedup$M.in) {
    ## Names are dropped so the survivors are renumbered eta_1 .. eta_M rather
    ## than keeping the gaps of the panel they came from, which would disagree
    ## with the column names of eta.grid.
    bounds <- unname(bounds[dedup$keep])
  }
  if (!is.null(init.grid) && !is.null(dim(init.grid)) &&
      ncol(init.grid) == dedup$M.in && dedup$applied) {
    init.grid <- init.grid[, dedup$keep, drop = FALSE]
  }

  # -- External reduction: computed once, reused every step --
  ext <- calcExtXtY(XtX, XtY, beta.external, multi.method, pca.var)
  XtY.external <- ext$XtY.external
  M <- ncol(XtY.external)

  fit.args <- list(...)
  fit.args$XtX <- XtX
  fit.args$XtY <- XtY
  fit.args$XtY.external <- XtY.external
  fit.args$family <- family
  fit.args$varnames <- varnames

  bounds <- .bopt_bounds(bounds, M, "external model(s)")
  init.grid <- .bopt_init_grid(init.grid, bounds, M)

  eval.eta <- function(eta) {
    args <- fit.args
    args$eta <- as.numeric(eta)
    fit <- do.call(BRIERs.eta, args)
    row <- .bopt_score_s(fit, eta, criteria, X.val, y.val, XtX.val, sumstats.val, TN, h2)
    list(fit = fit, row = row)
  }

  opt <- .bopt_engine(
    eval.eta, M, bounds, init.grid, init.points, n.iter,
    acq, kappa, acq.eps, kernel, verbose
  )

  out <- list(
    XtY              = XtY,
    XtY.external     = XtY.external,
    beta.external    = beta.external,
    family           = family,
    eta.list         = lapply(seq_len(M), function(k) sort(unique(opt$eta.grid[, k]))),
    eta.grid         = opt$eta.grid,
    res              = opt$res,
    null.dev         = 1,
    p                = nrow(XtX),
    M                = M,
    varnames         = varnames,
    external.dedup   = if (dedup$applied) dedup else NULL,
    external.pca     = ext$external.pca,
    criteria         = criteria,
    eta.min          = opt$eta.min,
    eta.min.index    = opt$eta.min.index,
    lambda.min       = opt$lambda.min,
    lambda.min.index = opt$lambda.min.index,
    eta.lambda       = opt$eta.lambda,
    bopt             = opt$bopt
  )
  class(out) <- c("BRIER.bopt", "BRIER.selection", "BRIER")
  out
}


#' Bayesian optimization of eta for BRIER.FULL
#'
#' Tunes the integration weights for \code{\link{BRIERfull}} with sequential
#' model-based optimization instead of an exhaustive grid. Each proposed eta
#' vector reweights the pooled external cohorts (target observations keep
#' weight 1, external cohort k receives weight \code{eta_k}) and triggers one
#' penalized fit over the full lambda path.
#'
#' Like \code{\link{BRIERfull.selection}}, this supports validation criteria
#' only: information criteria are not defined for the pooled fit.
#'
#' Note that \code{BRIERfull} requires at least one external cohort and cannot
#' fit a target-only model. Setting an eta bound to include 0 downweights an
#' external cohort to zero weight, which is not the same as a target-only fit;
#' use \code{BRIERi} with \code{eta = 0} for that baseline.
#'
#' @param X A numeric matrix of pooled predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param cohort An integer vector of length n: 0 for target samples, positive
#'   integers for external cohorts.
#' @param family A string: "gaussian", "binomial", or "poisson".
#' @param criteria Validation criterion, as in \code{\link{BRIERfull.selection}}.
#' @param X.val,y.val Held-out validation data. Required.
#' @param bounds A named list of length-2 numeric vectors giving the search box
#'   for each external cohort. Defaults to \code{c(0, 10)} per cohort.
#' @param init.grid Optional initial design, as in \code{\link{BRIERi.bopt}}.
#' @param init.points Number of random initial points when \code{init.grid} is
#'   not supplied.
#' @param n.iter Number of sequential optimization steps after initialization.
#' @param acq Acquisition function: "ucb", "ei", or "poi".
#' @param kappa Exploration parameter for "ucb".
#' @param acq.eps Exploration parameter for "ei" and "poi".
#' @param kernel Gaussian-process kernel specification.
#' @param verbose Logical. Print optimizer progress.
#' @param ... Further arguments passed to \code{\link{BRIERi.eta}}.
#'
#' @return An object of class \code{c("BRIER.bopt", "BRIER.selection", "BRIER")},
#'   shaped exactly like a \code{\link{BRIERfull.selection}} result with an
#'   added \code{bopt} element recording the search.
#'
#' @seealso \code{\link{BRIERfull}}, \code{\link{BRIERfull.selection}},
#'   \code{\link{BRIERi.bopt}}, \code{\link{BRIERs.bopt}}
#'
#' @examples
#' \dontrun{
#' sel <- BRIERfull.bopt(
#'   X, y, cohort,
#'   family      = "gaussian",
#'   criteria    = "gaussian.mspe",
#'   X.val       = X.val,
#'   y.val       = y.val,
#'   init.points = 10,
#'   n.iter      = 15
#' )
#' }
#'
#' @export
BRIERfull.bopt <- function(
  X, y, cohort, family = c("gaussian", "binomial", "poisson"),
  criteria = c(
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  X.val = NULL, y.val = NULL,
  bounds = NULL, init.grid = NULL, init.points = 10L, n.iter = 15L,
  acq = "ucb", kappa = 2.576, acq.eps = 0,
  kernel = list(type = "exponential", power = 2),
  verbose = TRUE,
  ...
) {

  family <- match.arg(family)
  criteria <- match.arg(criteria)

  # -- Validation of y (mirrors BRIERfull) --
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
  if (M == 0) {
    stop(
      "At least one external cohort (cohort > 0) is required. ",
      "Please use function BRIERi if you don't want to include external information.",
      call. = FALSE
    )
  }

  if (is.null(X.val) || is.null(y.val)) {
    stop("X.val and y.val are required: BRIER.FULL supports validation criteria only.", call. = FALSE)
  }
  X.val <- as.matrix(X.val)
  y.val <- as.numeric(as.vector(y.val))
  if (nrow(X.val) != length(y.val)) {
    stop("X.val and y.val do not have the same number of observations.", call. = FALSE)
  }

  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$family <- family
  fit.args$eta <- 0
  fit.args$y.external <- NULL

  penalty.factor <- if (!is.null(fit.args$penalty.factor)) fit.args$penalty.factor else rep(1, ncol(X))
  target_idx <- cohort == 0
  n_target <- sum(target_idx)
  null.dev <- calcNullDev(
    X[target_idx, , drop = FALSE], y[target_idx, , drop = FALSE],
    rep(1, n_target) / n_target, penalty.factor, family
  )

  bounds <- .bopt_bounds(bounds, M, "external cohort(s)")
  init.grid <- .bopt_init_grid(init.grid, bounds, M)

  eval.eta <- function(eta) {
    args <- fit.args
    w <- rep(1.0, length(cohort))
    for (k in seq_along(ext_cohorts)) {
      w[cohort == ext_cohorts[k]] <- eta[k]
    }
    args$weights <- w
    fit <- do.call(BRIERi.eta, args)
    fit$eta <- as.numeric(eta)
    row <- .bopt_score_i(fit, eta, criteria, X.val, y.val, NULL, NULL, y)
    list(fit = fit, row = row)
  }

  opt <- .bopt_engine(
    eval.eta, M, bounds, init.grid, init.points, n.iter,
    acq, kappa, acq.eps, kernel, verbose
  )

  out <- list(
    y                = y,
    y.external       = NULL,
    beta.external    = 0,
    family           = family,
    cohort           = cohort,
    eta.list         = lapply(seq_len(M), function(k) sort(unique(opt$eta.grid[, k]))),
    eta.grid         = opt$eta.grid,
    res              = opt$res,
    null.dev         = null.dev,
    n                = nrow(X),
    p                = ncol(X),
    M                = M,
    criteria         = criteria,
    eta.min          = opt$eta.min,
    eta.min.index    = opt$eta.min.index,
    lambda.min       = opt$lambda.min,
    lambda.min.index = opt$lambda.min.index,
    eta.lambda       = opt$eta.lambda,
    bopt             = opt$bopt
  )
  class(out) <- c("BRIER.bopt", "BRIER.selection", "BRIER")
  out
}
