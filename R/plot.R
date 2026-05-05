#' Plot and evaluate predictive performance as a function of integration weight eta
#'
#' Visualises predictive performance of a fitted \code{BRIER} model as a
#' function of the integration weight \eqn{\eta}. Operates in two modes:
#'
#' \describe{
#'   \item{Selection-criterion mode}{When \code{X} is not supplied, the
#'     function plots the per-eta optimum of the selection criterion stored
#'     in \code{object$eta.lambda} (e.g. CV error, BIC, GIC). This requires
#'     that \code{object} has been processed by \code{\link{BRIERi.cv}},
#'     \code{\link{BRIERi.selection}}, \code{\link{BRIERs.selection}}, or
#'     \code{\link{BRIERfull.selection}}.}
#'   \item{Validation mode}{When \code{X} and \code{covar.data} are supplied,
#'     the function evaluates the criterion on the held-out testing set,
#'     optionally with bootstrap CIs.}
#' }
#'
#' For models fitted with a single external model (\code{ncol(eta.grid) == 1}),
#' a \code{ggplot} object is returned. The black triangle marks the optimal
#' \eqn{\eta} (from \code{object$eta.min} when available). For models with
#' multiple external models (\code{ncol(eta.grid) > 1}), no plot is produced;
#' the function returns the summary \code{data.frame} only.
#'
#' For genetic applications such as polygenic risk score (PRS) evaluation,
#' supply \code{covar.data} as a \code{data.frame} or matrix with
#' \code{pheno.name} for the trait column and \code{adjust.covar} for
#' confounder columns. Confounder adjustment is performed by regressing the
#' phenotype on the confounders and using the residuals; when
#' \code{bootstrap = TRUE}, the regression is refit within each bootstrap
#' subsample so that uncertainty in the confounder adjustment is propagated
#' into the empirical CI band.
#'
#' Confounder adjustment and y-standardization are only valid for Gaussian
#' criteria. When \code{criteria} is a binomial or Poisson criterion, both
#' are automatically skipped with a warning. The design matrix \code{X} is
#' always standardizable when \code{standardize.data = TRUE}, regardless of
#' family.
#'
#' @param object An object of class \code{"BRIER"}. For selection-criterion
#'   mode, \code{object} must have an \code{eta.lambda} component (attached
#'   by \code{\link{BRIERi.cv}} or any \code{*.selection()} function).
#' @param X Optional numeric matrix of testing-set predictors (n_test x p).
#'   If supplied, the function evaluates \code{criteria} on the testing set;
#'   if missing, the function plots the selection criterion stored in
#'   \code{object$eta.lambda}.
#' @param covar.data A vector, matrix, or \code{data.frame} of testing-set
#'   outcomes and (optionally) confounders. Required when \code{X} is
#'   supplied; ignored otherwise. Vectors are treated as the phenotype and
#'   given column name \code{pheno.name}. Matrices and data.frames must have
#'   column names; the trait column is selected via \code{pheno.name} and
#'   confounders via \code{adjust.covar}.
#' @param criteria A string specifying the evaluation criterion when \code{X}
#'   is supplied. One of: \code{"gaussian.mspe"}, \code{"gaussian.rsq"},
#'   \code{"binomial.dev"}, \code{"binomial.mcfrsq"},
#'   \code{"binomial.tjursq"}, \code{"binomial.auc"}, \code{"poisson.dev"}.
#'   Ignored in selection-criterion mode (the criterion in
#'   \code{object$criteria} is used instead). See \code{\link{evalMetric}} for
#'   definitions.
#' @param pheno.name Character. Column name of the trait of interest in
#'   \code{covar.data}. Defaults to \code{"y"}.
#' @param adjust.covar Optional character vector of confounder column names
#'   in \code{covar.data} to regress out before evaluation. Ignored with a
#'   warning for binomial and Poisson criteria.
#' @param standardize.data Logical. If TRUE, the design matrix \code{X} is
#'   standardized to columns with mean 0 and variance 1. The phenotype is
#'   also standardized when \code{criteria} corresponds to a Gaussian
#'   outcome; for binomial and Poisson criteria, y-standardization is skipped
#'   with a warning.
#' @param bootstrap Logical. If TRUE, evaluates the criterion over
#'   \code{bootstrap.n} bootstrap subsamples of the testing set. Confounder
#'   adjustment, if requested, is refit within each subsample. Ignored in
#'   selection-criterion mode.
#' @param bootstrap.size Integer. Sample size for each bootstrap subsample.
#'   Defaults to \code{floor(0.8 * nrow(X))}.
#' @param bootstrap.n Integer. Number of bootstrap replicates. Defaults to 100.
#' @param seed Optional integer seed for reproducible bootstrap sampling.
#' @param ... Unused; present for S3 method compatibility.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{plot}{A \code{ggplot} object (only when \code{ncol(eta.grid) == 1}),
#'     otherwise \code{NULL}.}
#'   \item{summary.df}{A \code{data.frame} with one row per eta combination,
#'     containing eta values, the criterion name, the metric, and the
#'     empirical 95\% CI bounds (\code{NA} when no bootstrap or in
#'     selection-criterion mode).}
#'   \item{bootstrap.mat}{The bootstrap replicate matrix (validation mode with
#'     bootstrap), or \code{NULL} otherwise.}
#' }
#'
#' @seealso \code{\link{plot.box}}, \code{\link{BRIERi}}, \code{\link{BRIERi.selection}},
#'   \code{\link{BRIERi.cv}}, \code{\link{predict.BRIER}},
#'   \code{\link{evalMetric}}
#'
#' @examples
#' \dontrun{
#' fit <- BRIERi(
#'   X, y, family = "gaussian",
#'   eta.list = list(c(0, 0.25, 0.5, 1, 2, 4)),
#'   beta.external = beta.external
#' )
#' sel <- BRIERi.selection(fit, criteria = "BIC")
#'
#' # Selection-criterion mode (no X needed)
#' out <- plot.eta(sel)
#' out$plot
#'
#' # Validation mode with bootstrap
#' out <- plot.eta(
#'   sel, X = X_testing, covar.data = pheno_testing,
#'   criteria = "gaussian.rsq",
#'   pheno.name = "trait", adjust.covar = c("age", "sex"),
#'   standardize.data = TRUE,
#'   bootstrap = TRUE, bootstrap.n = 100, seed = 42
#' )
#' out$plot
#' head(out$summary.df)
#' }
#'
#' @export plot.eta
plot.eta <- function(
  object, X, covar.data,
  criteria = c(
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  pheno.name = "y", adjust.covar = NULL, standardize.data = FALSE,
  bootstrap = FALSE, bootstrap.size = NULL, bootstrap.n = 100L,
  seed = NULL,
  ...
) {

  if (!inherits(object, "BRIER")) {
    stop("object must be of class 'BRIER'.", call. = FALSE)
  }

  M <- ncol(object$eta.grid)
  do_plot <- (M == 1)

  if (do_plot && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plotting. ",
      "Install with: install.packages('ggplot2').", call. = FALSE
    )
  }

  # ===================================================================
  # Selection-criterion mode (no X supplied)
  # ===================================================================
  if (missing(X)) {
    if (is.null(object$eta.lambda)) {
      stop(
        "object$eta.lambda not found. Run BRIERi.cv(), BRIERi.selection(), ",
        "BRIERs.selection(), or BRIERfull.selection() first, ",
        "or supply X and covar.data to plot validation performance.",
        call. = FALSE)
    }

    sel.criteria <- if (!is.null(object$criteria)) object$criteria else "criterion"

    eta.df <- as.data.frame(object$eta.grid)
    colnames(eta.df) <- paste0("eta_", seq_len(M))
    summary.df <- cbind(
      eta.df,
      data.frame(
        criteria  = sel.criteria,
        metric    = object$eta.lambda$measure.min,
        metric.lo = NA_real_,
        metric.hi = NA_real_,
        stringsAsFactors = FALSE
      )
    )

    if (!do_plot) {
      return(list(plot = NULL, summary.df = summary.df, bootstrap.mat = NULL))
    }

    eta.values  <- as.numeric(object$eta.grid[, 1])
    metric.vals <- object$eta.lambda$measure.min

    best.eta <- if (!is.null(object$eta.min)) {
      as.numeric(object$eta.min)
    } else {
      eta.values[which.min(metric.vals)]
    }
    best.value <- metric.vals[which.min(abs(eta.values - best.eta))]

    df <- data.frame(eta = eta.values, metric = metric.vals)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$eta, y = .data$metric)) +
      ggplot2::geom_line(colour = "#2ca25f", linewidth = 1) +
      ggplot2::geom_point(colour = "#2ca25f", size = 2) +
      ggplot2::geom_point(
        data = data.frame(eta = best.eta, metric = best.value),
        ggplot2::aes(x = .data$eta, y = .data$metric),
        shape = 17, size = 4, colour = "black"
      ) +
      ggplot2::labs(
        x = expression(eta), y = metric_label(sel.criteria),
        title = NULL
      ) +
      ggplot2::theme_classic(base_size = 13, base_family = "serif")

    return(list(plot = p, summary.df = summary.df, bootstrap.mat = NULL))
  }

  # ===================================================================
  # Validation mode (X and covar.data supplied)
  # ===================================================================
  if (missing(covar.data)) {
    stop("covar.data must be supplied when X is supplied.", call. = FALSE)
  }
  criteria <- match.arg(criteria)
  X <- as.matrix(X)

  # -- Determine family from criterion --
  fam <- sub("\\..*$", "", criteria)

  # -- Family-aware standardization and confounder adjustment --
  standardize.X <- isTRUE(standardize.data)
  standardize.y <- isTRUE(standardize.data) && fam == "gaussian"
  if (isTRUE(standardize.data) && fam != "gaussian") {
    warning(
      "standardize.data = TRUE: standardizing X but not y for ", fam,
      " criteria (y standardization is invalid for non-Gaussian outcomes).",
      call. = FALSE
    )
  }
  if (!is.null(adjust.covar) && length(adjust.covar) > 0 && fam != "gaussian") {
    warning(
      "adjust.covar is ignored for ", fam, " criteria. ",
      "Confounder adjustment via lm() residuals is only valid for Gaussian outcomes. ",
      "If you need adjusted predictions, fit BRIER with the confounders included as covariates.",
      call. = FALSE
    )
    adjust.covar <- NULL
  }

  # -- Coerce covar.data to a matrix with column names --
  if (is.vector(covar.data) && !is.list(covar.data)) {
    covar.data <- matrix(as.numeric(covar.data), ncol = 1)
    colnames(covar.data) <- pheno.name
  } else if (is.data.frame(covar.data)) {
    covar.data <- as.matrix(covar.data)
  } else if (!is.matrix(covar.data)) {
    stop("covar.data must be a vector, matrix, or data.frame.", call. = FALSE)
  }
  if (is.null(colnames(covar.data))) {
    stop("covar.data must have column names.", call. = FALSE)
  }

  # -- Validate pheno.name and adjust.covar columns --
  if (!pheno.name %in% colnames(covar.data)) {
    stop("pheno.name '", pheno.name, "' not found in covar.data.", call. = FALSE)
  }
  if (!is.null(adjust.covar)) {
    missing.cols <- setdiff(adjust.covar, colnames(covar.data))
    if (length(missing.cols) > 0) {
      stop(
        "Columns not found in covar.data: ",
        paste(missing.cols, collapse = ", "), ".", call. = FALSE
      )
    }
  }

  # -- Validate dimensions --
  n.test <- nrow(covar.data)
  if (nrow(X) != n.test) {
    stop(
      "X and covar.data must have the same number of testing observations.",
      call. = FALSE
    )
  }

  n.eta <- nrow(object$eta.grid)

  # ----- Helper: evaluate metric for a given subset of rows -----
  # Standardizes X[idx, ] *within* the subset (no leakage), predicts at every
  # eta combination on the standardized subset, then evaluates the criterion.
  eval_on_subset <- function(idx) {
    X.b <- X[idx, , drop = FALSE]
    if (standardize.X) X.b <- standardize_X(X.b)[[1]]
    y.b <- resolve_y_subset(
      covar.data, idx, pheno.name, adjust.covar, standardize = standardize.y
    )
    vapply(seq_len(n.eta), function(i) {
      eta_row <- matrix(object$eta.grid[i, ], nrow = 1)
      pr <- drop(predict.BRIER(object, X = X.b, eta = eta_row, type = "response"))
      evalMetric(pr, y.b, criteria)
    }, numeric(1))
  }

  # -- Compute criterion (with optional bootstrap) --
  if (isTRUE(bootstrap)) {
    if (is.null(bootstrap.size)) bootstrap.size <- floor(0.8 * n.test)
    bootstrap.size <- as.integer(bootstrap.size)
    bootstrap.n    <- as.integer(bootstrap.n)
    if (bootstrap.size < 1L || bootstrap.size > n.test) {
      stop("bootstrap.size must be in [1, n_test].", call. = FALSE)
    }

    if (!is.null(seed)) { set.seed(seed) }

    boot.mat <- matrix(NA_real_, nrow = bootstrap.n, ncol = n.eta)
    for (b in seq_len(bootstrap.n)) {
      idx <- sample.int(n.test, bootstrap.size, replace = TRUE)
      boot.mat[b, ] <- eval_on_subset(idx)
    }

    metric.vals <- colMeans(boot.mat, na.rm = TRUE)
    metric.lo   <- apply(boot.mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
    metric.hi   <- apply(boot.mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  } else {
    metric.vals <- eval_on_subset(seq_len(n.test))
    metric.lo <- metric.hi <- NULL
    boot.mat <- NULL
  }

  # -- Build summary data.frame --
  eta.df <- as.data.frame(object$eta.grid)
  colnames(eta.df) <- paste0("eta_", seq_len(M))
  summary.df <- cbind(
    eta.df,
    data.frame(
      criteria  = criteria,
      metric    = metric.vals,
      metric.lo = if (!is.null(metric.lo)) metric.lo else NA_real_,
      metric.hi = if (!is.null(metric.hi)) metric.hi else NA_real_,
      stringsAsFactors = FALSE
    )
  )

  # -- Multiple external models: return summary only, no plot --
  if (!do_plot) {
    return(list(plot = NULL, summary.df = summary.df, bootstrap.mat = boot.mat))
  }

  # -- Single external model: build ggplot --
  eta.values <- as.numeric(object$eta.grid[, 1])
  best.eta <- if (!is.null(object$eta.min)) {
    as.numeric(object$eta.min)
  } else {
    is_lower_better <- criteria %in% c("gaussian.mspe", "binomial.dev", "poisson.dev")
    if (is_lower_better) eta.values[which.min(metric.vals)]
    else                 eta.values[which.max(metric.vals)]
  }
  best.value <- metric.vals[which.min(abs(eta.values - best.eta))]

  df <- data.frame(eta = eta.values, metric = metric.vals)
  if (!is.null(metric.lo)) { df$lo <- metric.lo; df$hi <- metric.hi }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$eta, y = .data$metric))
  if (!is.null(metric.lo)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
      alpha = 0.15, fill = "#2ca25f"
    )
  }
  p <- p +
    ggplot2::geom_line(colour = "#2ca25f", linewidth = 1) +
    ggplot2::geom_point(colour = "#2ca25f", size = 2) +
    ggplot2::geom_point(
      data = data.frame(eta = best.eta, metric = best.value),
      ggplot2::aes(x = .data$eta, y = .data$metric),
      shape = 17, size = 4, colour = "black"
    ) +
    ggplot2::labs(x = expression(eta), y = metric_label(criteria), title = NULL) +
    ggplot2::theme_classic(base_size = 13, base_family = "serif")

  list(
    plot          = p,
    summary.df    = summary.df,
    bootstrap.mat = boot.mat
  )
}


#' Boxplot of predictive performance across bootstrap resamples
#'
#' Visualises the distribution of predictive performance across bootstrap
#' resamples for the baseline (target-only) model, each external model alone,
#' and the optimal \code{BRIER} model. This complements
#' \code{\link{plot.eta}} by showing the variability introduced by sampling
#' randomness on the testing set.
#'
#' Three categories of models are evaluated:
#' \describe{
#'   \item{Baseline}{The target-only fit at \eqn{\eta = 0} (no external
#'     information). Requires \eqn{\eta = 0} to be present in the eta grid.}
#'   \item{External k}{The k-th external model alone, evaluated by
#'     applying \code{object$beta.external[, k]} directly to \code{X}.
#'     One box per external model.}
#'   \item{BRIER}{The optimal \code{BRIER} model at \code{object$eta.min}
#'     (and corresponding \code{lambda.min}). Requires that the object has
#'     been processed by a selection function.}
#' }
#'
#' For genetic applications such as polygenic risk score (PRS) evaluation,
#' supply \code{covar.data} as a \code{data.frame} or matrix with
#' \code{pheno.name} for the trait column and \code{adjust.covar} for
#' confounder columns. Confounder adjustment is performed by regressing the
#' phenotype on the confounders within each bootstrap subsample so that
#' uncertainty in the confounder adjustment is propagated into the
#' empirical distribution.
#'
#' Confounder adjustment and y-standardization are only valid for Gaussian
#' criteria. When \code{criteria} is a binomial or Poisson criterion, both
#' are automatically skipped with a warning. The design matrix \code{X} is
#' always standardizable when \code{standardize.data = TRUE}, regardless of
#' family.
#'
#' @param object An object of class \code{"BRIER"} with selection results
#'   attached (\code{object$eta.min} and \code{object$beta.external}
#'   required for the BRIER and external boxes).
#' @param X A numeric matrix of testing-set predictors (n_test x p).
#' @param covar.data A vector, matrix, or \code{data.frame} of testing-set
#'   outcomes and (optionally) confounders. Vectors are treated as the
#'   phenotype and given column name \code{pheno.name}.
#' @param criteria A string specifying the evaluation criterion. One of:
#'   \code{"gaussian.mspe"}, \code{"gaussian.rsq"},
#'   \code{"binomial.dev"}, \code{"binomial.mcfrsq"},
#'   \code{"binomial.tjursq"}, \code{"binomial.auc"}, \code{"poisson.dev"}.
#'   See \code{\link{evalMetric}} for definitions.
#' @param pheno.name Character. Column name of the trait of interest in
#'   \code{covar.data}. Defaults to \code{"y"}.
#' @param adjust.covar Optional character vector of confounder column names.
#'   Ignored with a warning for binomial and Poisson criteria.
#' @param standardize.data Logical. If TRUE, the design matrix \code{X} is
#'   standardized to columns with mean 0 and variance 1. The phenotype is
#'   also standardized for Gaussian criteria.
#' @param bootstrap.size Integer. Sample size for each bootstrap subsample.
#'   Defaults to \code{floor(0.8 * nrow(X))}.
#' @param bootstrap.n Integer. Number of bootstrap replicates. Defaults to 100.
#' @param seed Optional integer seed for reproducible bootstrap sampling.
#' @param ... Unused; present for S3 method compatibility.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{plot}{A \code{ggplot} boxplot.}
#'   \item{summary.df}{A long \code{data.frame} of bootstrap replicates with
#'     columns \code{model} (factor: "Baseline", "External 1", ..., "BRIER")
#'     and \code{metric}.}
#'   \item{bootstrap.mat}{The wide-form bootstrap matrix
#'     (\code{bootstrap.n x n.model}).}
#' }
#'
#' @seealso \code{\link{plot.eta}}, \code{\link{evalMetric}},
#'   \code{\link{BRIERi.selection}}, \code{\link{predict.BRIER}}
#'
#' @examples
#' \dontrun{
#' fit <- BRIERi(
#'   X, y, family = "gaussian",
#'   eta.list = list(c(0, 0.25, 0.5, 1, 2, 4)),
#'   beta.external = beta.external
#' )
#' sel <- BRIERi.selection(fit, criteria = "BIC")
#'
#' out <- plot.box(
#'   sel, X = X_testing, covar.data = pheno_testing,
#'   criteria = "gaussian.rsq",
#'   pheno.name = "trait", adjust.covar = c("age", "sex"),
#'   standardize.data = TRUE,
#'   bootstrap.size = 0.8 * nrow(X_testing),
#'   bootstrap.n = 100, seed = 42
#' )
#' out$plot
#' }
#'
#' @export plot.box
plot.box <- function(
  object, X, covar.data,
  criteria = c(
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  pheno.name = "y", adjust.covar = NULL, standardize.data = FALSE,
  bootstrap.size = NULL, bootstrap.n = 100L, seed = NULL,
  ...
) {

  if (!inherits(object, "BRIER")) {
    stop("object must be of class 'BRIER'.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plotting. ",
      "Install with: install.packages('ggplot2').", call. = FALSE
    )
  }
  if (is.null(object$eta.min) || is.null(object$lambda.min)) {
    stop(
      "object must have selection results attached. ",
      "Run BRIERi.selection(), BRIERs.selection(), BRIERfull.selection(), ",
      "or BRIERi.cv() first.", call. = FALSE
    )
  }
  if (is.null(object$beta.external)) {
    stop(
      "object$beta.external is required for plot.box(). ",
      "Available for BRIERi (with beta.external argument) and BRIERs models.",
      call. = FALSE
    )
  }
  criteria <- match.arg(criteria)

  X <- as.matrix(X)
  fam <- sub("\\..*$", "", criteria)

  # -- Family-aware standardization and confounder adjustment --
  standardize.X <- isTRUE(standardize.data)
  standardize.y <- isTRUE(standardize.data) && fam == "gaussian"
  if (isTRUE(standardize.data) && fam != "gaussian") {
    warning(
      "standardize.data = TRUE: standardizing X but not y for ", fam,
      " criteria (y standardization is invalid for non-Gaussian outcomes).",
      call. = FALSE
    )
  }
  if (!is.null(adjust.covar) && length(adjust.covar) > 0 && fam != "gaussian") {
    warning(
      "adjust.covar is ignored for ", fam, " criteria. ",
      "Confounder adjustment via lm() residuals is only valid for Gaussian outcomes.",
      call. = FALSE
    )
    adjust.covar <- NULL
  }

  # -- Coerce covar.data --
  if (is.vector(covar.data) && !is.list(covar.data)) {
    covar.data <- matrix(as.numeric(covar.data), ncol = 1)
    colnames(covar.data) <- pheno.name
  } else if (is.data.frame(covar.data)) {
    covar.data <- as.matrix(covar.data)
  } else if (!is.matrix(covar.data)) {
    stop("covar.data must be a vector, matrix, or data.frame.", call. = FALSE)
  }
  if (is.null(colnames(covar.data))) {
    stop("covar.data must have column names.", call. = FALSE)
  }
  if (!pheno.name %in% colnames(covar.data)) {
    stop("pheno.name '", pheno.name, "' not found in covar.data.", call. = FALSE)
  }
  if (!is.null(adjust.covar)) {
    missing.cols <- setdiff(adjust.covar, colnames(covar.data))
    if (length(missing.cols) > 0) {
      stop(
        "Columns not found in covar.data: ",
        paste(missing.cols, collapse = ", "), ".", call. = FALSE
      )
    }
  }

  n.test <- nrow(covar.data)
  if (nrow(X) != n.test) {
    stop(
      "X and covar.data must have the same number of testing observations.",
      call. = FALSE
    )
  }
  if (ncol(X) != object$p) {
    stop(
      "ncol(X) (= ", ncol(X), ") does not match object$p (= ", object$p, ").",
      call. = FALSE
    )
  }

  # -- Locate the eta = 0 row in the grid for the baseline --
  eta.grid <- object$eta.grid
  zero.row <- which(rowSums(abs(eta.grid)) < 1e-10)
  if (length(zero.row) == 0) {
    stop(
      "Cannot find eta = (0, ...) in object$eta.grid for the baseline. ",
      "Include 0 in your eta.list when fitting BRIER.", call. = FALSE
    )
  }
  zero.row <- zero.row[1]
  base.lam.idx <- object$eta.lambda$lambda.min.index[
    object$eta.lambda$eta.index == zero.row
  ]
  if (length(base.lam.idx) == 0L || is.na(base.lam.idx[1])) {
    stop(
      "No eta.lambda entry for the eta = 0 row (eta.index = ", zero.row,
      "). Did selection complete successfully?", call. = FALSE
    )
  }
  base.lam.idx <- base.lam.idx[1]

  # -- Determine M and model labels (using full-X dimensions, not per-bootstrap) --
  fam.obj <- object$family
  be      <- object$beta.external
  if (nrow(be) == object$p + 1) {
    has.intercept <- TRUE
  } else if (nrow(be) == object$p) {
    has.intercept <- FALSE
  } else {
    stop("beta.external dimensions do not match X.", call. = FALSE)
  }
  M <- ncol(be)
  model.labels <- c(
    "Baseline",
    if (M == 1) "External" else paste0("External", seq_len(M)),
    "BRIER"
  )
  n.model <- length(model.labels)

  # ----- Helper: compute predictions for a given subset of rows -----
  # Standardizes X[idx, ] *within* the subset (no leakage), then computes
  # predictions for Baseline, all External(s), and BRIER on the standardized
  # subset. Returns an (n.idx x n.model) matrix.
  predict_on_subset <- function(idx) {
    X.b <- X[idx, , drop = FALSE]
    if (standardize.X) X.b <- standardize_X(X.b)[[1]]

    # Baseline (eta = 0)
    base.pred <- drop(predict.BRIER(
      object, X = X.b,
      which.eta = zero.row,
      which.lambda = base.lam.idx,
      type = "response"
    ))

    # External-only predictions
    if (has.intercept) {
      ext.lin <- X.b %*% be[-1, , drop = FALSE]
      ext.lin <- sweep(ext.lin, 2, be[1, ], "+")
    } else {
      ext.lin <- X.b %*% be
    }
    ext.pred <- ginv_link(ext.lin, fam.obj)   # n.idx x M

    # Optimal BRIER prediction (uses object$eta.min / object$lambda.min)
    brier.pred <- drop(predict.BRIER(object, X = X.b, type = "response"))

    out <- cbind(base.pred, ext.pred, brier.pred)
    colnames(out) <- model.labels
    out
  }

  # -- Bootstrap loop --
  if (is.null(bootstrap.size)) bootstrap.size <- floor(0.8 * n.test)
  bootstrap.size <- as.integer(bootstrap.size)
  bootstrap.n    <- as.integer(bootstrap.n)
  if (bootstrap.size < 1L || bootstrap.size > n.test) {
    stop("bootstrap.size must be in [1, n_test].", call. = FALSE)
  }

  if (!is.null(seed)) { set.seed(seed) }

  boot.mat <- matrix(NA_real_, nrow = bootstrap.n, ncol = n.model)
  colnames(boot.mat) <- model.labels

  for (b in seq_len(bootstrap.n)) {
    idx <- sample.int(n.test, bootstrap.size, replace = TRUE)
    pred.b <- predict_on_subset(idx)
    y.b <- resolve_y_subset(
      covar.data, idx, pheno.name, adjust.covar, standardize = standardize.y
    )
    for (j in seq_len(n.model)) {
      boot.mat[b, j] <- evalMetric(pred.b[, j], y.b, criteria)
    }
  }

  # -- Build long-form data.frame --
  summary.df <- data.frame(
    model  = factor(rep(model.labels, each = bootstrap.n), levels = model.labels),
    metric = as.numeric(boot.mat),
    stringsAsFactors = FALSE
  )

  # -- Build ggplot --
  fill.colors <- c(
    "Baseline" = "#bdbdbd",
    "BRIER"    = "#2ca25f"
  )
  ext.fill <- "#74a9cf"
  for (lab in model.labels) {
    if (!lab %in% names(fill.colors)) fill.colors[lab] <- ext.fill
  }

  p <- ggplot2::ggplot(summary.df, ggplot2::aes(
        x = .data$model, y = .data$metric, fill = .data$model
      )) +
    ggplot2::geom_boxplot(
      width = 0.6, outlier.size = 0.8, colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = fill.colors, guide = "none") +
    ggplot2::labs(
      x = NULL,
      y = metric_label(criteria),
      title = NULL
    ) +
    ggplot2::theme_classic(base_size = 13, base_family = "serif") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )

  list(
    plot          = p,
    summary.df    = summary.df,
    bootstrap.mat = boot.mat
  )
}


#' Plot variable importance via bootstrap stability selection
#'
#' Provides a stability-based assessment of variable importance by repeatedly
#' re-running selection on the fitted \code{BRIER} model using bootstrap
#' resamples of the validation set as the held-out data. For each replication,
#' \code{\link{BRIERi.selection}} is re-run on the bootstrap subsample to
#' pick a new (\eqn{\eta^*}, \eqn{\lambda^*}); the procedure then records
#' which predictors enter the model (have nonzero coefficient) at that
#' choice via \code{\link{coef.BRIER}}. Selection frequencies are summarised
#' across all replications.
#'
#' This works for all three model classes (\code{BRIERi}, \code{BRIERs},
#' \code{BRIERfull}) because they share the \code{"BRIER"} class and
#' \code{\link{BRIERi.selection}} dispatches uniformly. It re-uses the
#' existing fit and only re-runs lambda/eta selection on each bootstrap
#' subsample - it does not refit BRIER. Standardization of \code{X}, when
#' requested, is applied per bootstrap subsample to avoid information
#' leakage from out-of-subsample observations. The intercept (when present
#' in the fitted coefficient matrix) is automatically excluded from the
#' importance summary, since it is not a predictor of interest.
#'
#' The resulting plot displays the top \code{n.top} predictors ranked by
#' selection frequency. High-frequency predictors indicate stable,
#' reproducible signals; low-frequency predictors may reflect noise or
#' dataset-specific fluctuations.
#'
#' @param object An object of class \code{"BRIER"} from \code{\link{BRIERi}},
#'   \code{\link{BRIERs}}, or \code{\link{BRIERfull}}.
#' @param X A numeric matrix of validation-set predictors (n_val x p).
#' @param covar.data A vector, matrix, or \code{data.frame} of validation-set
#'   outcomes and (optionally) confounders. Vectors are treated as the
#'   phenotype and given column name \code{pheno.name}.
#' @param criteria A string specifying the validation criterion passed to
#'   \code{\link{BRIERi.selection}} on each bootstrap subsample. One of:
#'   \code{"gaussian.mspe"}, \code{"gaussian.rsq"}, \code{"binomial.dev"},
#'   \code{"binomial.mcfrsq"}, \code{"binomial.tjursq"},
#'   \code{"binomial.auc"}, \code{"poisson.dev"}.
#' @param pheno.name Character. Column name of the trait of interest in
#'   \code{covar.data}. Defaults to \code{"y"}.
#' @param adjust.covar Optional character vector of confounder column names.
#'   Ignored with a warning for binomial and Poisson criteria.
#' @param standardize.data Logical. If TRUE, the design matrix \code{X} is
#'   standardized per bootstrap subsample; the phenotype is also
#'   standardized for Gaussian criteria.
#' @param n.top Integer. Number of top predictors to display, ranked by
#'   selection frequency. Defaults to 20.
#' @param replications Integer. Number of bootstrap replications. Defaults
#'   to 100.
#' @param bootstrap.size Integer. Sample size for each bootstrap subsample.
#'   Defaults to \code{floor(0.8 * nrow(X))}.
#' @param seed Optional integer seed for reproducible bootstrap sampling.
#' @param ... Unused; present for S3 method compatibility.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{plot}{A \code{ggplot} bar plot of the top \code{n.top} predictors
#'     by selection frequency.}
#'   \item{summary.df}{A \code{data.frame} with one row per predictor
#'     (intercept excluded), containing the variable name, selection count,
#'     and selection frequency, sorted in descending order of frequency.}
#'   \item{selection.mat}{A logical matrix (\code{replications x p}) where
#'     each row indicates which predictors were selected (nonzero
#'     coefficient) in that replication.}
#' }
#'
#' @seealso \code{\link{plot.eta}}, \code{\link{plot.box}},
#'   \code{\link{BRIERi.selection}}, \code{\link{coef.BRIER}}
#'
#' @examples
#' \dontrun{
#' fit <- BRIERi(
#'   X, y, family = "gaussian",
#'   eta.list = list(c(0, 0.25, 0.5, 1, 2, 4)),
#'   beta.external = beta.external
#' )
#' sel <- BRIERi.selection(fit, criteria = "BIC")
#'
#' out <- plot.importance(
#'   sel, X = X_validation, covar.data = pheno_validation,
#'   criteria = "gaussian.mspe",
#'   pheno.name = "trait", adjust.covar = c("age", "sex"),
#'   n.top = 10, replications = 100, seed = 1
#' )
#' out$plot
#' head(out$summary.df, 10)
#' }
#'
#' @export plot.importance
plot.importance <- function(
  object, X, covar.data,
  criteria = c(
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  pheno.name = "y", adjust.covar = NULL, standardize.data = FALSE,
  n.top = 20L, replications = 100L,
  bootstrap.size = NULL, seed = NULL,
  ...
) {

  if (!inherits(object, "BRIER")) {
    stop("object must be of class 'BRIER'.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plotting. ",
      "Install with: install.packages('ggplot2').", call. = FALSE
    )
  }
  criteria <- match.arg(criteria)

  X <- as.matrix(X)
  fam <- sub("\\..*$", "", criteria)

  # -- Family-aware standardization and confounder adjustment --
  standardize.X <- isTRUE(standardize.data)
  standardize.y <- isTRUE(standardize.data) && fam == "gaussian"
  if (isTRUE(standardize.data) && fam != "gaussian") {
    warning(
      "standardize.data = TRUE: standardizing X but not y for ", fam,
      " criteria (y standardization is invalid for non-Gaussian outcomes).",
      call. = FALSE
    )
  }
  if (!is.null(adjust.covar) && length(adjust.covar) > 0 && fam != "gaussian") {
    warning(
      "adjust.covar is ignored for ", fam, " criteria. ",
      "Confounder adjustment via lm() residuals is only valid for Gaussian outcomes.",
      call. = FALSE
    )
    adjust.covar <- NULL
  }

  # -- Coerce covar.data --
  if (is.vector(covar.data) && !is.list(covar.data)) {
    covar.data <- matrix(as.numeric(covar.data), ncol = 1)
    colnames(covar.data) <- pheno.name
  } else if (is.data.frame(covar.data)) {
    covar.data <- as.matrix(covar.data)
  } else if (!is.matrix(covar.data)) {
    stop("covar.data must be a vector, matrix, or data.frame.", call. = FALSE)
  }
  if (is.null(colnames(covar.data))) {
    stop("covar.data must have column names.", call. = FALSE)
  }
  if (!pheno.name %in% colnames(covar.data)) {
    stop("pheno.name '", pheno.name, "' not found in covar.data.", call. = FALSE)
  }
  if (!is.null(adjust.covar)) {
    missing.cols <- setdiff(adjust.covar, colnames(covar.data))
    if (length(missing.cols) > 0) {
      stop(
        "Columns not found in covar.data: ",
        paste(missing.cols, collapse = ", "), ".", call. = FALSE
      )
    }
  }

  n.val <- nrow(covar.data)
  if (nrow(X) != n.val) {
    stop(
      "X and covar.data must have the same number of validation observations.",
      call. = FALSE
    )
  }

  # -- Detect intercept by matrix dimensions, then resolve predictor names --
  if (length(object$res) == 0L) {
    stop("object$res is empty.", call. = FALSE)
  }
  beta.template <- object$res[[1]]$beta
  n.beta        <- nrow(beta.template)
  p             <- object$p

  if (ncol(X) != p) {
    stop(
      "ncol(X) (= ", ncol(X), ") does not match the number of predictors in ",
      "the fitted model (= ", p, "). ",
      "Make sure X is the validation design matrix with the same predictors ",
      "as the training data.", call. = FALSE
    )
  }

  if (n.beta == p + 1L) {
    has.intercept <- TRUE
  } else if (n.beta == p) {
    has.intercept <- FALSE
  } else {
    stop(
      "Number of rows in beta (", n.beta, ") does not match the number of ",
      "predictors stored in object$p (", p, ") or p + 1.", call. = FALSE
    )
  }

  all.names <- rownames(beta.template)
  if (has.intercept) {
    var.names <- if (!is.null(all.names)) all.names[-1L] else paste0("V", seq_len(p))
  } else {
    var.names <- if (!is.null(all.names)) all.names else paste0("V", seq_len(p))
  }

  # -- Bootstrap loop --
  if (is.null(bootstrap.size)) bootstrap.size <- floor(0.8 * n.val)
  bootstrap.size <- as.integer(bootstrap.size)
  replications   <- as.integer(replications)
  if (bootstrap.size < 1L || bootstrap.size > n.val) {
    stop("bootstrap.size must be in [1, n_val].", call. = FALSE)
  }

  if (!is.null(seed)) { set.seed(seed) }

  selection.mat <- matrix(FALSE, nrow = replications, ncol = p)
  colnames(selection.mat) <- var.names

  for (b in seq_len(replications)) {
    idx <- sample.int(n.val, bootstrap.size, replace = TRUE)

    # Standardize subset (no leakage)
    X.b <- X[idx, , drop = FALSE]
    if (standardize.X) X.b <- standardize_X(X.b)[[1]]

    # Confounder-adjusted, optionally standardized y
    y.b <- resolve_y_subset(
      covar.data, idx, pheno.name, adjust.covar, standardize = standardize.y
    )

    # Re-run selection on the bootstrap subsample
    sel.b <- BRIERi.selection(
      object, criteria = criteria,
      X.val = X.b, y.val = y.b
    )

    # Extract coefficients at the bootstrap-selected (eta*, lambda*)
    beta.b <- coef.BRIER(sel.b)

    # Strip intercept (deterministic, based on whether the original beta has it)
    if (has.intercept) beta.b <- beta.b[-1L]

    selection.mat[b, ] <- as.numeric(beta.b) != 0
  }

  # -- Summarise selection frequencies --
  freq  <- colMeans(selection.mat)
  count <- colSums(selection.mat)
  summary.df <- data.frame(
    variable  = var.names,
    count     = count,
    frequency = freq,
    stringsAsFactors = FALSE
  )
  summary.df <- summary.df[order(-summary.df$frequency, summary.df$variable), ]
  rownames(summary.df) <- NULL

  # -- Build ggplot of top n.top --
  n.top <- min(as.integer(n.top), nrow(summary.df))
  top.df <- summary.df[seq_len(n.top), ]
  top.df$variable <- factor(top.df$variable, levels = rev(top.df$variable))

  p.plot <- ggplot2::ggplot(
      top.df,
      ggplot2::aes(x = .data$frequency, y = .data$variable)
    ) +
    ggplot2::geom_col(fill = "#2ca25f", width = 0.7) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x = "Selection frequency",
      y = NULL,
      title = NULL
    ) +
    ggplot2::theme_classic(base_size = 13, base_family = "serif")

  list(
    plot          = p.plot,
    summary.df    = summary.df,
    selection.mat = selection.mat
  )
}


# -- Internal helper: resolve outcome on a subset --

resolve_y_subset <- function(covar.data, idx, pheno.name, adjust.covar, standardize) {

  pheno <- as.numeric(covar.data[idx, pheno.name])

  if (!is.null(adjust.covar) && length(adjust.covar) > 0) {
    covar.df <- as.data.frame(covar.data[idx, adjust.covar, drop = FALSE])
    lm.fit <- stats::lm(
      pheno ~ .,
      data = data.frame(pheno = pheno, covar.df)
    )
    pheno <- as.numeric(stats::residuals(lm.fit))
  }

  if (isTRUE(standardize)) {
    pheno <- as.numeric(standardize_X(matrix(pheno, ncol = 1))[[1]])
  }
  pheno
}



# -- Internal helper: axis label --

metric_label <- function(criteria) {
  switch(criteria,
    gaussian.mspe        = "MSPE",
    gaussian.rsq         = expression("Pearson's " * R^2),
    binomial.dev         = "Binomial Deviance",
    binomial.mcfrsq      = expression("McFadden's " * R^2),
    binomial.tjursq      = expression("Tjur's " * R^2),
    binomial.auc         = "AUC",
    poisson.dev          = "Poisson Deviance",
    cve                  = "Cross-Validated Error",
    AIC                  = "AIC",
    BIC                  = "BIC",
    Cp                   = expression(C[p]),
    GIC                  = "GIC",
    gcv                  = "Generalized CV",
    pseu.val             = "Pseudo-Validation Score",
    criteria
  )
}