#' Model selection for BRIER.S fits
#'
#' Selects the optimal eta combination and lambda value from a fitted
#' \code{BRIERs} object using either summary-statistic-based information
#' criteria or held-out validation data.
#'
#' @param object An object of class \code{"BRIER"} from \code{\link{BRIERs}}.
#' @param criteria A string specifying the selection criterion. One of:
#'   \itemize{
#'     \item Summary-stat criteria:
#'       \code{"Cp"}, \code{"GIC"} (require \code{TN} and optionally \code{h2});
#'       \code{"pseu.val"} (requires \code{XtX} and \code{sumstats}).
#'     \item Validation criteria (held-out data):
#'       \code{"gaussian.mspe"}, \code{"gaussian.rsq"},
#'       \code{"binomial.dev"}, \code{"binomial.mcfrsq"},
#'       \code{"binomial.tjursq"}, \code{"binomial.auc"},
#'       \code{"poisson.dev"}.
#'   }
#' @param X.val A numeric matrix of validation predictors (n_val x p).
#'   Required for validation criteria.
#' @param y.val A numeric response vector for validation (length n_val).
#'   Required for validation criteria.
#' @param XtX A p x p sparse LD matrix. Required for \code{"pseu.val"}.
#' @param sumstats A data.frame of GWAS summary statistics (with \code{corr}
#'   column). Required for \code{"pseu.val"}.
#' @param TN Integer. Training sample size. Required for \code{"Cp"} and
#'   \code{"GIC"}.
#' @param h2 Optional numeric. Heritability estimate for adjusting Cp/GIC.
#'   Defaults to 0 with a warning if not provided.
#'
#' @return An object of class \code{c("BRIER.selection", "BRIER")} which extends
#'   the input \code{object} with the following added elements:
#' \describe{
#'   \item{criteria}{The selection criterion used.}
#'   \item{eta.min}{The eta combination minimising the criterion.}
#'   \item{eta.min.index}{The corresponding row index in \code{eta.grid}.}
#'   \item{lambda.min}{The lambda value minimising the criterion at \code{eta.min}.}
#'   \item{lambda.min.index}{The corresponding lambda index.}
#'   \item{eta.lambda}{A data.frame summarising the optimal lambda for each eta.}
#' }
#'
#' @seealso \code{\link{BRIERs}}, \code{\link{BRIERi.selection}},
#'   \code{\link{BRIERfull.selection}}, \code{\link{coef.BRIER}},
#'   \code{\link{predict.BRIER}}
#'
#' @examples
#' \dontrun{
#' fit <- BRIERs(
#'   sumstats, XtX, family = "gaussian",
#'   eta.list = list(c(0, 0.5, 1, 2)),
#'   beta.external = beta.external
#' )
#'
#' # summary-stat selection
#' sel <- BRIERs.selection(fit, criteria = "GIC", TN = 50000, h2 = 0.3)
#'
#' # pseudo-validation
#' sel.pv <- BRIERs.selection(
#'   fit, criteria = "pseu.val",
#'   XtX = XtX, sumstats = sumstats
#' )
#'
#' # validation
#' sel.val <- BRIERs.selection(
#'   fit, criteria = "gaussian.mspe",
#'   X.val = X.val, y.val = y.val
#' )
#'
#' coef(sel)
#' }
#'
#' @export
BRIERs.selection <- function(
  object,
  criteria = c(
    "Cp", "GIC", "pseu.val",
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  X.val = NULL, y.val = NULL,
  XtX = NULL, sumstats = NULL, TN = NULL, h2 = NULL
) {

  if (!inherits(object, "BRIER")) {
    stop("Object must be of class 'BRIER', got '", class(object)[1], "'.", call. = FALSE)
  }

  criteria <- match.arg(criteria)
  n.fits <- length(object$res)

  # -- Validation-based criteria (need held-out data) --
  if (criteria %in% c(
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  )) {
    if (is.null(X.val) || is.null(y.val)) {
      stop("X.val and y.val are required for validation-based criteria.", call. = FALSE)
    }
    X.val <- as.matrix(X.val)
    y.val <- as.numeric(as.vector(y.val))
    if (nrow(X.val) != length(y.val)) {
      stop("X.val and y.val do not have the same number of observations.", call. = FALSE)
    }

    # reuse the same validation() function as BRIERi.selection
    eta.lambda <- do.call(rbind, lapply(seq_len(n.fits), function(i) {
      validation(i, object, X.val, y.val, criteria)
    }))

  # -- Summary-based criteria --
  } else {
    if (criteria %in% c("Cp", "GIC") && is.null(TN)) {
      stop("Training sample size TN must be provided for ", criteria, ".", call. = FALSE)
    }
    if (criteria == "pseu.val" && (is.null(XtX) || is.null(sumstats))) {
      stop("XtX and sumstats must be provided for pseudo-validation.", call. = FALSE)
    }

    eta.lambda <- do.call(rbind, lapply(seq_len(n.fits), function(i) {
      ic_selection_S(i, object, XtX, sumstats, TN, h2, criteria)
    }))
  }

  # -- Find optimal eta and lambda --
  eta.min.index    <- which.min(eta.lambda$measure.min)
  eta.min          <- object$eta.grid[eta.min.index, ]
  lambda.min.index <- eta.lambda$lambda.min.index[eta.min.index]
  lambda.min       <- eta.lambda$lambda.min[eta.min.index]

  if (ncol(object$eta.grid) == 1) {
    cat("Best eta:", round(eta.min, 3), "\n")
  } else {
    cat("Best eta: (", paste(round(eta.min, 3), collapse = ", "), ")\n")
  }
  cat("Best lambda:", round(lambda.min, 3), "\n")

  # -- Build output: inherit BRIER object and add selection results --
  out <- object
  out$criteria         <- criteria
  out$eta.min          <- eta.min
  out$eta.min.index    <- eta.min.index
  out$lambda.min       <- lambda.min
  out$lambda.min.index <- lambda.min.index
  out$eta.lambda       <- eta.lambda
  class(out) <- c("BRIER.selection", "BRIER")
  out
}


# -- Internal helper: summary-stat information criteria --

ic_selection_S <- function(i, object, XtX, sumstats, TN, h2, criteria) {

  fit    <- object$res[[i]]
  lambda <- fit$lambda
  eta    <- fit$eta
  dev    <- fit$deviance

  # sanity check
  if (!all(abs(eta - object$eta.grid[i, ]) < 1e-10)) {
    stop(paste0(
      "Mismatch: fit$eta = (", paste(round(eta, 4), collapse = ", "),
      ") but eta.grid[", i, ", ] = (",
      paste(round(object$eta.grid[i, ], 4), collapse = ", "), ")."
    ), call. = FALSE)
  }

  if (is.null(dev)) { stop("fit$deviance is required for model selection.", call. = FALSE) }

  df <- fit$k / (1 + sum(eta))

  measure <- if (criteria == "Cp") {
    if (is.null(h2)) {
      h2 <- 0
      warning("h2 not provided; defaulting to 0 for Cp.", call. = FALSE)
    }
    dev + (1 - h2) * 2 * df / TN

  } else if (criteria == "GIC") {
    if (is.null(h2)) {
      h2 <- 0
      warning("h2 not provided; defaulting to 0 for GIC.", call. = FALSE)
    }
    dev + (1 - h2) * log(TN) * df / TN

  } else if (criteria == "pseu.val") {
    corr <- as.numeric(as.vector(sumstats$corr))
    -pseudo_validation(fit$beta, XtX, corr)

  } else {
    stop("criteria must be one of 'Cp', 'GIC', or 'pseu.val'.", call. = FALSE)
  }

  min.idx <- which.min(measure)

  # build output with eta values
  eta_vals <- as.data.frame(t(eta))
  colnames(eta_vals) <- paste0("eta_", seq_along(eta))
  cbind(
    data.frame(eta.index = i),
    eta_vals,
    data.frame(
      criteria         = criteria,
      measure.min      = measure[min.idx],
      lambda.min.index = min.idx,
      lambda.min       = lambda[min.idx]
    )
  )
}