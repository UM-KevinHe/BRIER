#' Restrict a BRIER fit or selection to some of its eta values
#'
#' Returns the same kind of object, restricted to the requested rows of the eta
#' grid, so \code{predict()}, \code{coef()} and the plots work on it unchanged.
#' The typical use is the no-transfer comparison: \code{subsetEta(sel, eta = 0)}
#' is the target-only model of the same fit, with its own best lambda.
#'
#' \itemize{
#'   \item On a tuned object (from \code{BRIERi.selection},
#'     \code{BRIERs.selection}, \code{BRIERfull.selection}, \code{BRIERi.cv} or a
#'     \code{*.bopt} function), the result is re-selected WITHIN the requested
#'     etas: \code{eta.min} is the best of them by the original criterion, and
#'     \code{lambda.min} is that eta's own best lambda. With a single eta, the
#'     result is pinned to it.
#'   \item On an untuned fit, the result is the fit restricted to those etas;
#'     run a selection function on it as usual.
#' }
#'
#' Only fitted etas can be returned. Values given in \code{eta} are matched to the
#' nearest fitted row and accepted when the relative difference
#' \code{|a - b| / max(|a|, |b|)} is within \code{tol} for every component; a value
#' with no fitted row that close is an error listing the fitted grid (refit with
#' it in \code{eta.list} to use it).
#'
#' @param object An object of class \code{"BRIER"} (a fit, or a tuned
#'   \code{"BRIER.selection"} / \code{"BRIER.cv"} object).
#' @param eta Eta values to keep: a numeric vector when there is one external
#'   model (or one stacked predictor), or a matrix with one column per external
#'   model, one row per eta combination.
#' @param which.eta Alternatively, integer row indices into
#'   \code{object$eta.grid}. Give \code{eta} or \code{which.eta}, not both.
#' @param tol Relative tolerance for matching \code{eta} to the fitted grid.
#'   Default \code{0.01}.
#'
#' @return An object of the same class as \code{object}, with \code{res},
#'   \code{eta.grid}, \code{eta.list} and (if tuned) \code{eta.lambda},
#'   \code{eta.min}, \code{eta.min.index}, \code{lambda.min},
#'   \code{lambda.min.index} restricted and re-indexed, and \code{eta.subset}
#'   recording the rows of the original grid that were kept.
#'
#' @seealso \code{\link{predict.BRIER}}, \code{\link{coef.BRIER}},
#'   \code{\link{BRIERi.selection}}, \code{\link{BRIERs.selection}}
#'
#' @examples
#' \dontrun{
#' fit <- BRIERi(X, y, family = "gaussian", beta.external = beta.external)
#' sel <- BRIERi.selection(fit, criteria = "gaussian.mspe", X.val = X.val, y.val = y.val)
#'
#' # target-only (no transfer) model of the same fit, with its own best lambda
#' target.only <- subsetEta(sel, eta = 0)
#' evalMetric(predict(target.only, X.test), y.test, "gaussian.rsq")
#'
#' # the best model with eta restricted to (0, 1]
#' small <- subsetEta(sel, eta = sel$eta.grid[sel$eta.grid[, 1] <= 1, 1])
#' }
#'
#' @export
subsetEta <- function(object, eta = NULL, which.eta = NULL, tol = 1e-2) {

  if (!inherits(object, "BRIER")) {
    stop("Object must be of class 'BRIER', got '", class(object)[1], "'.", call. = FALSE)
  }
  if (!is.null(eta) && !is.null(which.eta)) {
    stop("Give 'eta' or 'which.eta', not both.", call. = FALSE)
  }
  if (is.null(eta) && is.null(which.eta)) {
    stop("Give the eta values to keep, as 'eta' or 'which.eta'.", call. = FALSE)
  }

  grid <- as.matrix(object$eta.grid)
  if (!is.null(eta)) {
    eta <- if (is.matrix(eta) || is.data.frame(eta)) as.matrix(eta) else
      matrix(as.numeric(eta), ncol = if (ncol(grid) == 1) 1 else length(eta))
    idx <- .eta_rows(object, eta, tol = tol)
  } else {
    idx <- as.integer(which.eta)
    if (anyNA(idx) || any(idx < 1) || any(idx > nrow(grid))) {
      stop("which.eta must be between 1 and ", nrow(grid), ".", call. = FALSE)
    }
  }
  idx <- sort(unique(idx))

  out <- object
  out$res <- object$res[idx]
  out$eta.grid <- grid[idx, , drop = FALSE]
  out$eta.list <- if (ncol(out$eta.grid) == 1) {
    as.numeric(out$eta.grid[, 1])
  } else {
    lapply(seq_len(ncol(out$eta.grid)), function(j) unique(out$eta.grid[, j]))
  }
  out$eta.subset <- if (!is.null(object$eta.subset)) object$eta.subset[idx] else idx

  # A tuned object: keep each kept eta's own best lambda, and re-select the best eta
  # among them by the criterion the selection used (eta.lambda$measure.min, lower is
  # better in every selection function).
  el <- object$eta.lambda
  if (!is.null(el)) {
    if (is.null(el$eta.index)) {
      stop("object$eta.lambda has no eta.index column; cannot subset it.", call. = FALSE)
    }
    rows <- match(idx, el$eta.index)
    if (anyNA(rows)) {
      stop(
        "eta row(s) ", paste(idx[is.na(rows)], collapse = ", "),
        " have no entry in eta.lambda.", call. = FALSE
      )
    }
    el <- el[rows, , drop = FALSE]
    el$eta.index <- seq_along(idx)
    rownames(el) <- NULL
    best <- which.min(el$measure.min)
    out$eta.lambda <- el
    out$eta.min <- out$eta.grid[best, ]
    out$eta.min.index <- best
    out$lambda.min.index <- el$lambda.min.index[best]
    out$lambda.min <- el$lambda.min[best]
  }

  out
}
