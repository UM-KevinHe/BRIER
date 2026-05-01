#' Model selection for BRIER.FULL fits
#'
#' Selects the optimal eta combination and lambda value from a fitted
#' \code{BRIERfull} object using held-out validation data. Information criteria
#' (AIC, BIC, etc.) are not supported for \code{BRIERfull} because external
#' observations are integrated via weights, which complicates the effective
#' degrees-of-freedom adjustment.
#'
#' @param object An object of class \code{"BRIER"} from \code{\link{BRIERfull}}.
#' @param criteria A string specifying the validation criterion. One of:
#'   \itemize{
#'     \item \code{"gaussian.mspe"}: mean squared prediction error.
#'     \item \code{"gaussian.rsq"}: negative R-squared (minimised).
#'     \item \code{"binomial.dev"}: binomial deviance.
#'     \item \code{"binomial.mcfrsq"}: McFadden's R-squared (negated, minimised).
#'     \item \code{"binomial.tjursq"}: Tjur's R-squared (negated, minimised).
#'     \item \code{"binomial.auc"}: AUC (negated, minimised).
#'     \item \code{"poisson.dev"}: Poisson deviance.
#'   }
#' @param X.val A numeric matrix of validation predictors (n_val x p).
#' @param y.val A numeric response vector for validation (length n_val).
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
#' @seealso \code{\link{BRIERfull}}, \code{\link{BRIERi.selection}},
#'   \code{\link{coef.BRIER}}, \code{\link{predict.BRIER}}
#'
#' @examples
#' \dontrun{
#' fit <- BRIERfull(
#'   X, y, cohort, family = "gaussian",
#'   eta.list = list(c(0, 0.5, 1, 2)),
#'   penalty = "LASSO"
#' )
#' sel <- BRIERfull.selection(
#'   fit, criteria = "gaussian.mspe",
#'   X.val = X.val, y.val = y.val
#' )
#' coef(sel)
#' }
#'
#' @export
BRIERfull.selection <- function(
  object, 
  criteria = c(
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  X.val = NULL, y.val = NULL
) {

  if (!inherits(object, "BRIER")) {
    stop("Object must be of class 'BRIER', got '", class(object)[1], "'.", call. = FALSE)
  }

  criteria <- match.arg(criteria)
  n.fits <- length(object$res)

  # -- Validation-based criteria (need held-out data) --
  if (is.null(X.val) || is.null(y.val)) {
    stop("X.val and y.val are required for validation-based criteria.", call. = FALSE)
  }
  X.val <- as.matrix(X.val)
  y.val <- as.numeric(as.vector(y.val))
  if (nrow(X.val) != length(y.val)) {
    stop("X.val and y.val do not have the same number of observations.", call. = FALSE)
  }

  eta.lambda <- do.call(rbind, lapply(seq_len(n.fits), function(i) {
    validation(i, object, X.val, y.val, criteria)
  }))


  # -- Find optimal eta and lambda --
  eta.min.index <- which.min(eta.lambda$measure.min)
  eta.min <- object$eta.grid[eta.min.index, ]
  lambda.min.index <- eta.lambda$lambda.min.index[eta.min.index]
  lambda.min <- eta.lambda$lambda.min[eta.min.index]

  if (ncol(object$eta.grid) == 1) {
    cat("Best eta: ", round(eta.min, 3), "\n")
  } else {
    cat("Best eta: (", paste(round(eta.min, 3), collapse = ", "), ")\n")
  }
  cat("Best lambda:", round(lambda.min, 3), "\n")

  # -- Build output: inherit BRIER object and add selection results --
  out <- object
  out$criteria <- criteria
  out$eta.min <- eta.min
  out$eta.min.index <- eta.min.index
  out$lambda.min <- lambda.min
  out$lambda.min.index <- lambda.min.index
  out$eta.lambda <- eta.lambda
  class(out) <- c("BRIER.selection", "BRIER")
  out
}

