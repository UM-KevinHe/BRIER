#' Compute a predictive performance metric
#'
#' Evaluate the predictive performance of a vector of predictions on
#' observed outcomes under a specified criterion. Used internally by
#' \code{\link{plot.eta}} and \code{\link{plot.box}}.
#'
#' @param pred A numeric vector of predicted values on the response scale.
#'   For binomial criteria these should be probabilities in \eqn{[0, 1]};
#'   for Poisson, non-negative rates.
#' @param y A numeric vector of observed outcomes (same length as \code{pred}).
#' @param criteria A string specifying the criterion. One of:
#'   \itemize{
#'     \item \code{"gaussian.mspe"}: mean squared prediction error.
#'     \item \code{"gaussian.rsq"}: coefficient of determination.
#'     \item \code{"binomial.dev"}: binomial deviance.
#'     \item \code{"binomial.mcfrsq"}: McFadden's R-squared.
#'     \item \code{"binomial.tjursq"}: Tjur's R-squared (mean prediction
#'       gap between positives and negatives).
#'     \item \code{"binomial.auc"}: area under the ROC curve.
#'     \item \code{"poisson.dev"}: Poisson deviance.
#'   }
#'
#' @return A numeric scalar: the value of the criterion.
#'
#' @seealso \code{\link{plot.eta}}, \code{\link{BRIERi.selection}}
#'
#' @examples
#' \dontrun{
#'  set.seed(1)
#'  y <- rnorm(100)
#'  pred <- y + rnorm(100, sd = 0.5)
#'  evalMetric(pred, y, "gaussian.mspe")
#'  evalMetric(pred, y, "gaussian.rsq")
#' }
#' @export
evalMetric <- function(pred, y, criteria) {

  criteria <- match.arg(criteria, c(
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ))

  pred <- as.numeric(pred)
  y    <- as.numeric(y)
  if (length(pred) != length(y)) {
    stop("pred and y must have the same length.", call. = FALSE)
  }
  tol <- 1e-8

  if (criteria == "gaussian.mspe") {
    return(mean((pred - y)^2))
  }
  if (criteria == "gaussian.rsq") {
    r <- suppressWarnings(stats::cor(pred, y))
    return(if (is.na(r)) 0 else r^2)
  }
  if (criteria == "binomial.dev") {
    pred <- pmin(pmax(pred, tol), 1 - tol)
    return(-2 * mean(y * log(pred) + (1 - y) * log(1 - pred)))
  }
  if (criteria == "binomial.mcfrsq") {
    pred <- pmin(pmax(pred, tol), 1 - tol)
    p0   <- pmin(pmax(mean(y), tol), 1 - tol)
    null.ll <- sum(log(p0) * y + log(1 - p0) * (1 - y))
    ll      <- sum(log(pred) * y + log(1 - pred) * (1 - y))
    return(1 - ll / null.ll)
  }
  if (criteria == "binomial.tjursq") {
    return(mean(pred[y == 1]) - mean(pred[y == 0]))
  }
  if (criteria == "binomial.auc") {
    if (!requireNamespace("pROC", quietly = TRUE)) {
      stop(
        "Package 'pROC' is required for AUC. ",
        "Install with: install.packages('pROC').", call. = FALSE
      )
    }
    return(as.numeric(
      pROC::roc(response = y, predictor = pred, levels = c(0, 1),
                direction = "<", quiet = TRUE)$auc
    ))
  }
  if (criteria == "poisson.dev") {
    pred <- pmax(pred, tol)
    term <- ifelse(y == 0, 0, y * log(y / pred))
    return(2 * mean(term - (y - pred)))
  }
  NA_real_
}
