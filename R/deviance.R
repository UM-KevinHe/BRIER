#' Forward link function
#'
#' Maps the mean (mu) to the linear predictor (eta) under a specified
#' GLM family.
#'
#' @param mu A numeric vector or matrix of mean values.
#' @param family A string: "gaussian", "binomial", or "poisson".
#'
#' @return A numeric vector or matrix of linear predictors on the link scale.
#'
#' @seealso \code{\link{ginv_link}}
#'
#' @export
g_link <- function(mu, family) {
  if (family == "gaussian") {
    val <- mu
  } else if (family == "binomial") {
    val <- log(mu / (1 - mu))
  } else if (family == "poisson") {
    val <- log(mu)
  } else {
    stop("family must be either 'gaussian', 'binomial', or 'poisson'.", call. = FALSE)
  }
  val
}


#' Inverse link function
#'
#' Maps the linear predictor (eta) to the mean (mu) under a specified GLM
#' family.
#'
#' @param val A numeric vector or matrix of linear predictors.
#' @param family A string: "gaussian", "binomial", or "poisson".
#'
#' @return A numeric vector or matrix of means on the response scale.
#'
#' @seealso \code{\link{g_link}}
#'
#' @export
ginv_link <- function(val, family) {
  if (family == "gaussian") {
    mu <- val
  } else if (family == "binomial") {
    mu <- 1 / (1 + exp(-val))
  } else if (family == "poisson") {
    mu <- exp(val)
  } else {
    stop("family must be either 'gaussian', 'binomial', or 'poisson'.", call. = FALSE)
  }
  mu
}


#' Null deviance
#'
#' Compute the average deviance of the intercept-only (or unpenalized-only)
#' model for a specified GLM family. Variables with \code{penalty.factor == 0}
#' are included as unpenalized covariates in the null model.
#'
#' @param X A numeric matrix of predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param wt A numeric vector of observation weights (length n).
#' @param penalty.factor A numeric vector of per-variable penalty multipliers
#'   (length p). Variables with multiplier 0 are unpenalized and included in
#'   the null model.
#' @param family A string: "gaussian", "binomial", or "poisson".
#'
#' @return A numeric scalar: the average null deviance (\code{sum(dev) / sum(wt)}).
#'
#' @seealso \code{\link{calcDev}}
#'
#' @export
calcNullDev <- function(X, y, wt, penalty.factor, family) {

  idx0 <- penalty.factor == 0
  if (any(idx0)) {
    df <- data.frame(y = y, as.data.frame(X[, idx0, drop = FALSE]))
    fit <- glm(y ~ ., data = df, weights = wt, family = family)
  } else {
    df <- data.frame(y = y)
    fit <- glm(y ~ 1, data = df, weights = wt, family = family)
  }
  yh <- predict(fit, type = "response")

  ## Compute and return the average deviance
  dev <- calcDev(y, yh, wt, family)
  sum(dev) / sum(wt)
}


#' Element-wise deviance
#'
#' Compute element-wise (weighted) deviance contributions for Gaussian,
#' binomial, or Poisson families. The result has the same dimensions as
#' \code{yhat}, allowing column-wise summation when \code{yhat} is a matrix
#' (e.g. one column per lambda value).
#'
#' @param y A numeric response vector or single-column matrix of length n.
#' @param yhat A numeric matrix of fitted means (n x L), or a vector of
#'   length n.
#' @param wt A numeric vector of observation weights (length n).
#'   Defaults to uniform weights.
#' @param family A string: "gaussian", "binomial", or "poisson".
#'
#' @return A numeric matrix (n x L) of element-wise deviance contributions.
#'   Sum over rows to get total deviance per lambda.
#'
#' @seealso \code{\link{calcNullDev}}
#'
#' @export
calcDev <- function(y, yhat, wt = rep(1, nrow(y)) / nrow(y), family) {

  yhat <- as.matrix(yhat)
  y <- as.matrix(y)
  if (ncol(y) > 1) { stop("y should contain only one column.", call. = FALSE) }
  Y <- matrix(as.numeric(y), nrow = nrow(yhat), ncol = ncol(yhat))

  wt <- as.numeric(wt)
  if (length(wt) != nrow(yhat)) { stop("wt must have length nrow(yhat).", call. = FALSE) }
  WT <- matrix(wt, nrow = nrow(yhat), ncol = ncol(yhat))

  dev <- matrix(NA, nrow = nrow(yhat), ncol = ncol(yhat))
  if (family == "gaussian") {
    dev <- WT * (Y - yhat)^2
  } else if (family == "binomial") {
    yhat[yhat < 1e-8] <- 1e-8
    yhat[yhat > 1 - 1e-8] <- 1 - 1e-8
    dev <- -2 * WT * (log(yhat) * Y + log(1 - yhat) * (1 - Y))
  } else if (family == "poisson") {
    yly <- Y * log(Y)
    yly[Y == 0] <- 0
    dev <- 2 * WT * (yly - Y + yhat - Y * log(yhat))
  }

  dev
}