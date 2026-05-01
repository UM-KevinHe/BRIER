# -- Internal helper: set up lambda.max and the lambda sequence --

setupLambda <- function(
  X, y, n, p, wt, penalty.factor, alpha,
  family, lambda.min, log.lambda, nlambda
) {

  ind <- which(penalty.factor != 0)
  if (length(ind) != p) {
    fit <- glm(y ~ X[, -ind], weights = wt, family = family)
  } else {
    fit <- glm(y ~ 1, weights = wt, family = family)
  }

  ## Determine lambda.max
  if (family == "gaussian") {
    r <- fit$residuals
    r <- r * fit$weights
  } else {
    w <- fit$weights
    # if (max(w) < 1e-4) stop("Unpenalized portion of model is already saturated; exiting...", call. = FALSE)
    r <- residuals(fit, "working") * w
  }
  lambda.max <- maxprod(X, r, n, p, penalty.factor)
  lambda.max <- lambda.max / alpha

  if (log.lambda) { # lambda sequence on log-scale
    if (lambda.min == 0) {
      lambda <- c(exp(seq(log(lambda.max), log(.001 * lambda.max), length = nlambda - 1)), 0)
    } else {
      lambda <- exp(seq(log(lambda.max), log(lambda.min * lambda.max), length = nlambda))
    }
  } else { # lambda sequence on linear-scale
    if (lambda.min == 0) {
      lambda <- c(seq(lambda.max, 0.001 * lambda.max, length = nlambda - 1), 0)
    } else {
      lambda <- seq(lambda.max, lambda.min * lambda.max, length = nlambda)
    }
  }

  lambda
}

