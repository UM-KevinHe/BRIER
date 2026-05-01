#' Model selection for BRIER.I fits
#'
#' Selects the optimal eta combination and lambda value from a fitted
#' \code{BRIERi} object using either information criteria computed on the
#' training data, or held-out validation data.
#'
#' @param object An object of class \code{"BRIER"} from \code{\link{BRIERi}}.
#' @param criteria A string specifying the selection criterion. One of:
#'   \itemize{
#'     \item Information criteria (training data):
#'       \code{"gcv"}, \code{"AIC"}, \code{"BIC"}, \code{"Cp"}.
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
#' @param n Optional integer. Sample size to use for information criteria
#'   (AIC, BIC, Cp, gcv). Defaults to \code{length(object$y)}.
#' @param var.y Optional numeric. Response variance to use for Gaussian Cp.
#'   Defaults to \code{var(object$y)}.
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
#' @seealso \code{\link{BRIERi}}, \code{\link{BRIERi.cv}},
#'   \code{\link{BRIERs.selection}}, \code{\link{BRIERfull.selection}},
#'   \code{\link{coef.BRIER}}, \code{\link{predict.BRIER}}
#'
#' @examples
#' \dontrun{
#' fit <- BRIERi(
#'   X, y, family = "gaussian",
#'   eta.list = list(c(0, 0.5, 1, 2)),
#'   beta.external = beta.external,
#'   penalty = "LASSO"
#' )
#'
#' # information criterion
#' sel <- BRIERi.selection(fit, criteria = "BIC")
#'
#' # validation
#' sel.val <- BRIERi.selection(
#'   fit, criteria = "gaussian.mspe",
#'   X.val = X.val, y.val = y.val
#' )
#'
#' coef(sel)
#' }
#'
#' @export
BRIERi.selection <- function(
  object, 
  criteria = c(
    "gcv", "AIC", "BIC", "Cp",
    "gaussian.mspe", "gaussian.rsq",
    "binomial.dev", "binomial.mcfrsq", "binomial.tjursq",
    "binomial.auc", "poisson.dev"
  ),
  X.val = NULL, y.val = NULL, 
  n = NULL, var.y = NULL
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

    eta.lambda <- do.call(rbind, lapply(seq_len(n.fits), function(i) {
      validation(i, object, X.val, y.val, criteria)
    }))

  # -- Information criteria (use training data) --
  } else {
    n <- if (!is.null(n)) n else length(object$y)
    var.y  <- if (!is.null(var.y)) var.y else var(object$y)
    eta.lambda <- do.call(rbind, lapply(seq_len(n.fits), function(i) {
      ic_selection(i, object, n, criteria, dispersion = NULL, var_y = var.y)
    }))
  }

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


# -- Internal helper: validation (held-out predictive measures) --

validation <- function(i, object, X.val, y.val, criteria) {

  fit <- object$res[[i]]
  eta <- fit$eta
  lambda <- fit$lambda

  # sanity check: eta stored in fit should match eta.grid
  if (!all(abs(eta - object$eta.grid[i, ]) < 1e-10)) {
    stop(paste0(
      "Mismatch: fit$eta = (", paste(round(eta, 4), collapse = ", "),
      ") but eta.grid[", i, ", ] = (",
      paste(round(object$eta.grid[i, ], 4), collapse = ", "), ")."
    ), call. = FALSE)
  }

  mu.pred <- predict.BRIER.fit(fit, X = X.val, type = "response")
  tol <- 1e-8

  if (grepl("binomial", criteria)) {
    mu.pred <- pmin(pmax(mu.pred, tol), 1 - tol)
    null_prob <- mean(y.val)
    null_prob <- pmin(pmax(null_prob, tol), 1 - tol)
  } else if (grepl("poisson", criteria)) {
    mu.pred <- pmax(mu.pred, tol)
  }

  measure <- if (criteria == "gaussian.mspe") {
    colMeans((mu.pred - as.vector(y.val))^2)

  } else if (criteria == "gaussian.rsq") {
    rsq <- as.vector(-cor(mu.pred, as.vector(y.val))^2)
    rsq[is.na(rsq)] <- 0
    rsq

  } else if (criteria == "binomial.dev") {
    apply(mu.pred, 2, function(pred) {
      -2 * sum(y.val * log(pred) + (1 - y.val) * log(1 - pred))
    })

  } else if (criteria == "binomial.mcfrsq") {
    null_loglik <- sum(log(null_prob) * y.val + log(1 - null_prob) * (1 - y.val))
    loglik <- colSums(log(mu.pred) * y.val + log(1 - mu.pred) * (1 - y.val))
    -(1 - loglik / null_loglik)

  } else if (criteria == "binomial.tjursq") {
    m1 <- mu.pred[which(y.val == 1), , drop = FALSE]
    m2 <- mu.pred[which(y.val == 0), , drop = FALSE]
    -(colMeans(m1) - colMeans(m2))

  } else if (criteria == "binomial.auc") {
    -apply(mu.pred, 2, function(x) {
      if (all(is.na(x))) { return(NA) }
      pROC::roc(response = y.val, predictor = x, levels = c(0, 1), direction = "<")$auc
    })

  } else if (criteria == "poisson.dev") {
    apply(mu.pred, 2, function(mu) {
      term <- ifelse(y.val == 0, 0, y.val * log(y.val / mu))
      2 * sum(term - (y.val - mu))
    })
  }

  measure[!is.finite(measure)] <- NA
  if (all(is.na(measure))) { stop("All validation measures are NA for fit ", i, ".", call. = FALSE) }

  min.idx <- which.min(replace(measure, is.na(measure), Inf))

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


# -- Internal helper: information criteria on training data --

ic_selection <- function(i, object, n, criteria, dispersion = NULL, var_y = NULL) {

  fit <- object$res[[i]]
  lambda <- fit$lambda
  family <- fit$family
  dev <- fit$deviance
  eta <- fit$eta

  # sanity check: eta stored in fit should match eta.grid
  if (!all(abs(eta - object$eta.grid[i, ]) < 1e-10)) {
    stop(paste0(
      "Mismatch: fit$eta = (", paste(round(eta, 4), collapse = ", "),
      ") but eta.grid[", i, ", ] = (",
      paste(round(object$eta.grid[i, ], 4), collapse = ", "), ")."
    ), call. = FALSE)
  }

  if (is.null(dev)) { stop("fit$deviance is required for model selection.", call. = FALSE) }

  # effective df adjusted for external information
  df <- fit$k / (1 + sum(eta))

  # dispersion
  if (is.null(dispersion)) {
    if (family %in% c("binomial", "poisson")) {
      dispersion <- rep(1, length(lambda))
    } else if (family == "gaussian") {
      dispersion <- dev / pmax(n - fit$k, 1)
    }
  }
  if (length(dispersion) == 1L) { dispersion <- rep(dispersion, length(lambda)) }
  if (length(dispersion) != length(lambda)) {
    stop("Length of dispersion does not match length of lambda.", call. = FALSE)
  }

  measure <- if (criteria == "AIC") {
    dev + 2 * df

  } else if (criteria == "BIC") {
    dev + log(n) * df

  } else if (criteria == "Cp") {
    if (family == "gaussian") {
      if (is.null(var_y)) {
        stop("var_y must be provided for Gaussian Cp.", call. = FALSE)
      }
      dev / var_y + 2 * df
    } else {
      dev / dispersion + 2 * df - n
    }

  } else if (criteria == "gcv") {
    denom <- (1 - df / n)^2
    denom[denom <= 0] <- NA_real_
    (dev / n) / denom

  } else {
    stop("criteria must be one of 'AIC', 'BIC', 'Cp', or 'gcv'.", call. = FALSE)
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
