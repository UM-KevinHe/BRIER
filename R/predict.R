#' Coefficients from a single BRIER fit
#'
#' Extract coefficients from a fitted \code{BRIER.eta} object at one or more
#' lambda values. Lambda values not in the fitted path are obtained by linear
#' interpolation between adjacent fitted lambdas.
#'
#' @param object An object of class \code{"BRIER.eta"}.
#' @param lambda Optional numeric vector of lambda values at which to extract
#'   coefficients. Must lie within the fitted path. If \code{NULL} (default),
#'   coefficients are returned at the indices given by \code{which}.
#' @param which Integer vector specifying which lambda indices to return.
#'   Defaults to all fitted lambdas.
#' @param drop Logical. If TRUE, drop singleton dimensions from the output.
#' @param ... Unused; present for S3 method compatibility.
#'
#' @return A numeric matrix or vector of coefficients.
#'
#' @seealso \code{\link{predict.BRIER.eta}}, \code{\link{coef.BRIER}}
#'
#' @export
coef.BRIER.eta <- function(object, lambda = NULL, which = seq_along(object$lambda), drop = TRUE, ...) {
  if (!inherits(object, "BRIER.eta")) {
    stop("Object must be of class 'BRIER.eta', got '", class(object)[1], "'.", call. = FALSE)
  }
  if (!is.null(lambda)) {
    if (any(lambda > max(object$lambda) | lambda < min(object$lambda))) {
      stop("lambda must lie within the range of the fitted coefficient path.", call. = FALSE)
    }
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    # beta <- (1 - w) * object$beta[, l, drop = FALSE] + w * object$beta[, r, drop = FALSE]
    beta <- sweep(object$beta[, l, drop = FALSE], 2, (1 - w), "*") +
        sweep(object$beta[, r, drop = FALSE], 2, w, "*")
    colnames(beta) <- round(lambda, 4)
  } else {
    beta <- object$beta[, which, drop = FALSE]
  }
  if (drop) return(drop(beta)) else return(beta)
}


#' Predict from a single BRIER fit
#'
#' Generate predictions from a fitted \code{BRIER.eta} object at one or more
#' lambda values. Handles both individual-level fits (with intercept) and
#' summary-statistic fits (without intercept) automatically.
#'
#' @param object An object of class \code{"BRIER.eta"}.
#' @param X A numeric matrix of predictors. Required for \code{type} in
#'   \code{c("link", "response")}.
#' @param type A string specifying the prediction type:
#'   \itemize{
#'     \item \code{"link"}: linear predictor on the link scale.
#'     \item \code{"response"}: prediction on the response scale (after inverse link).
#'     \item \code{"coefficients"}: extract coefficients (delegates to \code{coef.BRIER.eta}).
#'     \item \code{"vars"}: indices of selected variables at each lambda.
#'     \item \code{"nvars"}: number of selected variables at each lambda.
#'   }
#' @param lambda Optional numeric vector of lambda values for prediction.
#'   If \code{NULL} (default), \code{which} is used to select lambdas
#'   directly.
#' @param which Integer vector specifying which lambda indices to use.
#'   Defaults to all fitted lambdas. Ignored when \code{lambda} is supplied.
#' @param ... Unused; present for S3 method compatibility.
#'
#' @return A numeric vector or matrix of predictions, depending on \code{type}.
#'
#' @seealso \code{\link{coef.BRIER.eta}}, \code{\link{predict.BRIER}}
#'
#' @export
predict.BRIER.eta <- function(
  object, X,
  type = c("link", "response", "coefficients", "vars", "nvars"),
  lambda = NULL, which = seq_along(object$lambda),
  ...
) {
  if (!inherits(object, "BRIER.eta")) {
    stop("Object must be of class 'BRIER.eta'.", call. = FALSE)
  }

  type <- match.arg(type)
  beta <- coef.BRIER.eta(object, lambda = lambda, which = which, drop = FALSE)

  if (type == "coefficients") return(beta)
  if (type == "nvars")        return(object$k[which])
  if (type == "vars")         return(drop(apply(abs(beta) >= 1e-8, 2, base::which)))

  # all remaining types need X
  if (missing(X) || is.null(X)) {
    stop("X must be supplied for type = '", type, "'.", call. = FALSE)
  }
  if (!is.matrix(X)) {
    X <- tryCatch(
      model.matrix(~ 0 + ., data = X),
      error = function(e) stop("X must be a matrix or coercible to one.", call. = FALSE)
    )
  }

  # extract intercept
  if (nrow(beta) == ncol(X) + 1) {
    alpha <- beta[1, ]
    beta  <- beta[-1, , drop = FALSE]
  } else if (nrow(beta) == ncol(X)) {
    alpha <- 0
  } else {
    stop(paste0(
      "beta has ", nrow(beta), " rows but X has ", ncol(X), " columns. Dimensions do not match."
    ), call. = FALSE)
  }

  eta <- sweep(X %*% beta, 2, alpha, "+")
  if (type == "link") return(drop(eta))

  resp <- switch(object$family,
    gaussian = eta,
    binomial = exp(eta) / (1 + exp(eta)),
    poisson  = exp(eta)
  )
  drop(resp)
}


#' Coefficients from a BRIER object
#'
#' Extract coefficients from a fitted \code{BRIER} object (output from
#' \code{\link{BRIERi}}, \code{\link{BRIERs}}, \code{\link{BRIERfull}}, or
#' \code{\link{BRIERi.cv}}, optionally after running selection).
#' Defaults to the optimal eta and lambda if selection has been performed.
#'
#' @param object An object of class \code{"BRIER"}.
#' @param eta Optional matrix of eta combinations to look up in
#'   \code{object$eta.grid} (each row matched to the nearest fitted eta within a
#'   relative tolerance of 1\%; see \code{\link{subsetEta}}). If \code{NULL}
#'   (default), \code{which.eta} is used to select etas directly.
#' @param which.eta Optional integer vector of row indices into
#'   \code{object$eta.grid}. Defaults to \code{object$eta.min.index}.
#' @param lambda Optional numeric vector of lambda values at which to extract
#'   coefficients. If \code{NULL} (default), \code{which.lambda} is used.
#'   Only meaningful when a single \code{which.eta} is selected.
#' @param which.lambda Optional integer vector of lambda indices. If omitted, the
#'   best lambda for the requested eta: \code{object$lambda.min.index} for the
#'   selected eta, and that eta's own entry in \code{object$eta.lambda} for any
#'   other eta of a tuned object.
#' @param drop Logical. If TRUE, drop singleton dimensions from the output.
#' @param ... Unused; present for S3 method compatibility.
#'
#' @return If a single eta is requested, a numeric matrix or vector of
#'   coefficients. If multiple etas are requested, a named list of coefficient
#'   matrices, one per eta combination.
#'
#' @seealso \code{\link{coef.BRIER.eta}}, \code{\link{predict.BRIER}}
#'
#' @export
coef.BRIER <- function(
  object,
  eta = NULL, which.eta = object$eta.min.index,
  lambda = NULL, which.lambda = object$lambda.min.index,
  drop = TRUE,
  ...
) {

  if (!inherits(object, "BRIER")) {
    stop("Object must be of class 'BRIER', got '", class(object)[1], "'.", call. = FALSE)
  }

  # -- Resolve eta --
  if (!is.null(eta)) {
    eta <- as.matrix(eta)
    if (ncol(eta) != ncol(object$eta.grid)) {
      stop(paste0(
        "eta must have ", ncol(object$eta.grid), " columns (one per external model). ",
        "Got ", ncol(eta), "."
      ), call. = FALSE)
    }
    which.eta <- .eta_rows(object, eta)
  }

  # -- Validate which.eta --
  if (is.null(which.eta) || length(which.eta) == 0) {
    stop("No eta selected. Please supply 'eta' or 'which.eta', ",
         "or run selection/CV first.", call. = FALSE)
  }
  if (any(which.eta < 1) || any(which.eta > length(object$res))) {
    stop(paste0("which.eta must be between 1 and ", length(object$res), "."), call. = FALSE)
  }

  # -- Capture lambda missingness --
  has.lambda <- !is.null(lambda)

  # -- Single eta: use lambda/which.lambda directly --
  if (length(which.eta) == 1) {
    # The default lambda belongs to the SELECTED eta. For any other eta, use that eta's
    # own best lambda from the selection (eta.lambda), as the multi-eta branch does:
    # otherwise eta = 0 would be predicted at the integrated model's lambda index.
    if (!has.lambda && missing(which.lambda)) {
      which.lambda <- .lambda_for_eta(object, which.eta, which.lambda)
    }
    if (!has.lambda && (is.null(which.lambda) || length(which.lambda) == 0)) {
      stop("No lambda selected. Please supply 'lambda' or 'which.lambda'.", call. = FALSE)
    }
    fit <- object$res[[which.eta]]
    if (has.lambda) {
      return(coef.BRIER.eta(fit, lambda = lambda, drop = drop, ...))
    } else {
      return(coef.BRIER.eta(fit, which = which.lambda, drop = drop, ...))
    }
  }

  # -- Multiple eta: look up best lambda per eta from eta.lambda --
  if (is.null(object$eta.lambda)) {
    stop(
      "Multiple eta requested but object$eta.lambda not found. ",
      "Run selection or CV first.", call. = FALSE
    )
  }

  out <- lapply(which.eta, function(i) {
    lam.idx <- object$eta.lambda$lambda.min.index[object$eta.lambda$eta.index == i]
    if (length(lam.idx) == 0) {
      stop(paste0("No eta.lambda entry found for eta.index = ", i, "."), call. = FALSE)
    }
    fit <- object$res[[i]]
    coef.BRIER.eta(fit, which = lam.idx[1], drop = drop, ...)
  })

  names(out) <- apply(object$eta.grid[which.eta, , drop = FALSE], 1, function(row) {
    paste0("(", paste(round(row, 4), collapse = ", "), ")")
  })
  out
}


#' Predict from a BRIER object
#'
#' Generate predictions from a fitted \code{BRIER} object. Defaults to the
#' optimal eta and lambda if selection has been performed.
#'
#' @param object An object of class \code{"BRIER"}.
#' @param X A numeric matrix of predictors.
#' @param eta Optional matrix of eta combinations to look up in
#'   \code{object$eta.grid} (each row matched to the nearest fitted eta within a
#'   relative tolerance of 1\%; see \code{\link{subsetEta}}). If \code{NULL}
#'   (default), \code{which.eta} is used to select etas directly.
#' @param which.eta Optional integer vector of row indices into
#'   \code{object$eta.grid}. Defaults to \code{object$eta.min.index}.
#' @param lambda Optional numeric vector of lambda values for prediction.
#'   If \code{NULL} (default), \code{which.lambda} is used. Only meaningful
#'   when a single \code{which.eta} is selected.
#' @param which.lambda Optional integer vector of lambda indices. If omitted, the
#'   best lambda for the requested eta: \code{object$lambda.min.index} for the
#'   selected eta, and that eta's own entry in \code{object$eta.lambda} for any
#'   other eta of a tuned object.
#' @param type A string: "link", "response", "coefficients", "vars", or "nvars".
#'   See \code{\link{predict.BRIER.eta}}.
#' @param drop Logical. If TRUE, drop singleton dimensions from the output.
#' @param ... Unused; present for S3 method compatibility.
#'
#' @return If a single eta is requested, a numeric vector or matrix. If
#'   multiple etas are requested, a named list of predictions.
#'
#' @seealso \code{\link{predict.BRIER.eta}}, \code{\link{coef.BRIER}}
#'
#' @export
predict.BRIER <- function(
  object, X,
  eta = NULL, which.eta = object$eta.min.index,
  lambda = NULL, which.lambda = object$lambda.min.index,
  type = c("link", "response", "coefficients", "vars", "nvars"),
  drop = TRUE,
  ...
) {

  if (!inherits(object, "BRIER")) {
    stop("Object must be of class 'BRIER', got '", class(object)[1], "'.", call. = FALSE)
  }

  type <- match.arg(type)

  # -- Resolve eta --
  if (!is.null(eta)) {
    eta <- as.matrix(eta)
    if (ncol(eta) != ncol(object$eta.grid)) {
      stop(paste0(
        "eta must have ", ncol(object$eta.grid), " columns (one per external model). ",
        "Got ", ncol(eta), "."
      ), call. = FALSE)
    }
    which.eta <- .eta_rows(object, eta)
  }

  # -- Validate which.eta --
  if (is.null(which.eta) || length(which.eta) == 0) {
    stop(
      "No eta selected. Please supply 'eta' or 'which.eta', ",
      "or run selection/CV first.", call. = FALSE
    )
  }
  if (any(which.eta < 1) || any(which.eta > length(object$res))) {
    stop(paste0("which.eta must be between 1 and ", length(object$res), "."), call. = FALSE)
  }

  # -- Capture lambda missingness --
  has.lambda <- !is.null(lambda)

  # -- Single eta: use lambda/which.lambda directly --
  if (length(which.eta) == 1) {
    # The default lambda belongs to the SELECTED eta. For any other eta, use that eta's
    # own best lambda from the selection (eta.lambda), as the multi-eta branch does:
    # otherwise eta = 0 would be predicted at the integrated model's lambda index.
    if (!has.lambda && missing(which.lambda)) {
      which.lambda <- .lambda_for_eta(object, which.eta, which.lambda)
    }
    if (!has.lambda && (is.null(which.lambda) || length(which.lambda) == 0)) {
      stop("No lambda selected. Please supply 'lambda' or 'which.lambda'.", call. = FALSE)
    }
    fit <- object$res[[which.eta]]
    if (has.lambda) {
      return(predict.BRIER.eta(fit, X = X, lambda = lambda, type = type, drop = drop, ...))
    } else {
      return(predict.BRIER.eta(fit, X = X, which = which.lambda, type = type, drop = drop, ...))
    }
  }

  # -- Multiple eta: look up best lambda per eta from eta.lambda --
  if (is.null(object$eta.lambda)) {
    stop(
      "Multiple eta requested but object$eta.lambda not found. ",
      "Run selection or CV first.", call. = FALSE
    )
  }

  out <- lapply(which.eta, function(i) {
    lam.idx <- object$eta.lambda$lambda.min.index[object$eta.lambda$eta.index == i]
    if (length(lam.idx) == 0) {
      stop(paste0("No eta.lambda entry found for eta.index = ", i, "."), call. = FALSE)
    }
    fit <- object$res[[i]]
    predict.BRIER.eta(fit, X = X, which = lam.idx[1], type = type, drop = drop, ...)
  })

  names(out) <- apply(object$eta.grid[which.eta, , drop = FALSE], 1, function(row) {
    paste0("(", paste(round(row, 4), collapse = ", "), ")")
  })
  out
}

# Rows of object$eta.grid matching each row of `eta`. A typed value need not equal the
# fitted one to the last digit (3.59 for 3.5938...): each row is matched to the NEAREST
# fitted row, accepted when every component's relative difference |a - b| / max(|a|, |b|)
# is within `tol`. So 0 matches only 0, and a small value never snaps to a different one
# (0.005 is not 0.01). Nothing unfitted is ever returned: a value with no fitted row
# within tolerance is an error that lists the fitted grid.
.eta_rows <- function(object, eta, tol = 1e-2) {
  grid <- as.matrix(object$eta.grid)
  eta <- as.matrix(eta)
  if (ncol(eta) != ncol(grid)) {
    if (ncol(grid) == 1) {
      eta <- matrix(as.numeric(eta), ncol = 1)
    } else {
      stop(paste0(
        "eta must have ", ncol(grid), " columns (one per external model). ",
        "Got ", ncol(eta), "."
      ), call. = FALSE)
    }
  }
  apply(eta, 1, function(e) {
    diff <- abs(sweep(grid, 2, e, "-"))
    denom <- pmax(abs(grid), matrix(abs(e), nrow(grid), ncol(grid), byrow = TRUE))
    rel <- ifelse(denom == 0, 0, diff / denom)
    dist <- apply(rel, 1, max)
    i <- which.min(dist)
    if (length(i) == 0 || dist[i] > tol) {
      stop(
        "eta = (", paste(signif(e, 4), collapse = ", "), ") does not match any row in ",
        "eta.grid (it was not fitted). Fitted eta values: ",
        paste(apply(signif(grid, 4), 1, paste, collapse = ", "), collapse = "; "),
        ". Refit with it in eta.list, or use one of these.", call. = FALSE
      )
    }
    i
  })
}

# The best lambda index for one eta row of a tuned object: its own entry in eta.lambda.
# For the selected eta this equals object$lambda.min.index (every selection function
# sets lambda.min.index from eta.lambda), so the default is unchanged there. An untuned
# object (no eta.lambda) keeps the given default.
.lambda_for_eta <- function(object, which.eta, default) {
  el <- object$eta.lambda
  if (is.null(el) || is.null(el$eta.index)) {
    return(default)
  }
  idx <- el$lambda.min.index[el$eta.index == which.eta]
  if (length(idx) == 0) {
    warning(
      "eta row ", which.eta, " has no entry in eta.lambda; using lambda.min.index, ",
      "the best lambda of the SELECTED eta. Pass which.lambda to choose one.",
      call. = FALSE
    )
    return(default)
  }
  idx[1]
}
