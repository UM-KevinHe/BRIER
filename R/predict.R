#' Coefficients from a single BRIER fit
#'
#' Extract coefficients from a fitted \code{BRIER.eta} object at one or more
#' lambda values. Lambda values not in the fitted path are obtained by linear
#' interpolation between adjacent fitted lambdas.
#'
#' @param object An object of class \code{"BRIER.eta"}.
#' @param lambda Optional numeric vector of lambda values at which to extract
#'   coefficients. Must lie within the fitted path. If missing, coefficients are
#'   returned at the indices given by \code{which}.
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
coef.BRIER.eta <- function(object, lambda, which = seq_along(object$lambda), drop = TRUE, ...) {
  if (!inherits(object, "BRIER.eta")) {
    stop("Object must be of class 'BRIER.eta', got '", class(object)[1], "'.", call. = FALSE)
  }
  if (!missing(lambda)) {
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
#' @param which Integer vector specifying which lambda indices to use.
#'   Defaults to all fitted lambdas.
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
  lambda, which = seq_along(object$lambda),
  ...
) {
  if (!inherits(object, "BRIER.eta")) {
    stop("Object must be of class 'BRIER.eta'.", call. = FALSE)
  }

  type <- match.arg(type)
  beta <- coef.BRIER.eta(object, lambda = lambda, which = which, drop = FALSE)

  if (type == "coefficients") return(beta)
  if (type == "nvars")        return(object$k[which])
  if (type == "vars")         return(drop(apply(abs(beta) >= 1e-8, 2, which)))

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
#' @param eta Optional matrix of eta combinations to look up (one row per
#'   combination, one column per external model).
#' @param which.eta Optional integer vector of row indices into
#'   \code{object$eta.grid}. Defaults to \code{object$eta.min.index}.
#' @param lambda Optional numeric vector of lambda values at which to extract
#'   coefficients. Used only when a single \code{which.eta} is selected.
#' @param which.lambda Optional integer vector of lambda indices. Defaults to
#'   \code{object$lambda.min.index}.
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
  eta = object$eta.min, which.eta = object$eta.min.index,
  lambda = object$lambda.min, which.lambda = object$lambda.min.index,
  drop = TRUE,
  ...
) {

  if (!inherits(object, "BRIER")) {
    stop("Object must be of class 'BRIER', got '", class(object)[1], "'.", call. = FALSE)
  }

  # -- Resolve eta --
  if (!missing(eta)) {
    eta <- as.matrix(eta)
    if (ncol(eta) != ncol(object$eta.grid)) {
      stop(paste0(
        "eta must have ", ncol(object$eta.grid), " columns (one per external model). ",
        "Got ", ncol(eta), "."
      ), call. = FALSE)
    }
    which.eta <- apply(eta, 1, function(e) {
      idx <- which(apply(object$eta.grid, 1, function(row) { all(abs(row - e) < 1e-10) }))
      if (length(idx) == 0) {
        stop("eta = (", paste(round(e, 4), collapse = ", "),
             ") does not match any row in eta.grid.", call. = FALSE)
      }
      idx[1]
    })
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
  has.lambda <- !missing(lambda)

  # -- Single eta: use lambda/which.lambda directly --
  if (length(which.eta) == 1) {
    if (is.null(which.lambda) || length(which.lambda) == 0) {
      if (!has.lambda || is.null(lambda)) {
        stop("No lambda selected. Please supply 'lambda' or 'which.lambda'.", call. = FALSE)
      }
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
    stop("Multiple eta requested but object$eta.lambda not found. ",
         "Run selection or CV first.", call. = FALSE)
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
#' @param eta Optional matrix of eta combinations to look up.
#' @param which.eta Optional integer vector of row indices into
#'   \code{object$eta.grid}.
#' @param lambda Optional numeric vector of lambda values for prediction.
#'   Used only when a single \code{which.eta} is selected.
#' @param which.lambda Optional integer vector of lambda indices.
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
  eta = object$eta.min, which.eta = object$eta.min.index,
  lambda = object$lambda.min, which.lambda = object$lambda.min.index,
  type = c("link", "response", "coefficients", "vars", "nvars"),
  drop = TRUE,
  ...
) {

  if (!inherits(object, "BRIER")) {
    stop("Object must be of class 'BRIER', got '", class(object)[1], "'.", call. = FALSE)
  }

  type <- match.arg(type)

  # -- Resolve eta --
  if (!missing(eta)) {
    eta <- as.matrix(eta)
    if (ncol(eta) != ncol(object$eta.grid)) {
      stop(paste0(
        "eta must have ", ncol(object$eta.grid), " columns (one per external model). ",
        "Got ", ncol(eta), "."
      ), call. = FALSE)
    }
    which.eta <- apply(eta, 1, function(e) {
      idx <- which(apply(object$eta.grid, 1, function(row) {
        all(abs(row - e) < 1e-10)
      }))
      if (length(idx) == 0) {
        stop("eta = (", paste(round(e, 4), collapse = ", "),
             ") does not match any row in eta.grid.", call. = FALSE)
      }
      idx[1]
    })
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
  has.lambda <- !missing(lambda)

  # -- Single eta: use lambda/which.lambda directly --
  if (length(which.eta) == 1) {
    if (is.null(which.lambda) || length(which.lambda) == 0) {
      if (!has.lambda || is.null(lambda)) {
        stop("No lambda selected. Please supply 'lambda' or 'which.lambda'.", call. = FALSE)
      }
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