#' BRIER.I: integration of pretrained external models with target individual-level data
#'
#' Fits a penalized regression model that integrates one or more pretrained
#' external prediction models with target individual-level data. External
#' predictions are blended into the target response via a tunable eta parameter
#' (one per external model), and lambda is selected over a regularisation path.
#'
#' @param X A numeric matrix of predictors (n x p) for the target cohort.
#' @param y A numeric response vector of length n.
#' @param family A string specifying the response distribution: "gaussian",
#'   "binomial", or "poisson".
#' @param eta.list A list of numeric vectors (one per external model) specifying
#'   the eta grid for that model. A single numeric vector is replicated for
#'   all external models with a warning.
#' @param beta.external A numeric matrix of external model coefficients
#'   ((p+1) x M, where M is the number of external models). Must include the
#'   intercept as the first row.
#' @param multi.method A string specifying how to combine multiple external models:
#'   \itemize{
#'     \item \code{"ind"}: keep external models independent (default).
#'     \item \code{"PCA"}: aggregate via the leading principal component of
#'       \code{beta.external}.
#'     \item \code{"stacking"}: aggregate via stacking weights estimated from
#'       the target data.
#'   }
#' @param optim.args A list of arguments passed to \code{optim()} when using
#'   stacking with binomial or Poisson families.
#' @param ... Additional arguments passed to \code{BRIERi_fit} (e.g.
#'   \code{penalty}, \code{alpha}, \code{nlambda}, \code{penalty.factor}).
#' @param trace Logical. If TRUE, prints progress messages during fitting.
#' @param ncores Integer. Number of cores for parallel fitting.
#' @param parallel Logical. If TRUE and on a non-Windows platform, fits the eta
#'   grid in parallel using \code{parallel::mclapply}.
#'
#' @return An object of class \code{"BRIER"} containing:
#' \describe{
#'   \item{y}{The response vector.}
#'   \item{y.external}{A matrix of external predictions (n x M, or n x 1 after
#'     PCA/stacking aggregation).}
#'   \item{family}{The response family.}
#'   \item{eta.list}{The list of per-model eta grids.}
#'   \item{eta.grid}{The full combinatorial eta grid (matrix, rows are combinations).}
#'   \item{res}{A list of \code{BRIER.fit} objects, one per eta combination.}
#'   \item{null.dev}{The null deviance.}
#' }
#'
#' @seealso \code{\link{BRIERi_fit}}, \code{\link{BRIERi.cv}},
#'   \code{\link{BRIERi.selection}}, \code{\link{BRIERfull}},
#'   \code{\link{BRIERs}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 200
#' p <- 50
#' X <- matrix(rnorm(n * p), ncol = p)
#' beta_true <- c(rep(1, 5), rep(0, p - 5))
#' y <- X %*% beta_true + rnorm(n)
#'
#' # External model coefficients (intercept + p variables)
#' beta.external <- matrix(c(0, beta_true * 0.8), ncol = 1)
#'
#' fit <- BRIERi(
#'   X, y, family = "gaussian",
#'   eta.list = list(c(0, 0.5, 1, 2)),
#'   beta.external = beta.external,
#'   penalty = "LASSO"
#' )
#' }
#'
#' @export
BRIERi = function(
  X, y, family = c("gaussian", "binomial", "poisson"), eta.list, beta.external,
  multi.method = c("ind", "PCA", "stacking"), optim.args = list(),
  ...,
  trace = FALSE,
  ncores = max(1L, parallel::detectCores() - 1L),
  parallel = (ncores > 1L)
){

  family <- match.arg(family)
  multi.method <- match.arg(multi.method)

  # -- Validation of response y --

  ## Coerce y to matrix
  if (is.data.frame(y)){ y <- as.matrix(y) }
  if (is.vector(y)){ y <- matrix(y) }
  if (ncol(y) > 1) { stop("y must be a vector or a single column matrix", call. = FALSE) }
  if (any(is.na(y))) { stop("Missing data (NA's) detected in y.", call. = FALSE) }

  ## Convert fuzzy binomial data
  if (family == "binomial" && typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2) { stop("Attemping to use family = 'binomial' with non-binary data", call. = FALSE) }
    if (!identical(names(tab), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(y = ", names(tab)[2], ")"))
      y <- as.double(as.character(y) == names(tab)[2])
    } 
  }

  ## Convert to double, if necessary
  if (typeof(y) != "double") {
    tryCatch(
      storage.mode(y) <- "double", 
      warning = function(w){
        stop("y must be numeric or able to be coerced to numeric", call. = FALSE)
      }
    )
  }

  # -- Validation of design X --

  ## Coerce X to matrix
  if (!inherits(X, "matrix")) {
    tmp <- try(X = model.matrix( ~ 0 + ., data=X), silent = TRUE)
    if (inherits(tmp, "try-error")) { stop("X must be a matrix or able to be coerced to a matrix", call. = FALSE) }
  }
  if (storage.mode(X) == "integer") { storage.mode(X) <- "double"} 
  if (any(is.na(X))) { stop("Missing data (NA's) detected in X.", call. = FALSE) }

  ## Checking X and y dimensions
  if (nrow(X) != nrow(y)) { stop("X and y do not have the same number of observations", call. = FALSE) }

  # -- Validation of external info --

  ## Calculate external information from external models
  beta.external <- as.matrix(beta.external)
  if (any(is.na(beta.external))) { stop("Missing data (NA's) detected in beta.external.", call. = FALSE) }
  if (nrow(beta.external) != ncol(X) + 1) { stop("The dimension of external beta and X does not match. Please include intercept.") }
  ext <- calcExtY(X, y, beta.external, family, multi.method, optim.args)
  y.external <- ext$y.external
  M <- ncol(y.external)

  # -- Validate eta.list --
  if (!is.list(eta.list)) {
    if (!is.numeric(eta.list)) { stop("eta.list must be a numeric vector or a list of numeric vectors.", call. = FALSE) }
    warning("eta.list is not a list. Replicating the same eta grid for all ", M, " external model(s).", call. = FALSE)
    eta.list <- rep(list(eta.list), M)
  }
  if (length(eta.list) != M) {
    stop(paste0("eta.list has ", length(eta.list), " element(s) but there are ", M, " external model(s)."), call. = FALSE)
  }
  if (M >= 5) {
    warning(
      "5 or more external models provided. Please consider using ensemble ",
      "weighting to aggregate multiple external models.", call. = FALSE
    )
  }
  for (k in seq_len(M)) {
    if (!is.numeric(eta.list[[k]])) { stop(paste0("eta.list[[", k, "]] must be numeric."), call. = FALSE) }
  }

  # -- Build eta grid --
  eta.grid <- as.matrix(expand.grid(eta.list))
  colnames(eta.grid) <- paste0("eta_", seq_len(M))
  if (nrow(eta.grid) > 1000) {
    warning("Large eta grid (", nrow(eta.grid), " combinations). Consider coarser eta.list grids.", call. = FALSE)
  }
  eta.grid.list <- lapply(seq_len(nrow(eta.grid)), function(i) eta.grid[i, ])
  if (trace) { cat("Total eta combinations:", nrow(eta.grid), "\n") }

  # -- BRIER model fitting --

  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$family <- family 
  fit.args$y.external <- y.external

  ## Calculate null deviance
  penalty.factor <- if (!is.null(fit.args$penalty.factor)) {
    fit.args$penalty.factor
  } else {
    rep(1, ncol(X)) 
  }
  null.dev <- calcNullDev(X, y, rep(1, nrow(X)) / nrow(X), penalty.factor, family)

  ## Function to run single eta vec
  single_eta_fit <- function(eta, fit.args, trace) {
    if (trace) { cat("Fitting at eta =", paste(round(eta, 3), collapse = ", "), "\n") }
    fit.args$eta <- as.numeric(eta)  # vector of length K
    do.call(BRIERi_fit, fit.args)
  }

  ## Parallel fitting BRIER
  use_parallel <- isTRUE(parallel) && ncores > 1L && .Platform$OS.type != "windows"
  if (.Platform$OS.type == "windows" && ncores > 1L)
    message("Parallel execution not supported on Windows; using serial.")
  if (!use_parallel) {
    res <- lapply(eta.grid.list, FUN = single_eta_fit, fit.args = fit.args, trace = trace)
  } else {
    res <- parallel::mclapply(
      eta.grid.list,
      FUN      = single_eta_fit,
      fit.args = fit.args,
      trace    = trace,
      mc.cores = ncores
    )
  }

  out <- list(
    y = y,
    y.external = y.external,
    family = family,
    eta.list = eta.list,
    eta.grid = eta.grid,
    res = res,
    null.dev = null.dev
  )
  class(out) <- "BRIER"
  out
}


#' Core fitting function for BRIER.I
#'
#' Internal workhorse that fits a penalized regression model for a single eta
#' vector. Called by \code{\link{BRIERi}}, \code{\link{BRIERfull}}, and
#' \code{\link{BRIERi.cv}}. Performs response/predictor validation,
#' standardisation, lambda path setup, and dispatches to the appropriate
#' coordinate-descent C++ routine.
#'
#' @param X A numeric matrix of predictors (n x p).
#' @param y A numeric response vector of length n.
#' @param weights Optional numeric vector of observation weights (length n).
#'   Defaults to equal weights.
#' @param eta A numeric vector of external-model weights (length M).
#'   Use \code{eta = 0} to ignore external information.
#' @param y.external Optional numeric matrix of external predictions (n x M).
#' @param family A string: "gaussian", "binomial", or "poisson".
#' @param penalty.factor A numeric vector of per-variable penalty multipliers
#'   (length p). Variables with multiplier 0 are unpenalised.
#' @param penalty A string specifying the penalty: "LASSO", "MCP", or "SCAD".
#' @param alpha A numeric scalar for the elastic net mixing parameter.
#' @param gamma A numeric scalar for MCP and SCAD penalties.
#' @param n.sis Optional. Reserved for future SIS extensions.
#' @param lambda Optional user-specified lambda sequence.
#' @param nlambda An integer specifying the number of lambda values.
#' @param lambda.min A numeric scalar for the minimum lambda ratio.
#' @param log.lambda Logical. If TRUE, the lambda sequence is log-spaced.
#' @param max.iter An integer specifying the maximum total iterations.
#' @param eps A numeric scalar for the convergence tolerance.
#' @param nvar.max An integer specifying the maximum number of selected variables.
#' @param returnX Logical. If TRUE, the standardised data list is returned in
#'   the output for diagnostics.
#'
#' @return An object of class \code{"BRIER.fit"} containing model specification,
#'   coefficients, fit statistics, and data.
#'
#' @seealso \code{\link{BRIERi}}, \code{\link{BRIERfull}},
#'   \code{\link{coef.BRIER.fit}}, \code{\link{predict.BRIER.fit}}
#'
#' @keywords internal
#' @export
BRIERi_fit = function(
  X, y, weights = NULL, eta = 0, y.external = NULL, family = c("gaussian", "binomial", "poisson"), 
  penalty.factor = rep(1, ncol(X)), penalty = c("LASSO", "SCAD", "MCP"), 
  alpha = 1, gamma = ifelse(penalty == "SCAD", 3.7, 3), n.sis = NULL, 
  lambda, nlambda = 100, lambda.min = { if (nrow(X) > ncol(X)) 1e-4 else .05 }, log.lambda = TRUE, 
  max.iter = 1e6, eps = 1e-4, nvar.max = NULL, returnX = FALSE
){

  family <- match.arg(family)
  penalty <- match.arg(penalty)

  # -- Validation of y --
  if (is.data.frame(y)){ y <- as.matrix(y) }
  if (is.vector(y)){ y <- matrix(y) }
  if (ncol(y) > 1) stop("y must be a vector or a single column matrix", call. = FALSE)
  if (any(is.na(y))) stop("Missing data (NA's) detected in y.  Please eliminate missing data (e.g., by removing cases, removing features, or imputation)", call.=FALSE)  
  if (family == "binomial" && typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2) stop("Attemping to use family = 'binomial' with non-binary data", call. = FALSE)
    if (!identical(names(tab), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(y = ", names(tab)[2], ")"))
      y <- as.double(as.character(y) == names(tab)[2])
    } 
  } ## Convert fuzzy binomial data
  if (typeof(y) != "double") {
    tryCatch(
      storage.mode(y) <- "double", 
      warning = function(w){ stop("y must be numeric or able to be coerced to numeric", call. = FALSE) }
    )
  } ## Convert to double, if necessary

  # -- Validation of X --
  if (!inherits(X, "matrix")) {
    tmp <- try(X = model.matrix( ~ 0 + ., data = X), silent = TRUE)
    if (inherits(tmp, "try-error")) { stop("X must be a matrix or able to be coerced to a matrix", call. = FALSE) }
  }
  if (storage.mode(X) == "integer") storage.mode(X) <- "double"
  if (any(is.na(X))) { stop("Missing data (NA's) detected in X.", call.=FALSE) }
  if (nrow(X) != nrow(y)) { stop("X and y do not have the same number of observations", call. = FALSE) }

  # -- Validation of external info --
  eta <- as.vector(eta)
  eta <- as.numeric(eta)
  if (any(!is.finite(eta))) { stop("eta must be finite (no NA/Inf).", call. = FALSE) }
  if (any(eta < 0)) { stop("eta must be >= 0.", call. = FALSE) }
  if (!is.null(y.external)) {                              
    y.external <- as.matrix(y.external) 
    if (nrow(y) != nrow(y.external)) { stop("y and y.external do not have the same number of observations", call. = FALSE) }
    if (length(eta) != ncol(y.external)) { stop("Number of external models and length of eta do not match", call. = FALSE) }
  } else {
    if (any(eta != 0)) { stop("y.external must be provided when any eta != 0.", call. = FALSE) }
  }                                                        

  # -- Normalize weight --
  if(is.null(weights)) { weights <- rep(1, nrow(X)) }
  weights <- as.numeric(weights)
  if (length(weights) != nrow(X)){
    stop(paste0(
      "Number of elements in weights (", length(weights), 
      ") not equal to the number of rows of x (", nrow(X), ")"
    ))
  }
  if (any(!is.finite(weights))) { stop("Observation weights must be finite (no NA/Inf).") }
  if (any(weights < 0)) { stop("Observation weights must be >= 0.") }
  if (all(weights == 0)) { stop("At least one weight must be > 0.") }
  X <- X[weights != 0, , drop = FALSE]
  y <- y[weights != 0, , drop = FALSE]
  if (!is.null(y.external)) { y.external <- y.external[weights != 0, , drop = FALSE] }
  weights <- weights[weights != 0]
  weights <- weights / sum(weights)

  # -- Checking penalty factor --
  if (!is.double(penalty.factor)) { penalty.factor <- as.double(penalty.factor) }
  if (any(is.na(penalty.factor))) { stop("Missing data (NA's) detected in penalty.factor.", call. = FALSE) }
  if (any(penalty.factor < 0)) { stop("penalty.factor must be non-negative", call. = FALSE) }
  if (length(penalty.factor) != ncol(X)) { stop("Length of penalty.factor must be the same as the number of columns in X", call. = FALSE) }

  # -- Checking sparse parameters --
  if (nlambda < 2) {
    stop("nlambda must be at least 2", call. = FALSE)
  } else if (nlambda != round(nlambda)){
    stop("nlambda must be a positive integer", call. = FALSE)
  }
  if (alpha <= 0 || alpha > 1)  { stop("alpha must be between [0, 1]; choose a small positive number instead", call. = FALSE) }
  if (gamma <= 1 && penalty == "MCP") { stop("gamma must be greater than 1 for the MCP penalty", call. = FALSE) }
  if (gamma <= 2 && penalty == "SCAD") { stop("gamma must be greater than 2 for the SCAD penalty", call. = FALSE) }

  initial.p <- ncol(X)
  initial.colnames <- if (is.null(colnames(X))) {
    paste0("V", seq_len(ncol(X)))
  } else {
    colnames(X)
  }

  # -- Process X and y --
  # Remove constant covariates in X; 
  # Perform standardization of X and centering of y (gaussian)
  # Perform sure independence screening if necessary
  XX.list <- newData(X, y, eta, y.external, weights, family, penalty.factor)
  n <- nrow(XX.list$std.X)
  p <- ncol(XX.list$std.X)
  if (is.null(nvar.max)) nvar.max <- p

  # -- Main algorithm BRIERi_fit  --
  if (missing(lambda)) {
    lambda.seq <- setupLambda(
      XX.list$std.X, XX.list$yy, n, p, weights, XX.list$penalty.factor, alpha,
      family, lambda.min, log.lambda, nlambda
    )
    lam.max <- lambda.seq[1]
  } else {
    nlambda <- length(lambda)  # Note: lambda can be a single value
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
    lam.max <- -1
  }

  # -- Fit main algorithm  --
  if (family == "gaussian"){
    fit <- cd_wgaussian_fit_ssr(
      XX.list$std.X, XX.list$yy, weights, 
      XX.list$penalty.factor, penalty, lambda.seq, alpha, gamma,
      lam.max, max.iter, eps, nvar.max
    )
    b0 <- fit$beta0 + XX.list$yy.center
    b <- fit$beta
    df <- fit$df
    iter <- fit$iter
  } else {
    fit <- cd_wglm_fit_ssr(
      XX.list$std.X, XX.list$yy, weights, family, 
      XX.list$penalty.factor, penalty, lambda.seq, alpha, gamma,
      lam.max, max.iter, nvar.max, eps
    )
    b0 <- fit$beta0
    b <- fit$beta
    df <- fit$df
    iter <- fit$iter
  }

  # -- Eliminate saturated lambda values --
  ind <- !is.na(iter) & colSums(!is.finite(b)) == 0 
  b0 <- b0[ind, , drop = FALSE]
  b <- b[, ind, drop = FALSE]
  lambda <- lambda.seq[ind]
  df <- df[ind]
  iter <- iter[ind]
  if (iter[1] == max.iter){
    stop("Algorithm failed to converge for any values of lambda", call. = FALSE)
  }
  if (any(iter == max.iter)){
    warning("Algorithm failed to converge for some values of lambda", call. = FALSE)
  }

  # -- Unstandardize --
  beta <- matrix(0, nrow = (ncol(X) + 1), ncol = length(lambda))
  bb <- b/XX.list$scale[XX.list$nz]
  beta[XX.list$nz+1, ] <- bb
  beta[1, ] <- b0 - as.vector(crossprod(XX.list$center[XX.list$nz], bb))
  varname <- c("Intercept", initial.colnames)
  dimnames(beta) <- list(varname, round(lambda, digits = 4))

  # -- Deviance --
  linear.predictor <- sweep(X %*% beta[-1, , drop = FALSE], 2, beta[1, ], "+")
  colnames(linear.predictor) <- round(lambda, digits = 4)
  mu <- ginv_link(linear.predictor, family)
  deviance <- colSums(calcDev(y, mu, weights, family))

  # -- Degree of freedom --
  k <- colSums(abs(beta) >= 1e-8) 
  df <- df + 1

  out <- list(
    # Model Specification
    family = family,
    penalty = penalty,
    penalty.factor = penalty.factor,
    eta = eta,
    lambda = lambda,
    alpha  = alpha,
    gamma  = gamma,

    # Coefficients
    # beta0 = b0,
    beta  = beta,
    df = df, 
    k = k,

    # Fit Statistics
    deviance = deviance,
    linear.predictor = linear.predictor,
    iter = iter,

    # Data
    y = y,
    y.tilde = XX.list$y.tilde,
    y.external = XX.list$y.external,
    weights = weights
  )
  if (returnX){
    out$XX.list <- XX.list
  }
  class(out) <- "BRIER.fit"
  out
}


calcExtY = function(
  X, y, beta.external, family, 
  multi.method, 
  optim.args = list()
  ){
  if (multi.method == "PCA"){ 
    bb = apply(beta.external, 2, function(x) x / sqrt(sum(x^2)))
    w = prcomp(t(bb) %*% bb)$rotation[, 1]
    beta.external = beta.external %*% abs(w)
  } 
  pred = X %*% beta.external[-1, , drop = FALSE] + beta.external[1, ]
  y.external = ginv_link(pred, family)
  if (multi.method == "stacking"){
    optim.args$Y = y.external
    optim.args$z = y
    if (family == "gaussian"){
      w = stacking_gaussian(y.external, y)
    } else if (family == "binomial"){ 
      fit = do.call(stacking_binomial, optim.args)
      w = fit$weights
    } else if (family == "poisson"){ 
      fit = do.call(stacking_poisson, optim.args)
      w = fit$weights
    }
    y.external = y.external %*% w 
  }
  return(list(
    beta.external = beta.external, 
    y.external = y.external
  ))
}

stacking_gaussian = function(Y, z){ solve(as.matrix(t(Y) %*% Y), crossprod(Y, z))}

stacking_binomial = function(
  Y, z, 
  w_init = NULL, method = "BFGS",
  maxit = 1e6, reltol = 1e-10,
  eps = 1e-6, big = 1e20
) {
  Y = as.matrix(Y)
  z = as.numeric(z)

  n = nrow(Y)
  m = ncol(Y)

  if (length(z) != n) {
    stop("length(z) must equal nrow(Y)", call. = FALSE)
  }

  if (any(!(z %in% c(0, 1)))) {
    stop("z must be binary 0/1")
  }

  if (is.null(w_init)) {
    # equal weights is often a safer start
    w_init = rep(1 / m, m)
  }

  obj = function(w, Y, z, eps, big) {
    mu = as.vector(Y %*% w)

    if (any(mu <= eps) || any(mu >= 1 - eps)) {
      return(big)
    }

    -sum(z * log(mu) + (1 - z) * log(1 - mu))
  }

  grad = function(w, Y, z, eps, big) {
    mu = as.vector(Y %*% w)

    if (any(mu <= eps) || any(mu >= 1 - eps)) {
      return(rep(NA_real_, ncol(Y)))
    }

    r = -(z - mu) / (mu * (1 - mu))
    as.vector(crossprod(Y, r))
  }

  fit = optim(
    par = w_init,
    fn = obj,
    gr = grad,
    Y = Y,
    z = z,
    eps = eps,
    big = big,
    method = method,
    control = list(maxit = maxit, reltol = reltol)
  )

  w_hat = fit$par
  mu_hat = as.vector(Y %*% w_hat)

  list(
    weights = w_hat,
    fitted = mu_hat,
    value = fit$value,
    convergence = fit$convergence,
    message = fit$message,
    optim = fit
  )
}

stacking_poisson = function(
  Y, z, 
  w_init = NULL, method = "BFGS",
  maxit = 1e6, reltol = 1e-10,
  eps = 1e-6, big = 1e20
  ) {
  Y = as.matrix(Y)
  z = as.numeric(z)

  n = nrow(Y)
  m = ncol(Y)

  if (length(z) != n) {
    stop("length(z) must equal nrow(Y)")
  }

  if (any(z < 0) || any(abs(z - round(z)) > 1e-8)) {
    stop("z must be nonnegative integer counts")
  }

  if (is.null(w_init)) {
    w_init = rep(1 / m, m)
  }

  obj = function(w, Y, z, eps, big) {
    mu = as.vector(Y %*% w)

    if (any(mu <= eps)) {
      return(big)
    }

    sum(mu - z * log(mu))
  }

  grad = function(w, Y, z, eps, big) {
    mu = as.vector(Y %*% w)

    if (any(mu <= eps)) {
      return(rep(NA_real_, ncol(Y)))
    }

    r = -(z - mu) / mu
    as.vector(crossprod(Y, r))
  }

  fit = optim(
    par = w_init,
    fn = obj,
    gr = grad,
    Y = Y,
    z = z,
    eps = eps,
    big = big,
    method = method,
    control = list(maxit = maxit, reltol = reltol)
  )

  w_hat = fit$par
  mu_hat = as.vector(Y %*% w_hat)

  list(
    weights = w_hat,
    fitted = mu_hat,
    value = fit$value,
    convergence = fit$convergence,
    message = fit$message,
    optim = fit
  )
}



