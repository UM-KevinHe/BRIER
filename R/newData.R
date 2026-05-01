# -- Internal helper: process X and y --
#
# Standardize X (weighted), compute y.tilde (blended response), center y for
# Gaussian, and drop constant columns. Returns the standardized data along with
# centering/scaling info needed to unstandardize coefficients later.

newData <- function(X, y, eta, y.external, wt, family, penalty.factor) {

  ## Standardize X
  p <- ncol(X)
  XX.list <- wstandardize_X(X, wt)
  std.X <- XX.list[[1]]
  dimnames(std.X) <- dimnames(X)
  center <- XX.list[[2]]
  scale <- XX.list[[3]]

  ## Compute y.tilde
  if (any(eta < 0)) {
    stop("eta must be >= 0.", call. = FALSE)
  } else if (all(eta == 0)) {
    y.external <- NULL
    y.tilde    <- y
  } else {
    y.tilde <- (y + y.external %*% eta) / (1 + sum(eta))
  }

  if (family == "gaussian") {
    yy <- y.tilde - mean(y.tilde)
    yy.center <- mean(y.tilde)
  } else {
    yy <- y.tilde
    yy.center <- NULL
  }

  ## Remove constant columns according to standard deviation
  nz <- which(scale > 1e-6)  # non-constant columns
  if (length(nz) != ncol(X)) {
    std.X <- std.X[, nz, drop = FALSE]
    center <- center[nz]
    scale <- scale[nz]
    warning(paste("Constant columns are removed: ", paste(which(scale <= 1e-6), collapse = ", "), sep = ""))
  }

  list(
    std.X          = std.X,
    center         = center,
    scale          = scale,
    nz             = nz,
    yy             = yy,
    yy.center      = yy.center,
    y.tilde        = y.tilde,
    y.external     = y.external,
    penalty.factor = penalty.factor[nz]
  )
}


# -- Internal helper: sure independence screening (SIS) for GLM --

SIS_glm <- function(y, std.X, wt, family, penalty.factor, n.sis) {

  if (is.null(n.sis)) {
    return(list(idx = seq_len(ncol(std.X))))
  }

  penalized_idx <- which(penalty.factor > 0)
  unpenalized_idx <- which(penalty.factor == 0)

  if (family == "gaussian") {
    yy <- y - sum(wt * y) / sum(wt)
  } else {
    yy <- y
  }

  ## Calculate NULL log-likelihood
  if (length(unpenalized_idx) > 0) {
    fit0 <- glm(yy ~ std.X[, unpenalized_idx, drop = FALSE], weights = wt, family = family)
  } else {
    fit0 <- glm(yy ~ 1, weights = wt, family = family)
  }
  ll0 <- as.numeric(logLik(fit0))

  ## Calculate marginal log-likelihood for each penalized variable
  ll1 <- sapply(
    penalized_idx,
    function(j) {
      if (length(unpenalized_idx) > 0) {
        fit <- glm(yy ~ std.X[, c(unpenalized_idx, j), drop = FALSE], weights = wt, family = family)
      } else {
        fit <- glm(yy ~ std.X[, j, drop = FALSE], weights = wt, family = family)
      }
      as.numeric(logLik(fit))
    }
  )
  score <- 2 * (ll1 - ll0)
  k <- min(n.sis, length(penalized_idx))
  top_idx <- penalized_idx[order(score, decreasing = TRUE)][seq_len(k)]
  idx <- sort(c(unpenalized_idx, top_idx))

  list(
    score = score,
    idx   = sort(c(unpenalized_idx, top_idx))
  )
}