# =============================================================================
# df.R: degrees of freedom of a penalized fit, for the information criteria.
#
# THREE METHODS, chosen by `df.method` in BRIERi.eta / BRIERs.eta:
#
#   "active"      the number of nonzero coefficients. What Zhang, Li and Tsai (2010,
#                 JASA) use in their GIC, and what ncvreg uses; for the LASSO it is
#                 an unbiased estimate of df (Zou, Hastie and Tibshirani 2007).
#   "divergence"  the local (Stein) df with the active set held fixed: differentiating
#                 the fit's optimality conditions gives
#                   df = tr{ (H_A + diag p''(|b_A|))^{-1} H_A },
#                 with H_A the Hessian of the loss on the active set. p'' is 0 for the
#                 LASSO (so df = the active count exactly), -1/gamma for MCP below
#                 gamma*lambda, -1/(gamma - 1) for SCAD between lambda and
#                 gamma*lambda, plus the ridge part lambda*(1 - alpha). Nonconvex
#                 penalties shrink less, so their df exceeds the active count.
#   "lqa"         Fan and Li's (2001) local quadratic approximation, Zhang, Li and
#                 Tsai's (2010) df_L:  tr{ (H_A + diag p'(|b_A|)/|b_A|)^{-1} H_A }.
#                 Asymptotically the active count; below it in finite samples, the
#                 LASSO included.
#
# Everything is on the STANDARDIZED scale the fit was solved on. Where H_A + the
# curvature is not positive definite (a nonconvex penalty in its nonconvex region),
# the df is undefined there and that lambda falls back to the active count; the fit
# records how many did. BRIER's response adjustment is applied by the caller: the df
# with respect to y is this df divided by (1 + sum(eta)).
# =============================================================================

## The penalty's curvature for each active coefficient: p'' (divergence) or p'/|b|
## (lqa), plus the ridge part. lam1 = lambda * pf * alpha, lam2 = lambda * pf * (1 - alpha).
.df_curvature <- function(b, lam1, lam2, penalty, gamma, method) {
  a <- abs(b)
  pp <- if (identical(method, "divergence")) {
    switch(penalty,
      LASSO = rep(0, length(a)),
      MCP   = ifelse(a < gamma * lam1, -1 / gamma, 0),
      SCAD  = ifelse(a > lam1 & a < gamma * lam1, -1 / (gamma - 1), 0))
  } else {
    d1 <- switch(penalty,
      LASSO = lam1,
      MCP   = pmax(lam1 - a / gamma, 0),
      SCAD  = ifelse(a <= lam1, lam1, pmax(gamma * lam1 - a, 0) / (gamma - 1)))
    d1 / a
  }
  pp + lam2
}

## tr{ (G + diag(curv))^{-1} G }, or NA when G + diag(curv) is not positive definite.
.df_trace <- function(G, curv) {
  if (!length(curv)) return(0)
  M <- G
  diag(M) <- diag(M) + curv
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(ch)) return(NA_real_)
  sum(diag(backsolve(ch, forwardsolve(t(ch), G))))
}

## df along a BRIERi path. std.X: the standardized design; b: p x L standardized
## coefficients; b0: the L intercepts; w: the observation weights (sum to one). The
## intercept is unpenalized and always counted.
.df_path_i <- function(std.X, b, b0, lambda, penalty.factor, alpha, gamma, penalty,
                       family, w, method) {
  L <- ncol(b)
  out <- numeric(L); fallback <- 0L
  for (l in seq_len(L)) {
    A <- which(b[, l] != 0)
    k <- length(A) + 1
    X1 <- cbind(1, std.X[, A, drop = FALSE])
    vw <- if (identical(family, "gaussian")) {
      w
    } else {
      mu <- ginv_link(drop(b0[l] + std.X[, A, drop = FALSE] %*% b[A, l]), family)
      w * if (identical(family, "binomial")) mu * (1 - mu) else mu
    }
    G <- crossprod(X1, vw * X1)
    curv <- c(0, .df_curvature(b[A, l], lambda[l] * penalty.factor[A] * alpha,
                               lambda[l] * penalty.factor[A] * (1 - alpha),
                               penalty, gamma, method))
    d <- .df_trace(G, curv)
    if (!is.finite(d)) { d <- k; fallback <- fallback + 1L }
    out[l] <- d
  }
  list(df = out, fallback = fallback)
}

## df along a BRIERs path. XtX: the LD; b: p x L coefficients (standardized). No
## intercept on a summary target.
.df_path_s <- function(XtX, b, lambda, penalty.factor, alpha, gamma, penalty, eps,
                       method) {
  L <- ncol(b)
  out <- numeric(L); fallback <- 0L
  for (l in seq_len(L)) {
    A <- which(abs(b[, l]) >= eps)
    G <- as.matrix(XtX[A, A, drop = FALSE])
    curv <- .df_curvature(b[A, l], lambda[l] * penalty.factor[A] * alpha,
                          lambda[l] * penalty.factor[A] * (1 - alpha),
                          penalty, gamma, method)
    d <- .df_trace(G, curv)
    if (!is.finite(d)) { d <- length(A); fallback <- fallback + 1L }
    out[l] <- d
  }
  list(df = out, fallback = fallback)
}
