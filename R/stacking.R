# =============================================================================
# Stacking weights: one solver for every shape and family.
#
# Both entry points used to solve NORMAL EQUATIONS directly:
#
#   BRIERs   w = solve(crossprod(B, XtX %*% B), crossprod(B, XtY))   calcExtXtY
#   BRIERi   w = solve(crossprod(Y), crossprod(Y, z))                stacking_gaussian
#
# Normal equations square the condition number. External panels are collinear by
# construction, so even after deduplication the Gram can be near singular and the
# weights then blow up. Measured on the ten-trait height application: the raw Gram
# reached rcond 3.5e-11 on one trait and 7e-16 on another, while the equilibrated
# Gram stayed above 1e-3 in every set.
#
# The fix, and it is the same fix that application used:
#
#   1. EQUILIBRATE. Scale the Gram to unit diagonal, D = diag(G)^(-1/2), and solve
#      on D G D and D h. Jacobi preconditioning: the same solution in exact
#      arithmetic, a far better conditioned system in floating point.
#   2. SOLVE ON THE EIGENDECOMPOSITION, inverting only directions above a relative
#      tolerance. This reduces to the ordinary solve when the Gram is well
#      conditioned and returns a finite, minimum-norm answer when it is not,
#      instead of a wild one.
#   3. REPORT what happened: the condition number before and after equilibration,
#      and how many directions were dropped. A caller should be able to see that
#      truncation occurred rather than infer it from strange weights.
#
# No ridge for now (Ray, 2026-09-24): truncation is parameter free, and a ridge
# would change the estimator rather than the arithmetic.
#
# The constrained variant solves the same quadratic under the constraint
# w >= 0 and sum(w) = 1, by projected gradient with Armijo backtracking and a
# Barzilai-Borwein trial step, stopped on the KKT residual rather than on the
# objective. This is the routine that produced the application's constrained
# results, ported unchanged so the package and the manuscript agree.
# =============================================================================

#' Projection onto the constrained weight set
#'
#' The constrained set is w >= 0 with sum(w) = 1 (the probability simplex). This
#' is its Euclidean projection, by Duchi et al. (2008): sort descending, find the
#' largest index whose shifted value stays positive, subtract the implied
#' threshold and clamp at zero.
#'
#' @param v A numeric vector.
#' @return A numeric vector of the same length, non-negative and summing to one.
#' @keywords internal
.project_weights <- function(v) {
  K <- length(v)
  if (K == 1L) { return(1) }
  u <- sort(v, decreasing = TRUE)
  cs <- cumsum(u)
  ok <- u - (cs - 1) / seq_len(K) > 0
  rho <- if (any(ok)) max(which(ok)) else 1L
  theta <- (cs[rho] - 1) / rho
  w <- pmax(v - theta, 0)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) { rep(1 / K, K) } else { w / s }
}

#' KKT residual of the constrained problem
#'
#' For min f(w) subject to sum(w) = 1 and w >= 0 the conditions are: some lambda
#' with g_k = lambda wherever w_k > 0, and g_k >= lambda wherever w_k = 0. The
#' residual below is zero exactly at an optimum and is scale free, so it can
#' certify the answer where a stalled objective cannot.
#'
#' @param w A numeric weight vector on the simplex.
#' @param g The gradient at \code{w}.
#' @param active_tol Weights above this are treated as active.
#' @return A non-negative scalar, or NA when no weight is active.
#' @keywords internal
.kkt_residual <- function(w, g, active_tol = 1e-10) {
  on <- w > active_tol
  if (!any(on)) { return(NA_real_) }
  lam <- mean(g[on])
  spread_on <- if (sum(on) > 1L) max(g[on]) - min(g[on]) else 0
  short_off <- if (any(!on)) max(0, lam - min(g[!on])) else 0
  scale <- max(1e-12, max(abs(g)))
  max(spread_on, short_off) / scale
}

#' Projected gradient for the constrained problem
#'
#' Armijo backtracking with a Barzilai-Borwein trial step, stopped on the KKT
#' residual. Convex objective assumed.
#'
#' @param fn Objective, taking a weight vector.
#' @param gr Gradient, taking a weight vector.
#' @param K Number of weights.
#' @param w_init Optional starting point; defaults to equal weights.
#' @param maxit Maximum outer iterations.
#' @param tol Relative objective tolerance for the stall test.
#' @param patience Stop after this many consecutive iterations without a relative
#'   objective improvement above \code{tol}. The KKT residual can plateau just
#'   above \code{kkt_tol} while the weights no longer move, and without this the
#'   solver runs to \code{maxit} for no change in the answer.
#' @param kkt_tol Stop once the KKT residual falls below this.
#' @param step0 Initial step length.
#' @return A list with \code{w}, \code{value}, \code{iters} and \code{kkt}.
#' @keywords internal
.constrained_solve <- function(fn, gr, K, w_init = NULL, maxit = 20000L,
                           tol = 1e-15, kkt_tol = 1e-8, step0 = 1,
                           patience = 50L) {
  w <- if (is.null(w_init)) rep(1 / K, K) else .project_weights(w_init)
  f <- fn(w)
  g <- gr(w)
  step <- step0
  iters <- 0L
  stalled <- 0L
  for (it in seq_len(maxit)) {
    iters <- it
    moved <- FALSE
    for (ls in seq_len(80L)) {
      w_try <- .project_weights(w - step * g)
      d <- w_try - w
      f_try <- fn(w_try)
      if (is.finite(f_try) && f_try <= f + 1e-4 * sum(g * d)) {
        moved <- TRUE
        break
      }
      step <- step / 2
      if (step < 1e-18) { break }
    }
    if (!moved) { break }
    g_try <- gr(w_try)
    sv <- w_try - w
    yv <- g_try - g
    sy <- sum(sv * yv)
    step <- if (is.finite(sy) && sy > 0) {
      max(1e-14, min(1e14, sum(sv * sv) / sy))
    } else {
      step * 2
    }
    delta <- f - f_try
    w <- w_try
    f <- f_try
    g <- g_try
    kk <- .kkt_residual(w, g)
    if (is.finite(kk) && kk < kkt_tol) { break }
    stalled <- if (delta <= tol * max(1, abs(f))) stalled + 1L else 0L
    if (stalled >= patience) { break }
  }
  list(w = w, value = f, iters = iters, kkt = .kkt_residual(w, gr(w)))
}

#' Stacking weights from a Gram matrix and a cross-product vector
#'
#' Solves \eqn{\min_w w'Gw - 2w'h}, which is the stacking least-squares problem
#' in whatever metric the caller supplies: \eqn{G = B' X'X B} and \eqn{h = B'X'y}
#' for the summary-statistic shape, \eqn{G = Y'Y} and \eqn{h = Y'z} for the
#' individual-level shape. One routine so the two shapes cannot drift apart.
#'
#' The unconstrained solve equilibrates the Gram to unit diagonal and inverts
#' only the eigen-directions above a relative tolerance, so a near-singular panel
#' returns a finite minimum-norm answer rather than an arbitrary one. The
#' constrained solve minimises the same quadratic under the constraint.
#'
#' @param G A symmetric K x K Gram matrix.
#' @param h A numeric vector of length K.
#' @param constrain Logical. \code{FALSE} (default) for the ordinary solve,
#'   \code{TRUE} to constrain the weights to the probability simplex.
#' @param tol Relative eigenvalue tolerance: directions whose eigenvalue is below
#'   \code{tol} times the largest are dropped. Only used when
#'   \code{constrain = FALSE}.
#' @param kkt_tol Stopping tolerance on the KKT residual of the constrained
#'   solve. The default matches the height application, which settled on 1e-8:
#'   tighter values run to the iteration cap without changing the weights.
#'
#' @return A list with:
#' \describe{
#'   \item{weights}{The numeric weight vector, length K.}
#'   \item{constrained}{Whether the constraint was applied.}
#'   \item{rcond_raw}{Reciprocal condition number of \code{G} as supplied.}
#'   \item{rcond_equilibrated}{Reciprocal condition number after equilibration.}
#'   \item{n_dropped}{Eigen-directions dropped by the tolerance (0 when none).}
#'   \item{kkt}{KKT residual of the constrained solve, NA when unconstrained.}
#' }
#'
#' @keywords internal
stack_weights <- function(G, h, constrain = FALSE, tol = 1e-12, kkt_tol = 1e-8) {
  G <- as.matrix(G)
  h <- as.numeric(h)
  K <- ncol(G)
  if (K < 1L) { stop("The external panel is empty: nothing to stack.", call. = FALSE) }
  if (nrow(G) != K) { stop("G must be square.", call. = FALSE) }
  if (length(h) != K) { stop("h must have one entry per column of G.", call. = FALSE) }
  rc_raw <- tryCatch(rcond(G), error = function(e) NA_real_)

  if (constrain) {
    fn <- function(v) { as.numeric(crossprod(v, G %*% v)) - 2 * sum(h * v) }
    gr <- function(v) { 2 * (as.numeric(G %*% v) - h) }
    ft <- .constrained_solve(fn, gr, K, kkt_tol = kkt_tol)
    return(list(
      weights = ft$w, constrained = TRUE, rcond_raw = rc_raw,
      rcond_equilibrated = NA_real_, n_dropped = 0L, kkt = ft$kkt
    ))
  }

  # equilibrate: unit diagonal, so a fifty-fold spread in model size stops
  # dominating the conditioning
  dg <- diag(G)
  dg[!is.finite(dg) | dg <= 0] <- 1
  Dv <- 1 / sqrt(dg)
  Ge <- G * outer(Dv, Dv)
  Ge <- (Ge + t(Ge)) / 2                     # symmetrise away round-off
  rc_eq <- tryCatch(rcond(Ge), error = function(e) NA_real_)

  ev <- eigen(Ge, symmetric = TRUE)
  lam <- ev$values
  keep <- is.finite(lam) & lam > tol * max(lam, na.rm = TRUE)
  if (!any(keep)) {
    stop("The stacking Gram has no usable directions; the external panel is degenerate.",
         call. = FALSE)
  }
  he <- Dv * h
  w <- Dv * drop(ev$vectors[, keep, drop = FALSE] %*%
                   ((t(ev$vectors[, keep, drop = FALSE]) %*% he) / lam[keep]))
  list(
    weights = as.numeric(w), constrained = FALSE, rcond_raw = rc_raw,
    rcond_equilibrated = rc_eq, n_dropped = sum(!keep), kkt = NA_real_
  )
}

#' Constrained stacking weights from response-scale predictions
#'
#' The simplex-constrained counterpart of \code{\link{stacking_gaussian}},
#' \code{\link{stacking_binomial}} and \code{\link{stacking_poisson}}: it
#' minimises the same family loss under the constraint \eqn{w \ge 0},
#' \eqn{\sum_j w_j = 1},
#' where the stacked prediction is the convex combination
#' \eqn{\mu_i = \sum_j w_j y_{ij}} on the RESPONSE scale.
#'
#' Losses and their gradients, with \eqn{\partial Q/\partial w_k =
#' \sum_i (\partial Q/\partial \mu_i) y_{ik}}:
#' \describe{
#'   \item{gaussian}{\eqn{Q = \tfrac{1}{2}\sum_i (z_i - \mu_i)^2},
#'     \eqn{\partial Q/\partial \mu_i = \mu_i - z_i}}
#'   \item{binomial}{\eqn{Q = -\sum_i [z_i \log p_i + (1 - z_i)\log(1 - p_i)]},
#'     \eqn{\partial Q/\partial p_i = (p_i - z_i)/(p_i(1 - p_i))}}
#'   \item{poisson}{\eqn{Q = \sum_i [\mu_i - z_i \log \mu_i]},
#'     \eqn{\partial Q/\partial \mu_i = 1 - z_i/\mu_i}}
#' }
#'
#' Each loss is convex in \eqn{\mu} and \eqn{\mu} is linear in \eqn{w}, so the
#' problem is convex on the simplex and the projected-gradient solution is the
#' global optimum, certified by the KKT residual rather than by a stalled
#' objective.
#'
#' A softmax reparameterisation \eqn{w_k = e^{\alpha_k}/\sum_j e^{\alpha_j}}
#' reaches the same optimum and was used to derive these gradients independently;
#' it is not used here because it cannot put a weight exactly at zero and because
#' it carries a flat direction (adding a constant to every alpha leaves w
#' unchanged).
#'
#' @param Y An n x m matrix of external predictions on the response scale:
#'   probabilities for binomial, positive rates for poisson.
#' @param z The target response, length n.
#' @param family One of "gaussian", "binomial", "poisson".
#' @param w_init Optional starting weights; defaults to equal weights.
#' @param kkt_tol Stopping tolerance on the KKT residual.
#' @param eps Clamp keeping binomial probabilities inside \code{[eps, 1 - eps]}
#'   and poisson rates above \code{eps}, so the loss stays finite at the boundary.
#'
#' @return A list with \code{weights}, \code{constrained} (always \code{TRUE}),
#'   \code{fitted}, \code{value}, \code{kkt} and \code{iters}.
#'
#' @keywords internal
stack_weights_response <- function(Y, z, family = c("gaussian", "binomial", "poisson"),
                                   w_init = NULL, kkt_tol = 1e-8, eps = 1e-8) {
  family <- match.arg(family)
  Y <- as.matrix(Y)
  z <- as.numeric(z)
  m <- ncol(Y)
  if (nrow(Y) != length(z)) { stop("Y and z must have the same number of rows.", call. = FALSE) }
  if (m < 2L) {
    return(list(weights = 1, constrained = TRUE, fitted = drop(Y),
                value = NA_real_, kkt = NA_real_, iters = 0L))
  }

  mu_of <- function(w) { drop(Y %*% w) }
  if (family == "gaussian") {
    fn <- function(w) { r <- z - mu_of(w); 0.5 * sum(r * r) }
    gr <- function(w) { as.numeric(crossprod(Y, mu_of(w) - z)) }
  } else if (family == "binomial") {
    fn <- function(w) {
      p <- pmin(pmax(mu_of(w), eps), 1 - eps)
      -sum(z * log(p) + (1 - z) * log1p(-p))
    }
    gr <- function(w) {
      p <- pmin(pmax(mu_of(w), eps), 1 - eps)
      as.numeric(crossprod(Y, (p - z) / (p * (1 - p))))
    }
  } else {
    fn <- function(w) {
      mu <- pmax(mu_of(w), eps)
      sum(mu - z * log(mu))
    }
    gr <- function(w) {
      mu <- pmax(mu_of(w), eps)
      as.numeric(crossprod(Y, 1 - z / mu))
    }
  }

  ft <- .constrained_solve(fn, gr, m, w_init = w_init, kkt_tol = kkt_tol)
  ## `constrained` is carried on every report, including this one, so a caller
  ## can read $stack without knowing which solver produced it.
  list(weights = ft$w, constrained = TRUE, fitted = mu_of(ft$w),
       value = ft$value, kkt = ft$kkt, iters = ft$iters)
}

#' Constrained stacking weights by softmax reparameterisation of the constraint
#'
#' An independent route to the same optimum as \code{\link{stack_weights_response}},
#' kept as a cross-check rather than as the implementation. The constraint is
#' absorbed into the parameterisation
#' \eqn{w_k = e^{\alpha_k} / \sum_j e^{\alpha_j}}, so an ordinary unconstrained
#' optimiser can be used. With
#' \eqn{\partial w_k/\partial \alpha_j = w_k(\delta_{kj} - w_j)} the chain rule
#' gives \eqn{\partial \mu_i/\partial \alpha_j = w_j (y_{ij} - \mu_i)} and hence
#'
#' \deqn{\partial Q/\partial \alpha_j = w_j \sum_i (\partial Q/\partial \mu_i)(y_{ij} - \mu_i).}
#'
#' Two reasons this is the cross-check and not the implementation. The weights are
#' strictly positive, so a model that should be dropped only approaches zero as
#' \eqn{\alpha \to -\infty}, and the parameterisation is redundant: adding a
#' constant to every \eqn{\alpha} leaves \eqn{w} unchanged, leaving a flat
#' direction. The last coordinate is therefore pinned at zero as a reference.
#'
#' @inheritParams stack_weights_response
#' @param maxit Maximum iterations for \code{\link[stats]{optim}}.
#' @param reltol Relative convergence tolerance for \code{\link[stats]{optim}}.
#'
#' @return A list with \code{weights}, \code{fitted}, \code{value},
#'   \code{convergence} and \code{alpha}.
#'
#' @keywords internal
stack_weights_softmax <- function(Y, z, family = c("gaussian", "binomial", "poisson"),
                                  maxit = 5000L, reltol = 1e-12, eps = 1e-8) {
  family <- match.arg(family)
  Y <- as.matrix(Y)
  z <- as.numeric(z)
  m <- ncol(Y)
  if (m < 2L) { return(list(weights = 1, fitted = drop(Y), value = NA_real_, convergence = 0L, alpha = 0)) }

  # reference category: alpha_m = 0, removing the flat direction
  w_of <- function(a) { e <- exp(c(a, 0) - max(c(a, 0))); e / sum(e) }
  dQ_dmu <- function(mu) {
    if (family == "gaussian") { mu - z }
    else if (family == "binomial") { p <- pmin(pmax(mu, eps), 1 - eps); (p - z) / (p * (1 - p)) }
    else { m2 <- pmax(mu, eps); 1 - z / m2 }
  }
  obj <- function(a) {
    mu <- drop(Y %*% w_of(a))
    if (family == "gaussian") { r <- z - mu; 0.5 * sum(r * r) }
    else if (family == "binomial") { p <- pmin(pmax(mu, eps), 1 - eps); -sum(z * log(p) + (1 - z) * log1p(-p)) }
    else { m2 <- pmax(mu, eps); sum(m2 - z * log(m2)) }
  }
  grd <- function(a) {
    w <- w_of(a)
    mu <- drop(Y %*% w)
    d <- dQ_dmu(mu)
    # dQ/dalpha_j = w_j sum_i d_i (y_ij - mu_i), for j = 1..m-1
    full <- vapply(seq_len(m), function(j) w[j] * sum(d * (Y[, j] - mu)), numeric(1))
    full[seq_len(m - 1L)]
  }
  op <- stats::optim(
    par = rep(0, m - 1L),
    fn = obj,
    gr = grd,
    method = "BFGS",
    control = list(maxit = maxit, reltol = reltol)
  )
  w <- w_of(op$par)
  list(weights = w, fitted = drop(Y %*% w), value = op$value,
       convergence = op$convergence, alpha = c(op$par, 0))
}
