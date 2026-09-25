# =============================================================================
# screen.R: screening a panel of external models against a target.
#
# THE CASCADE, which is the application's cascade and the order matters:
#
#   non-degenerate
#     -> EXACTLY ONE screen, never both:
#          a validation split exists   -> k-fold CV screen, relative to the
#                                         target-only model, both legs at `bar`
#          no validation split         -> validation-free screen, from the
#                                         summary statistics and the LD alone
#       -> rank by THAT screen's own statistic
#         -> de-duplicate at `dedup.cor`, IN RANK ORDER
#           -> optionally keep the best K
#
# The two screens are EXCLUSIVE ALTERNATIVES. Chaining them, so that a model
# must clear both, is a different and stricter rule: measured on the ten-trait
# application it moves pooled sensitivity 0.834 -> 0.884 and specificity
# 0.496 -> 0.418, and leaves three of ten traits with one survivor or none. It
# is not what the application measured and it is not what this function does.
#
# DEDUP FOLLOWS RANKING, and that is not a detail. dedupExternals keeps the
# FIRST occurrence of each near-identical set, so which copy of a duplicate pair
# survives is decided by the order it is handed. Ranking first makes that the
# better-scoring copy; ranking afterwards would make it arbitrary.
# =============================================================================


#' Screen a panel of external models against a target
#'
#' Runs the screening cascade: remove degenerate models, apply exactly one
#' screen, rank the survivors, de-duplicate them in rank order and optionally
#' keep the best K.
#'
#' Which screen runs is decided by what is supplied, and the two are mutually
#' exclusive. Supplying \code{fit}, \code{X.val} and \code{y.val} selects the
#' CROSS-VALIDATED screen; supplying \code{corr}, \code{XtX} and \code{n}
#' selects the VALIDATION-FREE screen. Supplying both is an error rather than a
#' silent precedence rule, because the two answer different questions and a
#' caller who supplies both has not decided which one they want.
#'
#' \subsection{The cross-validated screen}{
#' Each external model is judged RELATIVE to the target-only model on the same
#' held-out samples, rather than against the null. The validation split is cut
#' into \code{nfolds} folds; in each fold the two remaining folds choose
#' \code{lambda} for the target-only path and fit the external's free intercept,
#' and everything is scored on the held-out fold. No refitting happens: the
#' coefficient path comes from \code{fit}, and the folds only re-select on it.
#'
#' Both legs must pass, and each relative bar is CLIPPED BY THE NULL MODEL:
#'
#' \deqn{\mathrm{loss}_{ext} \le \min(\mathrm{bar}^{-1}\,\mathrm{loss}_{target},\ \mathrm{loss}_{null})}
#' \deqn{\mathrm{acc}_{ext} \ge \max(\mathrm{bar}\,\mathrm{acc}_{target},\ \mathrm{acc}_{floor})}
#'
#' The clips are load bearing in both directions. Without the loss cap the bar
#' goes slack exactly when the target is weak, so a target at MSPE 1.0 would
#' admit scores predicting worse than the mean. Without the accuracy floor the
#' bar falls below chance whenever \code{bar < 1}, since the expected R-squared
#' of a signal-free predictor is \eqn{1/n}. Reported as the two ratios
#' \code{ratio_loss} and \code{ratio_acc}, both oriented so that a value at or
#' above 1 means the external beat its reference, and both averaged over folds
#' PAIRED WITHIN FOLD so that fold difficulty cancels.
#'
#' The accuracy floor uses the FOLD's n, not the split's. Each metric is
#' measured on roughly \code{1/nfolds} of the samples, and averaging folds cuts
#' the variance of the estimate but not its null bias, so the floor stays where
#' the measurement happened.
#' }
#'
#' \subsection{The validation-free screen}{
#' Computed from the coefficient vector, the target's marginal correlations and
#' the LD alone, so it needs no held-out samples and cannot leak any. With
#' \eqn{b} the external coefficients, \eqn{r} the marginal correlations and
#' \eqn{R} the LD, both legs must pass:
#'
#' \deqn{d^2 = 2 b'r - b'Rb > 0 \qquad\text{and}\qquad R^2 = (b'r)^2 / (b'Rb) > \mathrm{floor}}
#'
#' \eqn{d^2} is \eqn{1 - \mathrm{MSPE}} for a standardised outcome, so the first
#' leg asks that the score beat the mean at its published scale. The second is
#' the squared correlation the score would achieve after the best rescaling, so
#' it is scale free, and its floor is the no-signal level \eqn{1/n}, with
#' \code{n} the GWAS sample size of the target. The floor is \eqn{1/n} for EVERY
#' family, binary outcomes included: \eqn{R^2} here is a squared correlation on
#' the correlation scale BRIERs works in, not an AUC, so the 0.5 floor that belongs
#' to an AUC does not apply to it.
#' }
#'
#' @param beta.external A p x M matrix of external coefficients, no intercept
#'   row. A (p+1) x M matrix with an intercept row is accepted and the row is
#'   dropped.
#' @param fit A \code{BRIER} fit of the target at \code{eta = 0}, supplying the
#'   coefficient path the folds re-select on. Required by the CV screen.
#' @param X.val,y.val The validation split, on the RAW scale. Standardisation
#'   moments are taken from each fold's training folds and applied to the
#'   held-out fold, so a pre-standardised split would let the held-out samples
#'   set their own scale.
#' @param covariates An optional n x q matrix of covariates to residualise
#'   \code{y.val} against, refitted SEPARATELY within each piece of each fold so
#'   neither piece sees the other. Gaussian only, ignored otherwise. Omit it if
#'   \code{y.val} is already residualised, accepting that a single global
#'   residual has let every fold see every other.
#' @param corr A length-p vector of the target's marginal correlations.
#'   Required by the validation-free screen.
#' @param XtX A p x p LD matrix, sparse or dense. Required by the
#'   validation-free screen.
#' @param n The target's GWAS sample size, setting the no-signal floor.
#'   Required by the validation-free screen.
#' @param family One of \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}.
#' @param bar The relative bar for the CV screen, applied to both legs.
#'   Defaults to 1, meaning the external must match the target-only model.
#' @param nfolds Number of cross-validation folds. Defaults to 3.
#' @param seed Seed for the fold assignment. \code{NULL} leaves the RNG alone.
#' @param dedup.cor Absolute cosine similarity at which a survivor counts as a
#'   duplicate of one already kept. \code{NULL} or \code{NA} disables the step.
#' @param top.k Keep only the best K survivors after de-duplication.
#'   \code{NULL} keeps them all.
#' @param degenerate.tol A model whose coefficients are all zero, or whose PRS
#'   variance falls below this, carries no signal and is removed before any
#'   screen runs.
#' @param labels Optional names for the M external models.
#' @param trace Logical. If TRUE, reports each fold as it runs.
#'
#' @return A list of class \code{"BRIER.screen"}:
#' \describe{
#'   \item{externals}{A data frame with one row per external model: its index
#'     and label, its nonzero count, whether it was degenerate, the screen's own
#'     statistics, \code{keep} (passed the screen), \code{dedup} (survived
#'     de-duplication and the top-K cut) and \code{rank} (position among the
#'     final set, 1 is best, NA for a model that did not reach it).}
#'   \item{screen}{Which screen ran, \code{"cv"} or \code{"validation-free"}.}
#'   \item{survivors}{Indices passing the screen, before de-duplication.}
#'   \item{kept}{Indices of the final set, IN RANK ORDER.}
#'   \item{rank.statistic}{The column \code{kept} was ordered on.}
#'   \item{n.degenerate,n.screened.out,n.deduplicated}{What each stage removed.}
#'   \item{bar,nfolds,floor}{The settings the screen actually used.}
#' }
#'
#' @seealso \code{\link{dedupExternals}}, \code{\link{BRIERi}},
#'   \code{\link{BRIERs}}
#'
#' @examples
#' \dontrun{
#' # validation-free: summary statistics and LD only
#' scr <- screenExternals(beta.external = B, corr = r, XtX = R, n = 150000,
#'                        family = "gaussian", top.k = 10)
#' scr$externals[scr$externals$dedup, ]
#'
#' # cross-validated against a target-only fit
#' base <- BRIERs(sumstats, XtX, beta.external = B, eta.list = 0)
#' scr <- screenExternals(beta.external = B, fit = base,
#'                        X.val = X.val, y.val = y.val, family = "gaussian")
#' }
#'
#' @export
screenExternals <- function(
  beta.external,
  fit = NULL,
  X.val = NULL,
  y.val = NULL,
  covariates = NULL,
  corr = NULL,
  XtX = NULL,
  n = NULL,
  family = c("gaussian", "binomial", "poisson"),
  bar = 1,
  nfolds = 3L,
  seed = NULL,
  dedup.cor = 0.9,
  top.k = NULL,
  degenerate.tol = 1e-8,
  labels = NULL,
  trace = FALSE
) {
  family <- match.arg(family)
  B <- as.matrix(beta.external)
  M <- ncol(B)
  if (M < 1L) { stop("The external panel is empty: nothing to screen.", call. = FALSE) }
  if (is.null(labels)) {
    labels <- colnames(B)
    if (is.null(labels)) { labels <- paste0("external", seq_len(M)) }
  }
  if (length(labels) != M) {
    stop(sprintf("labels has length %d but the panel has %d models.", length(labels), M),
         call. = FALSE)
  }

  ## WHICH SCREEN. Exclusive by construction: supplying the inputs for both is
  ## refused rather than resolved by a precedence rule the caller cannot see.
  has.cv <- !is.null(fit) && !is.null(X.val) && !is.null(y.val)
  has.vf <- !is.null(corr) && !is.null(XtX)
  if (has.cv && has.vf) {
    stop(paste0("Both screens were supplied. They are alternatives, not stages: ",
                "pass fit/X.val/y.val for the cross-validated screen, or ",
                "corr/XtX/n for the validation-free screen, not both."),
         call. = FALSE)
  }
  if (!has.cv && !has.vf) {
    stop(paste0("No screen can run. Supply either fit, X.val and y.val (the ",
                "cross-validated screen) or corr, XtX and n (the ",
                "validation-free screen)."), call. = FALSE)
  }

  st <- if (has.cv) {
    .screen_cv(B, fit, X.val, y.val, covariates, family, bar, nfolds, seed,
               degenerate.tol, trace)
  } else {
    .screen_valfree(B, corr, XtX, n, family, degenerate.tol)
  }

  ## STAGE 2, ranking. It removes nothing, it only ORDERS, and each screen ranks
  ## on its OWN statistic: the CV screen on the loss ratio it was judged by, the
  ## validation-free screen on the scale-free R-squared. Ranking one screen by
  ## the other's statistic would order the survivors by a quantity that did not
  ## decide their admission.
  score <- st$stats[[st$rank.statistic]]
  keep <- st$stats$keep
  idx <- which(keep)
  if (length(idx)) {
    v <- score[idx]; v[!is.finite(v)] <- -Inf
    idx <- idx[order(-v, labels[idx])]
  }
  survivors <- idx

  ## STAGE 3, de-duplication IN RANK ORDER, so the surviving copy of a duplicate
  ## pair is the better-scoring one rather than whichever came first in the file.
  n.dedup <- 0L
  if (length(idx) >= 2L && !is.null(dedup.cor) && !is.na(dedup.cor)) {
    dd <- dedupExternals(B[, idx, drop = FALSE], cor.min = dedup.cor,
                         labels = labels[idx])
    n.dedup <- length(idx) - length(dd$keep)
    idx <- idx[dd$keep]
  }

  ## STAGE 4, the optional top-K cut. A study-design choice, not a property of
  ## the screen, which is why it is off by default.
  if (!is.null(top.k) && is.finite(top.k) && top.k > 0L && length(idx) > top.k) {
    idx <- idx[seq_len(as.integer(top.k))]
  }

  out.stats <- st$stats
  out.stats$index <- seq_len(M)
  out.stats$label <- labels
  out.stats$dedup <- FALSE
  out.stats$dedup[idx] <- TRUE
  out.stats$rank <- NA_integer_
  out.stats$rank[idx] <- seq_along(idx)
  lead <- c("index", "label", "nonzero", "degenerate")
  out.stats <- out.stats[, c(lead, setdiff(names(out.stats), lead)), drop = FALSE]

  structure(list(
    externals      = out.stats,
    screen         = st$screen,
    survivors      = survivors,
    kept           = idx,
    rank.statistic = st$rank.statistic,
    n.degenerate   = sum(out.stats$degenerate),
    n.screened.out = sum(!out.stats$degenerate & !out.stats$keep),
    n.deduplicated = n.dedup,
    bar            = if (has.cv) bar else NA_real_,
    nfolds         = if (has.cv) as.integer(nfolds) else NA_integer_,
    floor          = st$floor,
    family         = family
  ), class = "BRIER.screen")
}


#' @keywords internal
.acc_floor <- function(family, n) {
  if (identical(family, "binomial")) 0.5 else 1 / max(1L, n)
}

#' @keywords internal
.screen_metrics <- function(family) {
  ## The pair the screen judges on: a LOSS, lower better, and an ACCURACY,
  ## higher better. Poisson follows the convention the MCP screen already used
  ## (deviance and squared Spearman); it is the one combination in this function
  ## that no application result has measured.
  switch(family,
    gaussian = c(loss = "gaussian.mspe", acc = "gaussian.rsq"),
    binomial = c(loss = "binomial.dev",  acc = "binomial.auc"),
    poisson  = c(loss = "poisson.dev",   acc = "spearman.rsq")
  )
}

#' @keywords internal
.acc_value <- function(pred, y, family) {
  if (identical(family, "poisson")) {
    s <- suppressWarnings(stats::cor(pred, y, method = "spearman"))
    return(if (is.finite(s)) s^2 else NA_real_)
  }
  m <- .screen_metrics(family)[["acc"]]
  evalMetric(pred, y, m)
}


#' Validation-free screening statistics
#'
#' The leak-free legs, computed from the coefficients, the target's marginal
#' correlations and the LD alone.
#'
#' @keywords internal
.screen_valfree <- function(B, corr, XtX, n, family, tol) {
  p <- nrow(XtX)
  if (nrow(B) == p + 1L) { B <- B[-1, , drop = FALSE] }
  if (nrow(B) != p) {
    stop(sprintf("beta.external has %d rows but XtX is %d x %d.", nrow(B), p, p),
         call. = FALSE)
  }
  corr <- as.numeric(corr)
  if (length(corr) != p) {
    stop(sprintf("corr has length %d but XtX is %d x %d.", length(corr), p, p),
         call. = FALSE)
  }
  if (is.null(n) || !is.finite(n) || n <= 0) {
    stop(paste0("The validation-free screen needs the target's GWAS sample size `n`: ",
                "it sets the no-signal floor that the discrimination leg is judged ",
                "against."), call. = FALSE)
  }

  btr <- as.numeric(crossprod(B, corr))            # b'r
  bRb <- colSums(as.matrix(B * (XtX %*% B)))       # b'Rb
  nz <- colSums(B != 0)
  degenerate <- nz == 0L | !is.finite(bRb) | bRb < tol

  d2 <- 2 * btr - bRb                              # 1 - MSPE at the published scale
  r2 <- ifelse(is.finite(bRb) & bRb > 0, btr^2 / bRb, NA_real_)
  # 1/n for EVERY family. r2 is a squared correlation on the correlation scale
  # BRIERs works in, whatever the outcome's family: a binary trait enters through
  # how corr was derived, not through what r2 measures. The 0.5 of .acc_floor() is
  # an AUC floor and belongs to the CV screen, which measures a real AUC on held-out
  # samples. v1.4.0 applied 0.5 here, and since a PGS r2 is typically 0.01 to 0.1 it
  # dropped every external on a binary summary target (cad, t2d, ckd: 0 kept, where
  # the application's 1/N floor keeps 15, 84 and 11).
  fl <- 1 / n
  keep <- !degenerate &
    is.finite(d2) & d2 > 0 &
    is.finite(r2) & r2 > fl

  list(
    screen = "validation-free",
    rank.statistic = "r2",
    floor = fl,
    stats = data.frame(
      nonzero = nz, degenerate = degenerate,
      inner_product = btr, prs_variance = bRb,
      d2 = d2, r2 = r2, keep = keep,
      stringsAsFactors = FALSE
    )
  )
}


#' Cross-validated screening statistics
#'
#' Every external measured against the target-only model on the same held-out
#' samples, paired within fold.
#'
#' @keywords internal
.screen_cv <- function(B, fit, X.val, y.val, covariates, family, bar, nfolds,
                       seed, tol, trace) {
  X.val <- as.matrix(X.val)
  y.val <- as.numeric(y.val)
  N <- nrow(X.val)
  if (length(y.val) != N) {
    stop(sprintf("X.val has %d rows but y.val has length %d.", N, length(y.val)),
         call. = FALSE)
  }
  if (nrow(B) == ncol(X.val) + 1L) { B <- B[-1, , drop = FALSE] }
  if (nrow(B) != ncol(X.val)) {
    stop(sprintf("beta.external has %d rows but X.val has %d columns.",
                 nrow(B), ncol(X.val)), call. = FALSE)
  }
  nfolds <- as.integer(nfolds)
  if (!is.finite(nfolds) || nfolds < 2L) {
    stop("nfolds must be at least 2: one fold leaves nothing held out.", call. = FALSE)
  }
  if (nfolds > N) {
    stop(sprintf("nfolds is %d but the validation split has %d samples.", nfolds, N),
         call. = FALSE)
  }
  ## BRIERs fits carry XtY and no response; BRIERi fits carry the response. The
  ## two have different selection entry points and the same class.
  summary.module <- !is.null(fit$XtY)
  M <- ncol(B)
  bin <- identical(family, "binomial")
  mm <- .screen_metrics(family)
  sel.crit <- mm[["loss"]]

  if (!is.null(seed)) { set.seed(seed) }
  fold <- sample(rep_len(seq_len(nfolds), N))

  ## Residualise y WITHIN each piece of each fold, refitting separately, so
  ## neither piece has seen the other. A single global residual computed once on
  ## the whole split would let every fold see every other.
  within.resid <- function(i) {
    yi <- y.val[i]
    if (bin || is.null(covariates)) { return(yi) }
    Ci <- as.matrix(covariates)[i, , drop = FALSE]
    as.numeric(stats::resid(stats::lm(yi ~ Ci)))
  }

  ratio.loss <- matrix(NA_real_, M, nfolds)
  ratio.acc <- matrix(NA_real_, M, nfolds)
  for (kf in seq_len(nfolds)) {
    hold <- which(fold == kf); trn <- which(fold != kf)
    n.h <- length(hold)
    y.trn <- within.resid(trn); y.hold <- within.resid(hold)

    ## Standardise on the TRAINING folds and carry those constants into the
    ## held-out fold. X is a predictor: the model must meet the held-out samples
    ## on the scale it was trained on.
    sx <- standardize_X(X.val[trn, , drop = FALSE])
    ctr <- sx$center; scl <- sx$scale
    bad <- !is.finite(scl) | scl <= 0
    scl[bad] <- 1
    Xt <- sx$standardized
    Xh <- sweep(sweep(X.val[hold, , drop = FALSE], 2L, ctr, "-"), 2L, scl, "/")
    if (any(bad)) { Xt[, bad] <- 0; Xh[, bad] <- 0 }

    ## THE TARGET REFERENCE. lambda is chosen on the training folds and scored on
    ## the held-out fold. Nothing is refitted: the path came from `fit`.
    sel <- if (summary.module) {
      BRIERs.selection(fit, criteria = sel.crit, X.val = Xt, y.val = y.trn)
    } else {
      BRIERi.selection(fit, criteria = sel.crit, X.val = Xt, y.val = y.trn)
    }
    pt <- drop(stats::predict(sel, X = Xh, which.eta = 1L,
                              which.lambda = sel$lambda.min.index, type = "response"))
    t.loss <- evalMetric(pt, y.hold, mm[["loss"]])
    t.acc <- .acc_value(pt, y.hold, family)

    ## THE NULL MODEL on these same samples, measured rather than assumed, so the
    ## clip sits on the same scale as the thing it clips.
    null.loss <- if (bin) evalMetric(rep(0.5, n.h), y.hold, mm[["loss"]]) else 1 - 1 / n.h
    acc.fl <- .acc_floor(family, n.h)
    loss.ref <- min(t.loss, null.loss)
    acc.ref <- max(t.acc, acc.fl)

    ## EVERY EXTERNAL on the same held-out samples. The free intercept is fitted
    ## OFF-FOLD and applied to the held-out fold, so the external gets exactly the
    ## treatment the target got.
    Pt <- Xt %*% B
    Ph <- Xh %*% B
    for (m in seq_len(M)) {
      pm <- if (bin) stats::plogis(Ph[, m]) else mean(y.trn - Pt[, m]) + Ph[, m]
      if (all(is.na(pm))) { next }
      e.loss <- evalMetric(pm, y.hold, mm[["loss"]])
      e.acc <- .acc_value(pm, y.hold, family)
      ratio.loss[m, kf] <- loss.ref / e.loss
      ratio.acc[m, kf] <- e.acc / acc.ref
    }
    if (isTRUE(trace)) {
      message(sprintf("fold %d of %d: %d held out, target %s %.6f, null %.6f",
                      kf, nfolds, n.h, mm[["loss"]], t.loss, null.loss))
    }
  }

  ## Averaged over folds, PAIRED WITHIN FOLD: each ratio is formed inside its own
  ## fold before averaging, so fold difficulty cancels out of the comparison
  ## rather than being carried into it.
  rl <- rowMeans(ratio.loss, na.rm = TRUE)
  ra <- rowMeans(ratio.acc, na.rm = TRUE)
  rl[is.nan(rl)] <- NA_real_; ra[is.nan(ra)] <- NA_real_
  nz <- colSums(B != 0)
  degenerate <- nz == 0L
  keep <- !degenerate &
    is.finite(rl) & rl >= bar &
    is.finite(ra) & ra >= bar

  list(
    screen = "cv",
    rank.statistic = "ratio_loss",
    floor = NA_real_,
    stats = data.frame(
      nonzero = nz, degenerate = degenerate,
      ratio_loss = rl, ratio_acc = ra,
      sd_ratio_loss = apply(ratio.loss, 1L, stats::sd, na.rm = TRUE),
      sd_ratio_acc = apply(ratio.acc, 1L, stats::sd, na.rm = TRUE),
      keep = keep,
      stringsAsFactors = FALSE
    )
  )
}


#' @rdname screenExternals
#' @param x A \code{BRIER.screen} object.
#' @param ... Ignored.
#' @export
print.BRIER.screen <- function(x, ...) {
  cat(sprintf("BRIER screen: %s, %s outcome\n", x$screen, x$family))
  n <- nrow(x$externals)
  cat(sprintf("  %d candidates -> %d non-degenerate -> %d survivors -> %d after dedup",
              n, n - x$n.degenerate, length(x$survivors), sum(x$externals$dedup)))
  if (length(x$kept) < length(x$survivors) - x$n.deduplicated) {
    cat(" and the top-K cut")
  }
  cat("\n")
  if (identical(x$screen, "cv")) {
    cat(sprintf("  %d folds, both legs at bar %.4g, ranked on %s\n",
                x$nfolds, x$bar, x$rank.statistic))
  } else {
    cat(sprintf("  no-signal floor %.3g, ranked on %s\n", x$floor, x$rank.statistic))
  }
  invisible(x)
}
