# =============================================================================
# Preparation of the external panel, shared by BRIER.I and BRIER.S.
#
# Two steps act on beta.external before the aggregation in calcExtY() and
# calcExtXtY() sees it:
#
#   1. DE-DUPLICATION, which every multi-source fit runs first. A panel of
#      published external models routinely carries near-copies of itself. A
#      repeated source cannot add information, and it makes the stacking Gram
#      singular, so the second copy is dropped and the first occurrence kept.
#   2. PRINCIPAL-COMPONENT REDUCTION, requested with multi.method = "PCstacking".
#      The survivors still overlap heavily, so the panel is projected onto its
#      leading directions and those are combined instead of the raw models.
#
# Both work on the coefficient vectors alone. Neither touches y, X, or the LD
# matrix, so neither can leak held-out information into the fit.
# =============================================================================


# Residual norm below which a column counts as spanned by the columns already
# kept. Strict, because a false drop is worse than a missed one: on unit-norm
# columns this is the norm left after projecting out the kept set.
.DEDUP_RANK_TOL <- 1e-8


# A readable name per external column for the drop report. Positional indices
# alone are useless when the caller passed a hundred catalogue models.
.external_labels <- function(labels, B) {
  M <- ncol(B)
  if (!is.null(labels) && length(labels) == M) { return(as.character(labels)) }
  cn <- colnames(B)
  if (!is.null(cn) && length(cn) == M && all(nzchar(cn))) { return(cn) }
  paste("external", seq_len(M))
}


#' Drop redundant external models from a panel
#'
#' Removes external models that duplicate one already kept, keeping the first
#' occurrence of each near-identical set. Called automatically by
#' \code{\link{BRIERi}}, \code{\link{BRIERi.cv}}, \code{\link{BRIERs}},
#' \code{\link{BRIERi.bopt}} and \code{\link{BRIERs.bopt}} whenever more than
#' one external model is supplied, and exported so the same reduction can be
#' inspected or applied on its own.
#'
#' Two kinds of redundancy are removed. First, any model whose absolute cosine
#' similarity with an already-kept model reaches \code{cor.min}. Second, any
#' model that the kept set spans exactly, which pairwise similarity cannot see:
#' a model equal to the sum of two others is uncorrelated with either. Columns
#' that are entirely zero are also dropped, since they transfer nothing and put
#' a zero on the diagonal of the stacking Gram matrix.
#'
#' The test is on the coefficient vectors, not on the polygenic scores they
#' produce. If one coefficient vector is a multiple of another then their scores
#' agree in every metric, whatever the linkage disequilibrium is, so a drop here
#' is provable from the coefficients alone and needs no LD matrix. The converse
#' does not hold: two weakly correlated coefficient vectors can still give
#' collinear scores through LD. The rule is deliberately conservative in that
#' direction, and leaves LD-induced collinearity to the fitter.
#'
#' Similarity is the uncentred correlation. Centring would be wrong here: these
#' vectors are mostly zeros by construction, since a predictor an external model
#' does not cover enters as a zero coefficient, and the quantity stacking
#' inverts is an uncentred inner product.
#'
#' @param beta.external A numeric matrix of external coefficients, p x M for the
#'   summary-level convention or (p+1) x M with a leading intercept row for the
#'   individual-level one.
#' @param cor.min Similarity threshold in (0, 1]. A model is dropped when its
#'   absolute cosine similarity with an already-kept model reaches this value.
#'   \code{NULL} or \code{NA} disables de-duplication entirely.
#' @param intercept.row Logical. \code{TRUE} when the first row of
#'   \code{beta.external} is an intercept rather than a predictor coefficient,
#'   as \code{\link{BRIERi}} requires. The intercept is held out of the
#'   similarity computation and carried along with its column.
#' @param labels Optional character vector of length M naming the external
#'   models, used in the report. Defaults to the column names, then to
#'   positional labels.
#'
#' @return A list with elements:
#' \describe{
#'   \item{beta.external}{The panel with redundant columns removed.}
#'   \item{keep}{Integer indices of the surviving columns in the input.}
#'   \item{dropped}{A data.frame of the removed columns with the reason for each,
#'     or \code{NULL} when nothing was removed.}
#'   \item{M.in, M.out}{Column counts before and after.}
#'   \item{cor.min}{The threshold used.}
#'   \item{applied}{Logical, whether any column was removed.}
#'   \item{note}{Character report lines, empty when nothing was removed.}
#' }
#'
#' @seealso \code{\link{reduceExternalsPCA}}, \code{\link{calcExtY}},
#'   \code{\link{calcExtXtY}}
#'
#' @examples
#' \dontrun{
#' b1 <- c(rep(1, 5), rep(0, 45))
#' beta.external <- cbind(b1, b1 * 2, rnorm(50))
#' ded <- dedupExternals(beta.external, cor.min = 0.9)
#' ded$keep
#' ded$dropped
#' }
#'
#' @export
dedupExternals <- function(
  beta.external,
  cor.min = 0.9,
  intercept.row = FALSE,
  labels = NULL
) {

  B <- as.matrix(beta.external)
  storage.mode(B) <- "double"

  icpt <- NULL
  if (isTRUE(intercept.row)) {
    if (nrow(B) < 2L) {
      stop("beta.external has no predictor rows below the intercept.", call. = FALSE)
    }
    icpt <- B[1, , drop = FALSE]
    B <- B[-1, , drop = FALSE]
  }
  M <- ncol(B)

  out <- list(
    beta.external = as.matrix(beta.external),
    keep          = seq_len(M),
    dropped       = NULL,
    M.in          = M,
    M.out         = M,
    cor.min       = cor.min,
    applied       = FALSE,
    note          = character(0)
  )

  ## Disabled, or nothing to be redundant with. A single external has no partner
  ## to duplicate, so the whole question only exists from M = 2.
  if (is.null(cor.min) || length(cor.min) != 1L || is.na(cor.min)) { return(out) }
  if (!is.numeric(cor.min) || cor.min <= 0 || cor.min > 1) {
    stop("cor.min must be NULL, NA, or a similarity threshold in (0, 1].", call. = FALSE)
  }
  if (M < 2L) { return(out) }

  lab <- .external_labels(labels, B)
  reason <- rep(NA_character_, M)
  dup.of <- rep(NA_character_, M)

  ## An all-zero external is singular on its own, before any pair is compared.
  nrm <- sqrt(colSums(B * B))
  zero <- !is.finite(nrm) | nrm <= 0
  reason[zero] <- "all zero, so it transfers nothing"

  ## Pairwise near-duplicates. The inner loop compares only against columns that
  ## are still kept, which is what makes the rule transitive and deterministic:
  ## a column is only ever called a duplicate of a survivor, and the survivor is
  ## always the earliest member of its set.
  live <- which(!zero)
  Bn <- NULL
  if (length(live) >= 2L) {
    Bn <- sweep(B[, live, drop = FALSE], 2L, nrm[live], "/")
    G <- crossprod(Bn)
    for (a in seq_along(live)) {
      if (a < 2L) { next }
      ja <- live[a]
      for (b in seq_len(a - 1L)) {
        jb <- live[b]
        if (!is.na(reason[jb])) { next }
        s <- abs(G[b, a])
        if (is.finite(s) && s >= cor.min) {
          reason[ja] <- sprintf(
            "duplicate of %s (cosine similarity %.5f >= %.5f)", lab[jb], s, cor.min
          )
          dup.of[ja] <- lab[jb]
          break
        }
      }
    }
  }

  ## Backstop for exact spans, which pairwise similarity cannot detect. Columns
  ## are swept left to right so the earliest occurrence is again the one kept;
  ## a pivoted QR would find the same rank but choose its basis by column norm
  ## rather than by position.
  keep <- which(is.na(reason))
  if (length(keep) >= 2L && !is.null(Bn)) {
    idx <- match(keep, live)
    basis <- NULL
    dep <- integer(0)
    for (i in seq_along(keep)) {
      v <- Bn[, idx[i]]
      if (!is.null(basis)) {
        v <- v - basis %*% crossprod(basis, v)
        v <- v - basis %*% crossprod(basis, v)  # re-orthogonalise once
      }
      nv <- sqrt(sum(v^2))
      if (nv <= .DEDUP_RANK_TOL) {
        dep <- c(dep, keep[i])
      } else {
        basis <- cbind(basis, as.numeric(v) / nv)
      }
    }
    if (length(dep)) {
      reason[dep] <- sprintf(
        "linearly dependent on the %d external model(s) kept before it",
        length(keep) - length(dep)
      )
      keep <- which(is.na(reason))
    }
  }

  ## Every column redundant means every column is all zero. Returning a
  ## zero-column panel would break the fitter for no gain, so say so and change
  ## nothing: the fit is then equivalent to eta = 0.
  if (!length(keep)) {
    out$note <- sprintf(paste0(
      "All %d external model(s) are all zero, so none of them transfers ",
      "anything. The panel was left as-is and the fit is equivalent to eta = 0."
    ), M)
    return(out)
  }
  if (length(keep) == M) { return(out) }

  dropped <- setdiff(seq_len(M), keep)
  B <- B[, keep, drop = FALSE]
  if (!is.null(icpt)) { B <- rbind(icpt[, keep, drop = FALSE], B) }

  out$beta.external <- B
  out$keep    <- keep
  out$M.out   <- length(keep)
  out$applied <- TRUE
  out$dropped <- data.frame(
    index        = dropped,
    label        = lab[dropped],
    reason       = reason[dropped],
    duplicate.of = dup.of[dropped],
    stringsAsFactors = FALSE
  )
  out$note <- c(
    sprintf(paste0(
      "dedupExternals: dropped %d of %d external model(s) as redundant, %d ",
      "carried into the fit. A repeated source adds no information and makes ",
      "the stacking weights unidentifiable."
    ), length(dropped), M, length(keep)),
    sprintf("  dropped column %d (%s): %s", dropped, lab[dropped], reason[dropped])
  )
  out
}


#' Reduce an external panel to its leading principal directions
#'
#' Projects a panel of external coefficient vectors onto the leading directions
#' of its singular value decomposition, keeping components until \code{pca.var}
#' of the variance is covered. This is the first half of
#' \code{multi.method = "PCstacking"}, whose second half combines the retained
#' components by stacking.
#'
#' Screened external panels stay highly redundant with one another even after
#' exact duplicates are removed, so combining all of them means solving a system
#' dominated by noise in the directions that carry least. Projecting first and
#' combining the projections leaves the low-variance directions out of the
#' combination entirely.
#'
#' The decomposition is uncentred, deliberately. Centring across sources would
#' subtract a mean model, which is not a meaningful reference here: each column
#' is a model rather than a draw from a population, and a zero coefficient means
#' the variant has no effect, not that it is below average.
#'
#' What is returned is \code{B \%*\% V[, 1:k]}, the projection of the original
#' models onto the leading right singular directions, so the result stays on the
#' coefficient scale that eta shrinks toward rather than on an orthonormal basis
#' whose scale would be arbitrary relative to eta.
#'
#' @param beta.external A numeric matrix of external coefficients, p x M or
#'   (p+1) x M with a leading intercept row.
#' @param pca.var Fraction of the variance to retain, in (0, 1]. Components are
#'   added in order until the cumulative share reaches this value.
#' @param intercept.row Logical. \code{TRUE} when the first row is an intercept.
#'   That row is not a predictor direction: including it in the decomposition
#'   would let an intercept dominate a component, and the fitter would then
#'   reject the result on its (p+1) shape check. It is held out, the predictor
#'   rows are reduced, and the same combination of the input intercepts is put
#'   back on top.
#'
#' @return A list with elements:
#' \describe{
#'   \item{beta.external}{The reduced panel, p x k (or (p+1) x k), columns named
#'     \code{PC1 ... PCk}.}
#'   \item{n.pcs}{The number of components retained.}
#'   \item{variance.explained, cumulative}{Per-component and cumulative variance
#'     shares over all M components.}
#'   \item{rotation}{The M x k matrix of right singular vectors used, with the
#'     source labels as row names, so the reduction can be read as a set of
#'     weights on the input models.}
#'   \item{applied}{Logical, whether a reduction was performed.}
#' }
#'
#' @seealso \code{\link{dedupExternals}}, \code{\link{calcExtY}},
#'   \code{\link{calcExtXtY}}
#'
#' @examples
#' \dontrun{
#' red <- reduceExternalsPCA(beta.external, pca.var = 0.8)
#' red$n.pcs
#' round(red$cumulative, 3)
#' }
#'
#' @export
reduceExternalsPCA <- function(
  beta.external,
  pca.var = 0.8,
  intercept.row = FALSE
) {

  B <- as.matrix(beta.external)
  storage.mode(B) <- "double"

  icpt <- NULL
  if (isTRUE(intercept.row)) {
    if (nrow(B) < 2L) {
      stop("beta.external has no predictor rows below the intercept.", call. = FALSE)
    }
    icpt <- B[1, , drop = FALSE]
    B <- B[-1, , drop = FALSE]
  }
  M <- ncol(B)

  if (!is.numeric(pca.var) || length(pca.var) != 1L || !is.finite(pca.var) ||
      pca.var <= 0 || pca.var > 1) {
    stop("pca.var must be a single number in (0, 1].", call. = FALSE)
  }

  untouched <- list(
    beta.external      = as.matrix(beta.external),
    n.pcs              = M,
    variance.explained = rep(1 / max(M, 1L), M),
    cumulative         = seq_len(M) / max(M, 1L),
    rotation           = NULL,
    applied            = FALSE
  )

  ## One source has no redundancy to remove, and the leading direction of a
  ## single vector is the vector itself.
  if (M < 2L) { return(untouched) }

  sv <- svd(B)
  tot <- sum(sv$d^2)
  if (!is.finite(tot) || tot <= 0) { return(untouched) }

  ve  <- sv$d^2 / tot
  cum <- cumsum(ve)
  k <- which(cum >= pca.var)[1]
  if (is.na(k)) { k <- M }

  V <- sv$v[, seq_len(k), drop = FALSE]
  rownames(V) <- .external_labels(NULL, B)
  colnames(V) <- paste0("PC", seq_len(k))

  out <- B %*% V
  colnames(out) <- paste0("PC", seq_len(k))

  ## A component is a combination of models, so its intercept is that same
  ## combination of the inputs' intercepts. Computed rather than assumed, so a
  ## non-zero external intercept carries through correctly.
  if (!is.null(icpt)) { out <- rbind(icpt %*% V, out) }

  list(
    beta.external      = out,
    n.pcs              = k,
    variance.explained = ve,
    cumulative         = cum,
    rotation           = V,
    applied            = TRUE
  )
}


# -- Internal: the de-duplication step the fitters run, plus its report --
#
# One place so all five entry points behave identically and phrase the drop the
# same way. De-duplication changes what is fitted, so it is never silent.

.dedup_and_report <- function(beta.external, cor.min, intercept.row) {
  ded <- dedupExternals(beta.external, cor.min = cor.min, intercept.row = intercept.row)
  if (length(ded$note)) { message(paste(ded$note, collapse = "\n")) }
  ded
}


# Per-source settings (eta grids, search bounds) are sized to the panel the
# caller supplied. When de-duplication shortens the panel, they follow it, so a
# dropped duplicate does not turn into a length-mismatch error the caller has no
# way to anticipate. Settings already sized to the aggregated panel, as they are
# under PCA, stacking and PCstacking, are left alone.
.dedup_follow <- function(x, ded) {
  if (is.null(x) || !ded$applied) { return(x) }
  if (length(x) != ded$M.in) { return(x) }
  x[ded$keep]
}
