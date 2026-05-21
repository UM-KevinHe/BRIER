#' Preprocess inputs for BRIERi (genotype as predictors)
#'
#' Aligns external model coefficients to the SNP order of a target
#' individual-level genotype matrix for use with \code{\link{BRIERi}}.
#' \code{target.info} describes the columns of the target X matrix (in
#' order) and defines the canonical allele orientation: external
#' coefficients are matched by chromosome and position, then sign-flipped
#' where their alleles are swapped relative to \code{target.info}.
#'
#' \code{BRIERi} estimates an intercept internally from the target data,
#' so external model intercepts are not required as input. The
#' \code{beta.external} matrix passed to \code{BRIERi} must include an
#' intercept slot as the first row (typically zero, since no external
#' intercept is being imposed); see the example below.
#'
#' Inputs are validated for missing values: if \code{NA} is detected in
#' any of the required columns of \code{target.info} or \code{external.ss}
#' (chr, bp, ref, alt, or coefficient columns), the function stops with
#' an error. Clean the inputs before calling.
#'
#' Chromosome coding is normalized internally: \code{"chr"} prefix is
#' stripped, PLINK numeric codes (23/24/26) mapped to \code{"X"/"Y"/"MT"},
#' and \code{"M"} normalized to \code{"MT"}. SNPs with unrecognized CHR
#' values (anything not in \code{"1"}\eqn{-}\code{"22"}, \code{"X"},
#' \code{"Y"}, \code{"MT"}, \code{"XY"}) are dropped with a warning.
#' \code{varnames} in the output uses the canonical form.
#'
#' External models whose coefficient column is entirely zero after
#' alignment (no overlap with \code{target.info}) are dropped from the
#' output. If \emph{all} supplied external models are dropped, the
#' function emits a warning and returns \code{external.ss = NULL}.
#'
#' Multi-allelic variants (duplicate \code{CHR:BP} entries) are dropped
#' with a warning from every input independently. Strand-ambiguous SNPs
#' (alleles A/T or C/G) are dropped by default; set
#' \code{drop.ambiguous = FALSE} to keep them when strand is reliably
#' matched across studies.
#'
#' All three relevant data.frames (\code{target.info}, \code{external.ss})
#' share a \code{varnames} column of the form \code{"CHR:BP:REF:ALT"} (in
#' the target.info orientation) as the first column.
#'
#' @param target.info A data.frame describing the columns of the target
#'   genotype matrix in order. Must contain chromosome, basepair, REF, ALT.
#' @param external.ss Optional data.frame of external model coefficients.
#'   Must contain chromosome, basepair, REF, ALT, plus one or more
#'   coefficient columns.
#' @param target.info.cols Named character vector mapping \code{chr},
#'   \code{bp}, \code{ref}, \code{alt} to column names in
#'   \code{target.info}.
#' @param external.ss.cols Named character vector mapping \code{chr},
#'   \code{bp}, \code{ref}, \code{alt} to column names in
#'   \code{external.ss}.
#' @param external.coef.cols Character vector of column names in
#'   \code{external.ss} that hold model coefficients. Required when
#'   \code{external.ss} is supplied.
#' @param drop.ambiguous Logical. If TRUE (default), drops strand-ambiguous
#'   SNPs (A/T, T/A, C/G, G/C alleles) from every input before matching.
#' @param verbose Logical. If TRUE (default), prints SNP counts at each
#'   filter step.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{target.info}{Processed target SNP info as a data.frame with
#'     columns \code{varnames}, \code{CHR}, \code{BP}, \code{REF},
#'     \code{ALT}. Rows correspond to the surviving columns of the target
#'     X matrix.}
#'   \item{target.info.keep}{Integer vector of indices into the original
#'     \code{target.info} for the surviving SNPs. Use as
#'     \code{X[, target.info.keep]} to subset the target genotype matrix
#'     to the analysis set.}
#'   \item{external.ss}{Processed external coefficients (or \code{NULL} if
#'     no external input was supplied, or if every supplied external
#'     model had all-zero coefficients after alignment), with columns
#'     \code{varnames}, \code{CHR}, \code{BP}, \code{REF}, \code{ALT},
#'     \code{coef1}, \code{coef2}, ... where surviving coefficient columns
#'     are renumbered \code{coef1}, ..., \code{coefK} in input order.}
#'   \item{external.coef.names}{Named character vector mapping the
#'     standardized column names (\code{coef1}, \code{coef2}, ...) back to
#'     the original column names supplied in \code{external.coef.cols}, or
#'     \code{NULL} if no external input survived.}
#'   \item{summary}{Named integer vector of SNP counts:
#'     \code{n.input.tinfo}, \code{n.input.ess},
#'     \code{n.dropped.invalid.chr}, \code{n.dropped.ambiguous},
#'     \code{n.dropped.multiallelic}, \code{n.final},
#'     \code{n.external.matched}, \code{n.external.flipped},
#'     \code{n.external.conflict}, \code{n.external.zero}.}
#' }
#'
#' @seealso \code{\link{BRIERi}}, \code{\link{preprocessS}},
#'   \code{\link{mergeExternals}}
#'
#' @examples
#' \dontrun{
#' # target.info parallels the columns of X
#' out <- preprocessI(
#'   target.info        = target_info,
#'   external.ss        = prs_coefs,
#'   external.coef.cols = c("PRS_EUR", "PRS_EAS")
#' )
#'
#' # Subset X to the surviving SNPs
#' X.aligned <- X[, out$target.info.keep]
#'
#' # Extract beta.external matrix and prepend a zero row for the intercept
#' # slot expected by BRIERi (BRIERi estimates the intercept internally).
#' beta.ext <- as.matrix(
#'   out$external.ss[, grep("^coef", colnames(out$external.ss))]
#' )
#' beta.ext <- rbind(0, beta.ext)
#'
#' fit <- BRIERi(
#'   X             = X.aligned,
#'   y             = y,
#'   beta.external = beta.ext,
#'   ...
#' )
#' }
#'
#' @export
preprocessI <- function(
  target.info,
  external.ss = NULL,

  target.info.cols = c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT"),
  external.ss.cols = c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT"),
  external.coef.cols = NULL,

  drop.ambiguous = TRUE,
  verbose = TRUE
) {

  # ============================================================
  # 1. Coerce to data.frame and standardize column names
  # ============================================================
  target.info <- as.data.frame(target.info)
  if (!is.null(external.ss)) external.ss <- as.data.frame(external.ss)

  target.info <- .standardize_cols(
    target.info, target.info.cols,
    c("chr", "bp", "ref", "alt"),
    "target.info"
  )

  if (!is.null(external.ss)) {
    if (is.null(external.coef.cols) || length(external.coef.cols) == 0) {
      stop(
        "external.coef.cols must be supplied when external.ss is provided. ",
        "Specify the column names that hold model coefficients, e.g. ",
        "external.coef.cols = c('PRS1', 'PRS2').",
        call. = FALSE
      )
    }
    external.ss <- .standardize_cols(
      external.ss, external.ss.cols,
      c("chr", "bp", "ref", "alt"),
      "external.ss"
    )
    missing.coef <- setdiff(external.coef.cols, colnames(external.ss))
    if (length(missing.coef) > 0) {
      stop(
        "external.coef.cols not found in external.ss: ",
        paste(missing.coef, collapse = ", "), ".", call. = FALSE
      )
    }
  }

  # ============================================================
  # 2. NA check on required columns
  # ============================================================
  .check_na(target.info, c("CHR", "BP", "REF", "ALT"), "target.info")
  if (!is.null(external.ss)) {
    .check_na(
      external.ss,
      c("CHR", "BP", "REF", "ALT", external.coef.cols),
      "external.ss"
    )
  }

  # ============================================================
  # 3. Normalize allele case and CHR coding; drop invalid CHR
  # ============================================================
  target.info$REF <- toupper(target.info$REF)
  target.info$ALT <- toupper(target.info$ALT)
  if (!is.null(external.ss)) {
    external.ss$REF <- toupper(external.ss$REF)
    external.ss$ALT <- toupper(external.ss$ALT)
  }

  # Capture original-row indices BEFORE any row drops, so target.info.keep
  # always references the user-supplied target.info.
  target.info$.orig.idx <- seq_len(nrow(target.info))

  n.input.tinfo <- nrow(target.info)
  n.input.ess   <- if (is.null(external.ss)) 0L else nrow(external.ss)

  if (verbose) {
    message(
      "Input SNPs: ",
      "target.info=", n.input.tinfo,
      ", external.ss=", n.input.ess
    )
  }

  # Normalize CHR (strip "chr" prefix, map 23/24/26 -> X/Y/MT, "M" -> "MT").
  # Drop rows with unrecognized CHR values.
  pre.tinfo.chr <- nrow(target.info)
  target.info$CHR <- .normalize_chr(target.info$CHR)
  target.info     <- .drop_invalid_chr(target.info, "target.info", verbose)
  n.dropped.tinfo.chr <- pre.tinfo.chr - nrow(target.info)

  n.dropped.ess.chr <- 0L
  if (!is.null(external.ss)) {
    pre.ess.chr <- nrow(external.ss)
    external.ss$CHR <- .normalize_chr(external.ss$CHR)
    external.ss     <- .drop_invalid_chr(external.ss, "external.ss", verbose)
    n.dropped.ess.chr <- pre.ess.chr - nrow(external.ss)
  }

  n.dropped.invalid.chr <- n.dropped.tinfo.chr + n.dropped.ess.chr

  # ============================================================
  # 4. Drop strand-ambiguous (optional)
  # ============================================================
  n.dropped.ambiguous <- 0L
  if (isTRUE(drop.ambiguous)) {
    pre.tinfo <- nrow(target.info)
    pre.ess   <- if (is.null(external.ss)) 0L else nrow(external.ss)
    target.info <- .drop_ambiguous(target.info, "target.info", verbose)
    if (!is.null(external.ss)) {
      external.ss <- .drop_ambiguous(external.ss, "external.ss", verbose)
    }
    n.dropped.ambiguous <-
      (pre.tinfo - nrow(target.info)) +
      (pre.ess - (if (is.null(external.ss)) 0L else nrow(external.ss)))
  }

  # ============================================================
  # 5. Drop multi-allelic (duplicate CHR:BP)
  # ============================================================
  pre.tinfo <- nrow(target.info)
  pre.ess   <- if (is.null(external.ss)) 0L else nrow(external.ss)

  target.info <- .drop_multiallelic(target.info, "target.info", verbose)
  if (!is.null(external.ss)) {
    external.ss <- .drop_multiallelic(external.ss, "external.ss", verbose)
  }
  n.dropped.multiallelic <-
    (pre.tinfo - nrow(target.info)) +
    (pre.ess - (if (is.null(external.ss)) 0L else nrow(external.ss)))

  # ============================================================
  # 6. target.info is canonical — finalize and build varnames
  # ============================================================
  target.info.keep <- target.info$.orig.idx
  target.info$.orig.idx <- NULL

  varnames <- paste(
    target.info$CHR, target.info$BP,
    target.info$REF, target.info$ALT,
    sep = ":"
    )
  target.info$varnames <- varnames
  target.info <- target.info[, c("varnames", "CHR", "BP", "REF", "ALT"), drop = FALSE]
  rownames(target.info) <- NULL

  n.final <- nrow(target.info)

  # ============================================================
  # 7. Align external.ss to target.info
  # ============================================================
  n.external.matched  <- 0L
  n.external.flipped  <- 0L
  n.external.conflict <- 0L
  n.external.zero     <- 0L
  external.ss.out     <- NULL
  external.coef.names <- NULL

  if (!is.null(external.ss)) {
    key.tinfo <- paste(target.info$CHR, target.info$BP, sep = ":")
    key.ess   <- paste(external.ss$CHR, external.ss$BP, sep = ":")

    me <- match(key.tinfo, key.ess)
    tinfo.alleles <- paste(target.info$REF, target.info$ALT, sep = "/")
    ess.alleles   <- paste(external.ss$REF[me], external.ss$ALT[me], sep = "/")
    ess.flipped   <- paste(external.ss$ALT[me], external.ss$REF[me], sep = "/")

    e.match    <- !is.na(me) & (tinfo.alleles == ess.alleles)
    e.flip     <- !is.na(me) & (tinfo.alleles == ess.flipped) & !e.match
    e.conflict <- !is.na(me) & !e.match & !e.flip
    e.zero     <- is.na(me)

    M <- length(external.coef.cols)
    coef.mat <- matrix(0, nrow = n.final, ncol = M)

    if (any(e.match)) {
      coef.mat[e.match, ] <- as.matrix(
        external.ss[me[e.match], external.coef.cols, drop = FALSE]
      )
    }
    if (any(e.flip)) {
      coef.mat[e.flip, ] <- -as.matrix(
        external.ss[me[e.flip], external.coef.cols, drop = FALSE]
      )
    }
    # e.conflict and e.zero rows remain 0

    # Drop external models whose post-alignment coefficients are entirely 0.
    all.zero <- colSums(coef.mat != 0) == 0L

    if (all(all.zero)) {
      warning(
        "All external models have only zero coefficients after alignment ",
        "with target.info. This usually indicates: ",
        "(a) target.info and external.ss do not overlap on CHR:BP; ",
        "(b) chromosome or basepair coding differs between inputs ",
        "(e.g., hg19 vs hg38 build); or ",
        "(c) allele coding is inconsistent. ",
        "Returning external.ss = NULL. Consider rerunning preprocessI without ",
        "external models, or check the inputs.",
        call. = FALSE
      )
      external.ss.out     <- NULL
      external.coef.names <- NULL
    } else {
      if (any(all.zero)) {
        dropped.orig <- external.coef.cols[all.zero]
        warning(
          "Dropping external model(s) with all-zero coefficients after ",
          "alignment: ", paste(dropped.orig, collapse = ", "), ".",
          call. = FALSE
        )
      }
      keep <- !all.zero
      kept.coef.mat <- coef.mat[, keep, drop = FALSE]
      K <- ncol(kept.coef.mat)
      new.std.names <- paste0("coef", seq_len(K))
      colnames(kept.coef.mat) <- new.std.names
      external.coef.names <- setNames(external.coef.cols[keep], new.std.names)

      external.ss.out <- data.frame(
        varnames = varnames,
        CHR = target.info$CHR,
        BP  = target.info$BP,
        REF = target.info$REF,
        ALT = target.info$ALT,
        kept.coef.mat,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    n.external.matched  <- sum(e.match)
    n.external.flipped  <- sum(e.flip)
    n.external.conflict <- sum(e.conflict)
    n.external.zero     <- sum(e.zero)

    if (verbose) {
      message(
        "External: ", n.external.matched, " matched, ",
        n.external.flipped, " flipped, ",
        n.external.conflict, " conflict, ",
        n.external.zero, " zero-padded."
      )
    }
  }

  if (verbose) message("Final analysis set: ", n.final, " SNPs.")

  # ============================================================
  # 8. Return
  # ============================================================
  list(
    target.info         = target.info,
    target.info.keep    = target.info.keep,
    external.ss         = external.ss.out,
    external.coef.names = external.coef.names,
    summary             = c(
      n.input.tinfo          = n.input.tinfo,
      n.input.ess            = n.input.ess,
      n.dropped.invalid.chr  = n.dropped.invalid.chr,
      n.dropped.ambiguous    = n.dropped.ambiguous,
      n.dropped.multiallelic = n.dropped.multiallelic,
      n.final                = n.final,
      n.external.matched     = n.external.matched,
      n.external.flipped     = n.external.flipped,
      n.external.conflict    = n.external.conflict,
      n.external.zero        = n.external.zero
    )
  )
}