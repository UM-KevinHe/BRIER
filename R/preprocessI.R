#' Preprocess inputs for BRIERi genotype-based modeling
#'
#' `preprocessI()` prepares SNP-level input information for `BRIERi()` when
#' individual-level genotype data are available in the target cohort and one or
#' more external models are available as SNP-level coefficients. The function
#' performs genotype annotation quality control, harmonizes SNP representation
#' across data sources, and aligns external model coefficients to the SNP order
#' and allele orientation of the target genotype matrix.
#'
#' The argument `target.info` should be a data.frame describing the SNPs in the
#' target genotype matrix, with rows ordered to match the columns of the target
#' predictor matrix. The optional argument `external.ss` should be a data.frame
#' containing SNP annotation columns and one or more external coefficient
#' columns, specified through `external.coef.cols`.
#'
#' First, `preprocessI()` performs quality control and harmonization separately
#' for the target SNP information and the external coefficient file. SNP
#' annotation columns are standardized to `CHR`, `BP`, `REF`, and `ALT`, where
#' `CHR` is chromosome, `BP` is base-pair position, `REF` is the reference
#' allele, and `ALT` is the alternative/effect allele. Chromosome labels are
#' normalized by removing `"chr"` prefixes, converting PLINK chromosome codes
#' `23`, `24`, `25`, and `26` to `X`, `Y`, `XY`, and `MT`, and converting `M` to `MT`.
#' SNPs with unrecognized chromosome labels are removed with a warning. The
#' function also removes sites with duplicated CHR:BP entries 
#' (treated as multi-allelic or input duplicates), and optionally removes 
#' strand-ambiguous SNPs, including A/T and C/G allele pairs, when 
#' `drop.ambiguous = TRUE`.
#'
#' After quality control, the target cohort is used as the canonical reference
#' for SNP order and allele orientation. SNP identifiers in the processed output
#' are written as `CHR:BP:REF:ALT`, using the target-cohort allele orientation.
#' External SNPs are matched to target SNPs by `CHR:BP`. If the external alleles
#' match the target orientation, the external coefficient is retained unchanged.
#' If the external `REF` and `ALT` alleles are reversed relative to the target,
#' the external coefficient sign is flipped. If the alleles conflict, the
#' corresponding external coefficient is set to zero. SNPs present in the target
#' cohort but absent from an external model are retained and assigned an
#' external coefficient of zero, whereas SNPs present only in the external model
#' are not retained because the output is aligned to the target genotype matrix.
#'
#' External models whose aligned coefficient columns are entirely zero are
#' removed from the output. If all supplied external models are removed after
#' alignment, the function returns `external.ss = NULL` with a warning.
#' 
#' The function errors out if any required annotation or coefficient column
#' contains missing values.
#' 
#' The returned `target.info.keep` element gives the indices of SNPs retained
#' from the original target SNP information and should be used to subset the
#' columns of the target genotype matrix before calling `BRIERi()`. The returned
#' `external.ss` element contains external coefficients aligned to the retained
#' target SNP order and allele orientation.
#'
#' `BRIERi()` expects the external coefficient matrix to include an intercept
#' slot as the first row. If an external-model intercept is unavailable or is not
#' used, this intercept slot should be set to zero. The intercept row is required
#' for `BRIERi()` but is not needed for `BRIERs()`.
#'
#' @param target.info A data.frame describing the columns of the target
#'   genotype matrix in order. Must contain chromosome, base-pair position, REF, ALT.
#' @param external.ss Optional data.frame of external model coefficients.
#'   Must contain chromosome, base-pair position, REF, ALT, plus one or more
#'   coefficient columns.
#' @param target.info.cols Named character vector mapping the canonical keys
#'   \code{chr}, \code{bp}, \code{ref}, \code{alt} to column names in
#'   \code{target.info}. Defaults to
#'   \code{c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT")}; override
#'   any key whose column name differs in your input.
#' @param external.ss.cols Named character vector mapping canonical keys
#'   (\code{chr}, \code{bp}, \code{ref}, \code{alt}) to the column names in
#'   \code{external.ss}. Defaults to
#'   \code{c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT")}; 
#'   override any key whose column name differs in your input.
#' @param external.coef.cols Character vector of column names in
#'   \code{external.ss} that hold model coefficients. Required when
#'   \code{external.ss} is supplied.
#' @param drop.ambiguous Logical. If TRUE (default), drops strand-ambiguous
#'   SNPs (A/T and C/G allele pairs) from every input before matching.
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
#' @seealso \code{\link{BRIERi}}, \code{\link{mergeExternals}},
#'   \code{\link{preprocessS}}
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
#'   out$external.ss[, grep("^coef", colnames(out$external.ss)), drop = FALSE]
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
        "(b) chromosome or base-pair position coding differs between inputs ",
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