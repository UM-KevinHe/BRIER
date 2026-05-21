#' Preprocess summary statistics for BRIERs genotype-based modeling
#'
#' `preprocessS()` prepares SNP-level summary-statistic inputs for `BRIERs()`
#' when target-cohort individual-level genotypes are not available. The function
#' harmonizes target summary statistics, target LD reference-panel SNP
#' information, and optional external model coefficients, and then aligns all
#' inputs to the SNP order and allele orientation of the target LD reference
#' panel.
#'
#' The argument `target.ss` should be a data.frame containing target summary
#' statistics. The argument `target.ld` should be a data.frame describing the
#' SNPs in the target LD reference panel, with rows ordered to match the rows
#' and columns of the LD matrix used in downstream `BRIERs()` fitting. The
#' optional argument `external.ss` should be a data.frame containing SNP
#' annotation columns and one or more external coefficient columns, specified
#' through `external.coef.cols`.
#'
#' First, `preprocessS()` performs quality control and harmonization separately
#' for the target summary statistics, target LD reference panel, and external
#' coefficient file. SNP annotation columns are standardized to `CHR`, `BP`,
#' `REF`, and `ALT`, where `CHR` is chromosome, `BP` is base-pair position,
#' `REF` is the reference allele, and `ALT` is the alternative/effect allele.
#' Chromosome labels are normalized by removing `"chr"` prefixes, converting
#' PLINK chromosome codes `23`, `24`, `25`, and `26` to `X`, `Y`, `XY`, and `MT`, and
#' converting `M` to `MT`. SNPs with unrecognized chromosome labels are removed
#' with a warning. The function also removes sites with duplicated CHR:BP entries 
#' (treated as multi-allelic or input duplicates), and optionally removes 
#' strand-ambiguous SNPs, including A/T and C/G allele pairs, when 
#' `drop.ambiguous = TRUE`.
#'
#' When `target.ind = "gwas"`, `preprocessS()` converts GWAS summary
#' statistics into marginal SNP-outcome correlations. In this setting,
#' `target.ss` must contain p-values, sample sizes, and either effect signs or
#' beta coefficients. If beta coefficients are provided, effect directions are
#' inferred from the sign of beta. P-values equal to zero are replaced by the
#' smallest observed non-zero p-value before computing correlations. If all
#' p-values are zero, the function stops with an error. When
#' `target.ind = "corr"`, the user must directly provide a marginal
#' correlation column.
#'
#' After quality control, the target LD reference panel is used as the canonical
#' reference for SNP order and allele orientation. SNP identifiers in the
#' processed output are written as `CHR:BP:REF:ALT`, using the LD-reference
#' allele orientation. Target summary-statistic SNPs are matched to the LD
#' reference panel by `CHR:BP`. If the summary-statistic alleles match the LD
#' reference orientation, the marginal correlation is retained unchanged. If the
#' summary-statistic `REF` and `ALT` alleles are reversed relative to the LD
#' reference, the marginal correlation sign is flipped. SNPs without compatible
#' allele matches between `target.ss` and `target.ld` are removed from the final
#' analysis set.
#'
#' External model coefficients are then aligned to the surviving LD-reference
#' SNP set. External SNPs are matched by `CHR:BP`. If the external alleles match
#' the LD-reference orientation, the external coefficient is retained unchanged.
#' If the external `REF` and `ALT` alleles are reversed relative to the LD
#' reference, the external coefficient sign is flipped. If the alleles conflict,
#' the corresponding external coefficient is set to zero. SNPs present in the
#' final target/LD SNP set but absent from an external model are retained and
#' assigned an external coefficient of zero, whereas SNPs present only in the
#' external model are not retained because the output is aligned to the target
#' LD reference panel.
#'
#' External models whose aligned coefficient columns are entirely zero are
#' removed from the output. If all supplied external models are removed after
#' alignment, the function returns `external.ss = NULL` with a warning.
#' 
#' The function errors out if any required annotation, summary-statistic, 
#' LD-reference, or coefficient column contains missing values.
#'
#' The returned `target.ld.keep` element gives the indices of SNPs retained
#' from the original target LD reference panel and should be used to subset the
#' rows and columns of a precomputed LD matrix before calling `BRIERs()`. The
#' returned `target.ss` contains marginal correlations aligned to the retained
#' LD-reference SNP order and allele orientation. The returned `external.ss`
#' contains external coefficients aligned to the same SNP set.
#'
#' Unlike `BRIERi()`, `BRIERs()` does not require an intercept row in
#' `beta.external`.
#'
#' @param target.ss A data.frame of target summary statistics. Must contain
#'   chromosome, base-pair position, REF, ALT, plus value columns determined by
#'   \code{target.ind}.
#' @param target.ind One of \code{"gwas"} or \code{"corr"}. For
#'   \code{"gwas"}, \code{target.ss} must contain p-value, sample size, and
#'   either sign or beta columns. For \code{"corr"}, a correlation column is
#'   required.
#' @param target.ld A data.frame of the target LD reference panel SNP info,
#'   defining the canonical allele orientation. Must contain chromosome,
#'   base-pair position, REF, ALT. 
#' @param external.ss Optional data.frame of external model coefficients.
#'   Must contain chromosome, base-pair position, REF, ALT, plus one or more
#'   coefficient columns.
#' @param target.ss.cols Named character vector mapping canonical keys
#'   (\code{chr}, \code{bp}, \code{ref}, \code{alt}, \code{p}, \code{n},
#'   \code{sgn}, \code{beta}, \code{corr}) to the column names in
#'   \code{target.ss}. Only the keys relevant to \code{target.ind} need to
#'   be present. Defaults to
#'   \code{c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT", p = "pval",
#'   n = "n", sgn = "sgn", beta = "beta", corr = "corr")}; override any key
#'   whose column name differs in your input.
#' @param target.ld.cols Named character vector mapping canonical keys
#'   (\code{chr}, \code{bp}, \code{ref}, \code{alt}) to the column names in
#'   \code{target.ld}. Defaults to
#'   \code{c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT")}; 
#'   override any key whose column name differs in your input.
#' @param external.ss.cols Named character vector mapping canonical keys
#'   (\code{chr}, \code{bp}, \code{ref}, \code{alt}) to the column names in
#'   \code{external.ss}. Defaults to
#'   \code{c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT")}; 
#'   override any key whose column name differs in your input.
#' @param external.coef.cols Character vector of column names in
#'   \code{external.ss} that hold model coefficients. Required when
#'   \code{external.ss} is supplied; the function does not auto-detect.
#' @param drop.ambiguous Logical. If TRUE (default), drops strand-ambiguous
#'   SNPs (A/T and C/G allele pairs) from every input before matching.
#' @param verbose Logical. If TRUE (default), prints SNP counts at each
#'   filter step.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{target.ss}{Processed target sumstats as a data.frame with
#'     columns \code{varnames}, \code{CHR}, \code{BP}, \code{REF}, \code{ALT},
#'     \code{corr} (in the LD reference orientation).}
#'   \item{target.ld}{Processed target LD reference panel SNP info, subset
#'     to the surviving SNPs, with columns \code{varnames}, \code{CHR},
#'     \code{BP}, \code{REF}, \code{ALT}.}
#'   \item{target.ld.keep}{Integer vector of indices into the original
#'     \code{target.ld} for the surviving SNPs. Use as
#'     \code{LD[target.ld.keep, target.ld.keep]} to subset a pre-computed
#'     LD matrix.}
#'   \item{external.ss}{Processed external coefficients (or \code{NULL} if
#'     no external input was supplied, or if every supplied external model
#'     had all-zero coefficients after alignment), with columns
#'     \code{varnames}, \code{CHR}, \code{BP}, \code{REF}, \code{ALT},
#'     \code{coef1}, \code{coef2}, ... where surviving coefficient columns
#'     are renumbered \code{coef1}, ..., \code{coefK} in input order.}
#'   \item{external.coef.names}{Named character vector mapping the
#'     standardized column names (\code{coef1}, \code{coef2}, ...) back to
#'     the original column names supplied in \code{external.coef.cols}, or
#'     \code{NULL} if no external input survived.}
#'   \item{summary}{Named integer vector of SNP counts: \code{n.input.tss},
#'     \code{n.input.tld}, \code{n.input.ess}, \code{n.dropped.invalid.chr},
#'     \code{n.dropped.ambiguous}, \code{n.dropped.multiallelic},
#'     \code{n.dropped.unmatched.tss}, \code{n.aligned}, \code{n.flipped},
#'     \code{n.final}, \code{n.external.matched}, \code{n.external.flipped},
#'     \code{n.external.conflict}, \code{n.external.zero}.}
#' }
#'
#' @seealso \code{\link{BRIERs}}, \code{\link{calLD}}, \code{\link{p2cor}}, \code{\link{preprocessI}}, \code{\link{mergeExternals}}
#'
#' @examples
#' \dontrun{
#' out <- preprocessS(
#'   target.ss          = my_gwas,
#'   target.ind         = "gwas",
#'   target.ld          = hapmap3_panel,
#'   external.ss        = prs_coefs,
#'   external.coef.cols = c("PRS_EUR", "PRS_EAS"),
#'   target.ss.cols     = c(
#'     chr  = "Chromosome", bp = "Position",
#'     ref  = "Allele2",    alt = "Allele1",
#'     p    = "PVAL",       n  = "N",
#'     beta = "Effect"
#'   )
#' )
#'
#' LD.input <- precomputed_LD[out$target.ld.keep, out$target.ld.keep]
#' beta.ext <- as.matrix(
#'   out$external.ss[, grep("^coef", colnames(out$external.ss)), drop = FALSE]
#' )
#'
#' fit <- BRIERs(
#'   sumstats = out$target.ss, XtX = LD.input, beta.external = beta.ext, ...)
#' }
#'
#' @export
preprocessS <- function(
  target.ss,
  target.ind = c("gwas", "corr"),
  target.ld,
  external.ss = NULL,

  target.ss.cols = c(
    chr  = "CHR",  bp   = "BP",
    ref  = "REF",  alt  = "ALT",
    p    = "pval", n    = "n",
    sgn  = "sgn",  beta = "beta",
    corr = "corr"
  ),
  target.ld.cols = c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT"),
  external.ss.cols = c(chr = "CHR", bp = "BP", ref = "REF", alt = "ALT"),
  external.coef.cols = NULL,

  drop.ambiguous = TRUE,
  verbose = TRUE
) {

  target.ind <- match.arg(target.ind)

  # ============================================================
  # 1. Coerce to data.frame and standardize column names
  # ============================================================
  target.ss <- as.data.frame(target.ss)
  target.ld <- as.data.frame(target.ld)
  if (!is.null(external.ss)) external.ss <- as.data.frame(external.ss)

  # Required canonical keys for target.ss depend on input type
  req.target <- switch(
    target.ind,
    gwas = c("chr", "bp", "ref", "alt", "p", "n"),
    corr = c("chr", "bp", "ref", "alt", "corr")
  )

  target.ss <- .standardize_cols(target.ss, target.ss.cols, req.target, "target.ss")
  target.ld <- .standardize_cols(
    target.ld, target.ld.cols,
    c("chr", "bp", "ref", "alt"), "target.ld"
  )

  # GWAS: need sgn; derive from beta if not supplied. Then compute corr.
  if (target.ind == "gwas") {
    # Safely access optional keys — sgn and beta may be omitted from the mapping
    # (one is enough; sgn is derived from beta if only beta is supplied).
    sgn.col  <- if ("sgn"  %in% names(target.ss.cols)) target.ss.cols[["sgn"]]  else NA_character_
    beta.col <- if ("beta" %in% names(target.ss.cols)) target.ss.cols[["beta"]] else NA_character_

    if (!is.na(sgn.col) && sgn.col %in% colnames(target.ss)) {
      colnames(target.ss)[colnames(target.ss) == sgn.col] <- "sgn"
    } else if (!is.na(beta.col) && beta.col %in% colnames(target.ss)) {
      target.ss$sgn <- sign(target.ss[[beta.col]])
    } else {
      stop(
        "target.ind = 'gwas' requires either a 'sgn' or 'beta' column in target.ss.",
        call. = FALSE
      )
    }
  }

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
  # target.ss required cols depend on target.ind.
  tss.req <- switch(
    target.ind,
    gwas = c("CHR", "BP", "REF", "ALT", "p", "n", "sgn"),
    corr = c("CHR", "BP", "REF", "ALT", "corr")
  )
  .check_na(target.ss, tss.req, "target.ss")
  .check_na(target.ld, c("CHR", "BP", "REF", "ALT"), "target.ld")
  if (!is.null(external.ss)) {
    .check_na(
      external.ss,
      c("CHR", "BP", "REF", "ALT", external.coef.cols),
      "external.ss"
    )
  }

  # ============================================================
  # 3. GWAS: impute p == 0 and compute corr (done AFTER NA check so
  #    zero-vs-NA distinction is clean)
  # ============================================================
  if (target.ind == "gwas") {
    zero.p <- target.ss$p == 0
    if (any(zero.p)) {
      nonzero.p <- target.ss$p[target.ss$p > 0]
      if (length(nonzero.p) == 0) {
        stop(
          "All p-values in target.ss are zero. ",
          "Cannot compute correlations. Check the input.",
          call. = FALSE
        )
      }
      impute.p <- min(nonzero.p)
      target.ss$p[zero.p] <- impute.p
      if (isTRUE(verbose)) {
        message(sprintf(
          "Imputed %d p-value(s) equal to 0 with smallest non-zero p = %.3g.",
          sum(zero.p), impute.p
        ))
      }
    }
    target.ss$corr <- p2cor(target.ss$p, target.ss$n, target.ss$sgn)
  }

  # ============================================================
  # 4. Normalize allele case and CHR coding; drop invalid CHR
  # ============================================================
  target.ss$REF <- toupper(target.ss$REF)
  target.ss$ALT <- toupper(target.ss$ALT)
  target.ld$REF <- toupper(target.ld$REF)
  target.ld$ALT <- toupper(target.ld$ALT)
  if (!is.null(external.ss)) {
    external.ss$REF <- toupper(external.ss$REF)
    external.ss$ALT <- toupper(external.ss$ALT)
  }

  # Capture original-row indices for target.ld BEFORE any drops, so
  # target.ld.keep always references the user-supplied target.ld.
  target.ld$.orig.idx <- seq_len(nrow(target.ld))

  n.input.tss <- nrow(target.ss)
  n.input.tld <- nrow(target.ld)
  n.input.ess <- if (is.null(external.ss)) 0L else nrow(external.ss)

  if (verbose) {
    message(
      "Input SNPs: ",
      "target.ss=", n.input.tss,
      ", target.ld=", n.input.tld,
      ", external.ss=", n.input.ess
    )
  }

  # Normalize CHR (strip "chr" prefix, map 23/24/26 -> X/Y/MT, "M" -> "MT").
  # Drop rows with unrecognized CHR values.
  pre.tss.chr <- nrow(target.ss)
  target.ss$CHR <- .normalize_chr(target.ss$CHR)
  target.ss     <- .drop_invalid_chr(target.ss, "target.ss", verbose)
  n.dropped.tss.chr <- pre.tss.chr - nrow(target.ss)

  pre.tld.chr <- nrow(target.ld)
  target.ld$CHR <- .normalize_chr(target.ld$CHR)
  target.ld     <- .drop_invalid_chr(target.ld, "target.ld", verbose)
  n.dropped.tld.chr <- pre.tld.chr - nrow(target.ld)

  n.dropped.ess.chr <- 0L
  if (!is.null(external.ss)) {
    pre.ess.chr <- nrow(external.ss)
    external.ss$CHR <- .normalize_chr(external.ss$CHR)
    external.ss     <- .drop_invalid_chr(external.ss, "external.ss", verbose)
    n.dropped.ess.chr <- pre.ess.chr - nrow(external.ss)
  }

  n.dropped.invalid.chr <-
    n.dropped.tss.chr + n.dropped.tld.chr + n.dropped.ess.chr

  # ============================================================
  # 5. Drop strand-ambiguous (optional)
  # ============================================================
  n.dropped.ambiguous <- 0L
  if (isTRUE(drop.ambiguous)) {
    pre.tss <- nrow(target.ss); pre.tld <- nrow(target.ld)
    pre.ess <- if (is.null(external.ss)) 0L else nrow(external.ss)
    target.ss <- .drop_ambiguous(target.ss, "target.ss", verbose)
    target.ld <- .drop_ambiguous(target.ld, "target.ld", verbose)
    if (!is.null(external.ss)) {
      external.ss <- .drop_ambiguous(external.ss, "external.ss", verbose)
    }
    n.dropped.ambiguous <-
      (pre.tss - nrow(target.ss)) +
      (pre.tld - nrow(target.ld)) +
      (pre.ess - (if (is.null(external.ss)) 0L else nrow(external.ss)))
  }

  # ============================================================
  # 6. Drop multi-allelic (duplicate CHR:BP)
  # ============================================================
  pre.tss <- nrow(target.ss); pre.tld <- nrow(target.ld)
  pre.ess <- if (is.null(external.ss)) 0L else nrow(external.ss)

  target.ss <- .drop_multiallelic(target.ss, "target.ss", verbose)
  target.ld <- .drop_multiallelic(target.ld, "target.ld", verbose)
  if (!is.null(external.ss)) {
    external.ss <- .drop_multiallelic(external.ss, "external.ss", verbose)
  }
  n.dropped.multiallelic <-
    (pre.tss - nrow(target.ss)) +
    (pre.tld - nrow(target.ld)) +
    (pre.ess - (if (is.null(external.ss)) 0L else nrow(external.ss)))

  # ============================================================
  # 7. Match target.ss to target.ld
  # ============================================================
  key.ld  <- paste(target.ld$CHR, target.ld$BP, sep = ":")
  key.tss <- paste(target.ss$CHR, target.ss$BP, sep = ":")

  m <- match(key.ld, key.tss)  # for each LD row, index in target.ss (or NA)

  ld.alleles  <- paste(target.ld$REF, target.ld$ALT, sep = "/")
  tss.alleles <- paste(target.ss$REF[m], target.ss$ALT[m], sep = "/")
  tss.flipped <- paste(target.ss$ALT[m], target.ss$REF[m], sep = "/")

  is.match <- !is.na(m) & (ld.alleles == tss.alleles)
  is.flip  <- !is.na(m) & (ld.alleles == tss.flipped) & !is.match
  keep.ld  <- is.match | is.flip

  n.aligned   <- sum(is.match)
  n.flipped   <- sum(is.flip)
  n.unmatched <- sum(!keep.ld)

  if (verbose) {
    message(
      "Aligned ", n.aligned, " SNPs in same orientation, ",
      n.flipped, " in flipped orientation; ",
      n.unmatched, " dropped (no allele match)."
    )
  }

  # Subset LD panel to surviving SNPs, preserving original order
  target.ld.surv <- target.ld[keep.ld, , drop = FALSE]
  ld.keep.idx    <- target.ld.surv$.orig.idx
  target.ld.surv$.orig.idx <- NULL

  # Pull corresponding target.ss rows
  tss.idx   <- m[keep.ld]
  flip.flag <- is.flip[keep.ld]

  target.ss.surv <- target.ss[tss.idx, , drop = FALSE]
  # Re-orient REF/ALT to match LD reference, flip corr where alleles swapped
  target.ss.surv$REF <- target.ld.surv$REF
  target.ss.surv$ALT <- target.ld.surv$ALT
  target.ss.surv$corr[flip.flag] <- -target.ss.surv$corr[flip.flag]

  # Build shared varnames key from LD-canonical alleles.
  varnames <- paste(
    target.ld.surv$CHR, target.ld.surv$BP,
    target.ld.surv$REF, target.ld.surv$ALT,
    sep = ":"
  )
  target.ld.surv$varnames <- varnames
  target.ss.surv$varnames <- varnames

  target.ld.surv <- target.ld.surv[, c("varnames", "CHR", "BP", "REF", "ALT"), drop = FALSE]
  target.ss.surv <- target.ss.surv[, c("varnames", "CHR", "BP", "REF", "ALT", "corr"), drop = FALSE]
  rownames(target.ss.surv) <- NULL
  rownames(target.ld.surv) <- NULL

  n.final <- nrow(target.ld.surv)

  # ============================================================
  # 8. Align external.ss to surviving SNPs
  # ============================================================
  n.external.matched  <- 0L
  n.external.flipped  <- 0L
  n.external.conflict <- 0L
  n.external.zero     <- 0L
  external.ss.out     <- NULL
  external.coef.names <- NULL

  if (!is.null(external.ss)) {
    key.surv <- paste(target.ld.surv$CHR, target.ld.surv$BP, sep = ":")
    key.ess  <- paste(external.ss$CHR, external.ss$BP, sep = ":")

    me <- match(key.surv, key.ess)
    surv.alleles <- paste(target.ld.surv$REF, target.ld.surv$ALT, sep = "/")
    ess.alleles  <- paste(external.ss$REF[me], external.ss$ALT[me], sep = "/")
    ess.flipped  <- paste(external.ss$ALT[me], external.ss$REF[me], sep = "/")

    e.match    <- !is.na(me) & (surv.alleles == ess.alleles)
    e.flip     <- !is.na(me) & (surv.alleles == ess.flipped) & !e.match
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
        "with the surviving target SNP set. This usually indicates: ",
        "(a) target.ss / target.ld and external.ss do not overlap on CHR:BP; ",
        "(b) chromosome or base-pair position coding differs between inputs ",
        "(e.g., hg19 vs hg38 build); or ",
        "(c) allele coding is inconsistent. ",
        "Returning external.ss = NULL. Consider rerunning preprocessS without ",
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
        CHR = target.ld.surv$CHR,
        BP  = target.ld.surv$BP,
        REF = target.ld.surv$REF,
        ALT = target.ld.surv$ALT,
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
  # 9. Return
  # ============================================================
  list(
    target.ss           = target.ss.surv,
    target.ld           = target.ld.surv,
    target.ld.keep      = ld.keep.idx,
    external.ss         = external.ss.out,
    external.coef.names = external.coef.names,
    summary             = c(
      n.input.tss             = n.input.tss,
      n.input.tld             = n.input.tld,
      n.input.ess             = n.input.ess,
      n.dropped.invalid.chr   = n.dropped.invalid.chr,
      n.dropped.ambiguous     = n.dropped.ambiguous,
      n.dropped.multiallelic  = n.dropped.multiallelic,
      n.dropped.unmatched.tss = n.unmatched,
      n.aligned               = n.aligned,
      n.flipped               = n.flipped,
      n.final                 = n.final,
      n.external.matched      = n.external.matched,
      n.external.flipped      = n.external.flipped,
      n.external.conflict     = n.external.conflict,
      n.external.zero         = n.external.zero
    )
  )
}