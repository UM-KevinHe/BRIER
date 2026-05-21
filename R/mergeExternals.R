#' Merge multiple external model coefficients
#'
#' `mergeExternals()` combines multiple external SNP-level coefficient files
#' into a single wide-format data.frame. The merged output can be used as the
#' `external.ss` input for downstream preprocessing functions such as
#' `preprocessI()` or `preprocessS()`.
#'
#' Each element of `external.list` should be a data.frame containing SNP
#' annotation columns `CHR`, `BP`, `REF`, `ALT`, and a coefficient column
#' `coef`. Additional columns are ignored. Each data.frame represents one
#' external model.
#'
#' First, `mergeExternals()` performs quality control and harmonization
#' independently for each external model. Chromosome labels are normalized by
#' removing `"chr"` prefixes, converting PLINK chromosome codes `23`, `24`, `25`,
#' and `26` to `X`, `Y`, `XY`, and `MT`, and converting `M` to `MT`. SNPs with
#' unrecognized chromosome labels are removed with a warning. The function errors 
#' out if any required column contains missing values, removes sites with 
#' duplicated `CHR:BP` entries within each input (treated as multi-allelic 
#' or input duplicates), and optionally removes strand-ambiguous SNPs, 
#' including A/T and C/G variants, when `drop.ambiguous = TRUE`.
#'
#' After quality control, the function constructs the union of SNPs across all
#' external models. For each SNP, matched by `CHR:BP`, the first external model 
#' that survives QC for that SNP defines the canonical `REF`/`ALT` orientation in
#' the merged output. Coefficients from all other external models are then
#' aligned to this canonical orientation. If the alleles match the canonical
#' orientation, the coefficient is retained unchanged. If the alleles are
#' reversed relative to the canonical orientation, the coefficient sign is
#' flipped. If the alleles conflict, the coefficient is set to zero. If a SNP is
#' absent from a given external model, its coefficient for that model is also
#' set to zero.
#'
#' The returned data.frame contains columns `CHR`, `BP`, `REF`, `ALT`, followed
#' by standardized coefficient columns `coef1`, `coef2`, ..., `coefM`, where
#' `M` is the number of external models. Rows are sorted by chromosome and
#' base-pair position. The output also has a `coef.names` attribute mapping
#' standardized coefficient names back to the input model names, or to
#' automatically generated names `model1`, `model2`, ..., when the input list is
#' unnamed.
#'
#' @param external.list A non-empty list of data.frames. Each data.frame
#'   must contain columns \code{CHR}, \code{BP}, \code{REF}, \code{ALT},
#'   \code{coef}. The list may be named — those names are used as labels
#'   for the per-model verbose output and stored on the result via the
#'   \code{coef.names} attribute.
#' @param drop.ambiguous Logical. If TRUE (default), drops strand-ambiguous
#'   SNPs from each input before merging.
#' @param verbose Logical. If TRUE (default), prints SNP counts plus a
#'   per-model alignment breakdown.
#'
#' @return A data.frame with columns \code{CHR}, \code{BP}, \code{REF},
#'   \code{ALT}, \code{coef1}, \code{coef2}, ..., \code{coefM}, where M is
#'   the length of \code{external.list}. Rows are sorted by \code{(CHR, BP)}.
#'   The result has an attribute \code{coef.names} mapping the standardized
#'   \code{coef1}, \code{coef2}, ... back to the input list names (or
#'   \code{model1}, \code{model2}, ... if the list was unnamed).
#'
#' Pass this directly as `external.ss` to \code{\link{preprocessI}} or
#' \code{\link{preprocessS}} with
#' \code{external.coef.cols = paste0("coef", seq_len(M))}.
#'
#' @seealso \code{\link{preprocessI}}, \code{\link{preprocessS}}
#'
#' @examples
#' \dontrun{
#' eur <- data.frame(
#'   CHR = c(1,1,2), BP = c(100,200,300),
#'   REF = c("A","C","G"), ALT = c("G","T","A"),
#'   coef = c(0.5, -0.2, 0.1)
#' )
#' eas <- data.frame(
#'   CHR = c(1,2,2), BP = c(100,300,400),
#'   REF = c("G","A","T"), ALT = c("A","G","C"),   # # rows 1 & 2: REF/ALT flipped vs eur
#'   coef = c(0.3, 0.4, -0.1)
#' )
#'
#' merged <- mergeExternals(list(EUR = eur, EAS = eas))
#' merged
#' attr(merged, "coef.names")     # c(coef1 = "EUR", coef2 = "EAS")
#'
#' # Pass to preprocessS
#' out <- preprocessS(
#'   target.ss          = my_gwas,
#'   target.ind         = "gwas",
#'   target.ld          = my_ld_panel,
#'   external.ss        = merged,
#'   external.coef.cols = c("coef1", "coef2")
#' )
#' }
#'
#' @export
mergeExternals <- function(
  external.list,
  drop.ambiguous = TRUE,
  verbose = TRUE
) {

  # ============================================================
  # 1. Validate input
  # ============================================================
  if (!is.list(external.list) || length(external.list) == 0) {
    stop("external.list must be a non-empty list of data.frames.", call. = FALSE)
  }

  M <- length(external.list)
  list.names <- names(external.list)
  if (is.null(list.names)) list.names <- rep("", M)

  required.cols <- c("CHR", "BP", "REF", "ALT", "coef")

  # Friendly per-element label for messages
  make.lab <- function(i) {
    if (nzchar(list.names[i])) list.names[i] else sprintf("external.list[[%d]]", i)
  }

  # Coerce, validate columns, uppercase alleles
  std.list <- vector("list", M)
  for (i in seq_len(M)) {
    df <- as.data.frame(external.list[[i]])
    missing <- setdiff(required.cols, colnames(df))
    if (length(missing) > 0) {
      stop(sprintf(
        "%s is missing required column(s): %s. ",
        make.lab(i), paste(missing, collapse = ", ")),
        "Each data.frame must contain exactly: CHR, BP, REF, ALT, coef.",
        call. = FALSE
      )
    }
    df <- df[, required.cols, drop = FALSE]
    df$REF <- toupper(df$REF)
    df$ALT <- toupper(df$ALT)
    std.list[[i]] <- df
  }

  # ============================================================
  # 2. NA check on required columns
  # ============================================================
  for (i in seq_len(M)) {
    .check_na(std.list[[i]], required.cols, make.lab(i))
  }

  # ============================================================
  # 3. Normalize CHR coding and drop invalid CHR
  # ============================================================
  for (i in seq_len(M)) {
    std.list[[i]]$CHR <- .normalize_chr(std.list[[i]]$CHR)
    std.list[[i]]     <- .drop_invalid_chr(std.list[[i]], make.lab(i), verbose)
  }

  n.input <- vapply(std.list, nrow, integer(1))
  if (verbose) {
    message(
      "Input SNPs per model (after CHR validation): ",
      paste(sprintf(
        "%s=%d", vapply(seq_len(M), make.lab, character(1)), n.input
        ), collapse = ", ")
    )
  }

  # ============================================================
  # 4. Drop strand-ambiguous (optional) and multi-allelic
  # ============================================================
  if (isTRUE(drop.ambiguous)) {
    for (i in seq_len(M)) {
      std.list[[i]] <- .drop_ambiguous(std.list[[i]], make.lab(i), verbose)
    }
  }
  for (i in seq_len(M)) {
    std.list[[i]] <- .drop_multiallelic(std.list[[i]], make.lab(i), verbose)
  }

  # ============================================================
  # 5. Build union SNP set with first-wins canonical orientation
  # ============================================================
  # Walk the list in order. New CHR:BP keys are added with the first
  # model's REF/ALT defining canonical orientation. Repeated keys are
  # handled later in the per-model alignment step.
  union.chr <- vector(typeof(std.list[[1]]$CHR), 0)
  union.bp  <- vector(typeof(std.list[[1]]$BP),  0)
  union.ref <- character(0)
  union.alt <- character(0)
  union.key <- character(0)

  for (i in seq_len(M)) {
    df <- std.list[[i]]
    df.key <- paste(df$CHR, df$BP, sep = ":")
    new <- !(df.key %in% union.key)
    if (any(new)) {
      union.chr <- c(union.chr, df$CHR[new])
      union.bp  <- c(union.bp,  df$BP[new])
      union.ref <- c(union.ref, df$REF[new])
      union.alt <- c(union.alt, df$ALT[new])
      union.key <- c(union.key, df.key[new])
    }
  }

  # Sort by CHR, BP for deterministic ordering
  ord <- order(union.chr, union.bp)
  union.chr <- union.chr[ord]
  union.bp  <- union.bp[ord]
  union.ref <- union.ref[ord]
  union.alt <- union.alt[ord]
  union.key <- union.key[ord]
  n.union <- length(union.key)

  if (verbose) {
    message("Union SNP set: ", n.union, " variants across ", M, " model(s).")
  }

  # ============================================================
  # 6. Align each model to canonical, fill coef matrix
  # ============================================================
  std.names <- paste0("coef", seq_len(M))
  coef.mat  <- matrix(0, nrow = n.union, ncol = M)
  colnames(coef.mat) <- std.names

  union.alleles <- paste(union.ref, union.alt, sep = "/")

  n.matched  <- integer(M)
  n.flipped  <- integer(M)
  n.conflict <- integer(M)
  n.zero     <- integer(M)

  for (i in seq_len(M)) {
    df <- std.list[[i]]
    df.key <- paste(df$CHR, df$BP, sep = ":")

    m <- match(union.key, df.key)

    df.alleles <- paste(df$REF[m], df$ALT[m], sep = "/")
    df.flipped <- paste(df$ALT[m], df$REF[m], sep = "/")

    is.match    <- !is.na(m) & (union.alleles == df.alleles)
    is.flip     <- !is.na(m) & (union.alleles == df.flipped) & !is.match
    is.conflict <- !is.na(m) & !is.match & !is.flip
    is.zero     <- is.na(m)

    n.matched[i]  <- sum(is.match)
    n.flipped[i]  <- sum(is.flip)
    n.conflict[i] <- sum(is.conflict)
    n.zero[i]     <- sum(is.zero)

    if (any(is.match)) coef.mat[is.match, i] <-  df$coef[m[is.match]]
    if (any(is.flip))  coef.mat[is.flip, i]  <- -df$coef[m[is.flip]]
    # is.conflict and is.zero stay at 0
  }

  # ============================================================
  # 7. Resolve labels for the coef.names attribute
  # ============================================================
  raw.names <- vapply(seq_len(M), function(i) {
    if (nzchar(list.names[i])) list.names[i] else sprintf("model%d", i)
  }, character(1))
  raw.names <- make.unique(raw.names, sep = ".")
  coef.names.attr <- setNames(raw.names, std.names)

  if (verbose) {
    for (i in seq_len(M)) {
      message(sprintf(
        "  %s (-> coef%d): matched=%d, flipped=%d, conflict=%d, zero-padded=%d",
        raw.names[i], i,
        n.matched[i], n.flipped[i], n.conflict[i], n.zero[i]
      ))
    }
  }

  # ============================================================
  # 8. Build output data.frame
  # ============================================================
  out <- data.frame(
    CHR = union.chr,
    BP  = union.bp,
    REF = union.ref,
    ALT = union.alt,
    coef.mat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- NULL
  attr(out, "coef.names") <- coef.names.attr
  out
}