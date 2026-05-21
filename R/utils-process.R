#' @noRd
.standardize_cols <- function(df, mapping, required.keys, lab) {
  # Metadata keys (CHR/BP/REF/ALT) are uppercased; value keys stay lowercase.
  upper.keys <- c("chr", "bp", "ref", "alt")
  for (key in required.keys) {
    user.name <- mapping[[key]]
    if (is.null(user.name) || is.na(user.name)) {
      stop(sprintf(
        "'%s' column for %s is required but no mapping was provided.",
        key, lab), call. = FALSE
      )
    }
    if (!(user.name %in% colnames(df))) {
      stop(sprintf(
        "Column '%s' (mapped from key '%s') not found in %s.",
        user.name, key, lab), call. = FALSE
      )
    }
    new.name <- if (key %in% upper.keys) toupper(key) else key
    colnames(df)[colnames(df) == user.name] <- new.name
  }
  df
}

#' @noRd
.drop_multiallelic <- function(df, lab, verbose) {
  key <- paste(df$CHR, df$BP, sep = ":")
  dups <- key %in% key[duplicated(key)]
  if (any(dups)) {
    n.dropped <- sum(dups)
    if (isTRUE(verbose)) {
      warning(sprintf(
        "Dropped %d multi-allelic SNPs (duplicated CHR:BP) from %s.",
        n.dropped, lab), call. = FALSE
      )
    }
    df <- df[!dups, , drop = FALSE]
  }
  df
}

#' @noRd
.drop_ambiguous <- function(df, lab, verbose) {
  ambig <- (df$REF == "A" & df$ALT == "T") |
           (df$REF == "T" & df$ALT == "A") |
           (df$REF == "C" & df$ALT == "G") |
           (df$REF == "G" & df$ALT == "C")
  if (any(ambig)) {
    n.dropped <- sum(ambig)
    if (isTRUE(verbose)) {
      message(sprintf("Dropped %d strand-ambiguous SNPs from %s.", n.dropped, lab))
    }
    df <- df[!ambig, , drop = FALSE]
  }
  df
}

#' Normalize chromosome coding to canonical form.
#' Strips "chr" prefix, maps PLINK numeric (23/24/26) to letters,
#' normalizes "M" to "MT". Anything not recognized is passed through
#' uppercased — pair with .drop_invalid_chr() to filter those.
#' @noRd
.normalize_chr <- function(chr) {
  chr <- as.character(chr)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr <- toupper(chr)
  chr[chr == "23"] <- "X"
  chr[chr == "24"] <- "Y"
  chr[chr == "25"] <- "XY"
  chr[chr == "26"] <- "MT"
  chr[chr == "M"]  <- "MT"
  chr
}

#' Canonical set of recognized chromosomes (human-genome focus).
#' @noRd
.valid_chr <- c(as.character(1:22), "X", "Y", "MT", "XY")

#' Drop rows whose CHR is not in the recognized set, with a warning
#' that names the unrecognized values found.
#' @noRd
.drop_invalid_chr <- function(df, lab, verbose) {
  is.invalid <- !(df$CHR %in% .valid_chr)
  if (any(is.invalid)) {
    bad <- unique(df$CHR[is.invalid])
    n.dropped <- sum(is.invalid)
    if (isTRUE(verbose)) {
      warning(sprintf(
        "Dropped %d SNPs from %s with unrecognized CHR values: %s%s.",
        n.dropped, lab,
        paste(head(bad, 10), collapse = ", "),
        if (length(bad) > 10) sprintf(" (and %d more)", length(bad) - 10) else ""
      ), call. = FALSE)
    }
    df <- df[!is.invalid, , drop = FALSE]
  }
  df
}

#' @noRd
.check_na <- function(df, cols, lab) {
  for (col in cols) {
    if (any(is.na(df[[col]]))) {
      n <- sum(is.na(df[[col]]))
      stop(sprintf(
        "%d NA value(s) found in column '%s' of %s. ",
        n, col, lab),
        "Clean the input (e.g., na.omit() or impute) before processing.",
        call. = FALSE
      )
    }
  }
}