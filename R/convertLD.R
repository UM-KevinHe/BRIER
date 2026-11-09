#' Restrict an LD matrix to chromosome or LD-block structure
#'
#' Takes an LD matrix (dense or sparse, for one chromosome or several) and sets to
#' zero every entry between two SNPs on different chromosomes (\code{to = "chr"})
#' or in different LD blocks (\code{to = "block"}), returning the same kind of
#' object as \code{\link{calLD}}. Entries inside a chromosome or block are kept
#' unchanged, so the result is exactly the chromosome-wise or block-wise LD that
#' \code{\link{calLD}} would compute from the same genotypes.
#'
#' The conversion only removes entries. It goes from a dense LD to chromosome LD
#' to block LD, never back: LD between blocks that a panel does not store cannot
#' be recovered. In particular, a panel shipped as blocks (for example the PRS-CS
#' UK Biobank panels, see \code{\link{convertLDref}}) cannot be turned into
#' chromosome LD or into a different block partition.
#'
#' Keeping the diagonal blocks of a positive semidefinite matrix keeps it positive
#' semidefinite (the block mask is itself positive semidefinite), so a valid LD
#' matrix stays valid.
#'
#' @param LD A p x p symmetric LD matrix: a base \code{matrix} or a \code{Matrix}
#'   (for example the \code{XtX} of \code{\link{calLD}}).
#' @param SNP.info A data.frame describing the rows of \code{LD}, in order, with
#'   columns \code{CHR} (numeric, 1 to 22) and \code{BP}. SNPs must be sorted by
#'   chromosome, then position, as \code{\link{calLD}} expects.
#' @param to \code{"block"} (default) or \code{"chr"}.
#' @param LDB For \code{to = "block"}, an m x 3 matrix of LD block boundaries
#'   (\code{chr}, \code{start}, \code{stop}) in the build of \code{SNP.info$BP},
#'   with numeric chromosomes, as for \code{\link{calLD}}; for example
#'   \code{\link{getLDB}()} read with \code{read.table()}. Ignored for \code{"chr"}.
#'
#' @return An object of class \code{"LD"}: \code{XtX} (the restricted p x p
#'   \code{dgCMatrix}, dimnames kept), \code{blk} (0-based block start indices with
#'   the end marker p, as \code{calLD}; for \code{"chr"} one block per chromosome)
#'   and \code{nz} (\code{seq_len(p)}: no column is dropped).
#'
#' @seealso \code{\link{calLD}}, \code{\link{getLDB}}, \code{\link{convertLDref}}
#'
#' @examples
#' \dontrun{
#' ldb = read.table(getLDB("AFR", "hg38"), header = TRUE)
#' ldb[, 1] = as.numeric(sub("chr", "", ldb[, 1]))
#' dense = calLD(X_ref)                                # full LD
#' blk = convertLD(dense$XtX, SNP.info, to = "block", LDB = as.matrix(ldb))
#' chr = convertLD(dense$XtX, SNP.info, to = "chr")
#' }
#'
#' @export
convertLD <- function(LD, SNP.info, to = c("block", "chr"), LDB = NULL) {
  to <- match.arg(to)
  if (!is.matrix(LD) && !methods::is(LD, "Matrix")) {
    stop("LD must be a matrix or a Matrix.", call. = FALSE)
  }
  p <- ncol(LD)
  if (nrow(LD) != p) { stop("LD must be square.", call. = FALSE) }
  if (!is.data.frame(SNP.info) || is.null(SNP.info$CHR) || is.null(SNP.info$BP)) {
    stop("SNP.info must be a data.frame with CHR and BP columns.", call. = FALSE)
  }
  if (nrow(SNP.info) != p) {
    stop("SNP.info has ", nrow(SNP.info), " rows for an LD of dimension ", p, ".", call. = FALSE)
  }
  chr <- suppressWarnings(as.numeric(SNP.info$CHR))
  pos <- as.numeric(SNP.info$BP)
  if (anyNA(chr) || anyNA(pos)) {
    stop("SNP.info$CHR and SNP.info$BP must be numeric, with no NA.", call. = FALSE)
  }
  if (is.unsorted(order(chr, pos)) || any(duplicated(chr[c(TRUE, diff(chr) != 0)]))) {
    stop("SNP.info must be sorted by chromosome, then position.", call. = FALSE)
  }

  if (to == "chr") {
    starts <- which(c(TRUE, diff(chr) != 0)) - 1L
    blk <- c(as.integer(starts), as.integer(p))
  } else {
    if (is.null(LDB)) { stop("to = \"block\" needs LDB (see getLDB).", call. = FALSE) }
    LDB <- as.matrix(LDB)
    if (ncol(LDB) < 3 || !is.numeric(LDB)) {
      stop("LDB must be a numeric m x 3 matrix (chr, start, stop).", call. = FALSE)
    }
    LDB <- LDB[order(LDB[, 1], LDB[, 2]), 1:3, drop = FALSE]
    blk <- as.integer(LD_blocks(chr, pos, LDB))
  }
  id <- findInterval(seq_len(p) - 1L, blk[-length(blk)])

  M <- if (is.matrix(LD)) {
    Matrix::Matrix(LD, sparse = TRUE)
  } else {
    LD
  }
  M <- methods::as(methods::as(methods::as(M, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  i <- M@i + 1L
  j <- rep.int(seq_len(p), diff(M@p))
  keep <- id[i] == id[j]
  out <- Matrix::sparseMatrix(
    i = i[keep],
    j = j[keep],
    x = M@x[keep],
    dims = c(p, p),
    dimnames = dimnames(LD)
  )
  out <- list(
    XtX = Matrix::drop0(out),
    blk = blk,
    nz  = seq_len(p)
  )
  class(out) <- "LD"
  out
}

#' Convert a published LD reference panel onto your SNP set
#'
#' Reads an LD reference panel distributed by another method and converts it,
#' one chromosome at a time, onto the SNPs and allele coding of your analysis.
#' Currently \code{format = "prscs"}: the PRS-CS / PRS-CSx panels (for example
#' \code{ldblk_ukbb_afr}), one HDF5 file per chromosome
#' (\code{ldblk_<panel>_chr<c>.hdf5}, groups \code{blk_<k>} with \code{ldblk} and
#' \code{snplist}) and the panel's SNP table (\code{snpinfo_<panel>_hm3}, columns
#' \code{CHR SNP BP A1 A2 MAF}).
#'
#' The panel's own LD blocks are kept: PRS-CS stores only the LD within each block,
#' so the result is block LD on PRS-CS's partition and cannot be converted to
#' chromosome LD or to other blocks (see \code{\link{convertLD}}). Signs are
#' aligned to the \code{ALT} allele of \code{SNP.info}: an entry r_ij is multiplied
#' by s_i s_j, with s = +1 where the panel's A1 is your ALT and -1 where it is your
#' REF. A SNP of \code{SNP.info} that the panel lacks, or whose alleles differ, is
#' not filled in: it is reported in \code{missing}, and the LD covers the rest.
#' Nothing is shrunk.
#'
#' Requires the Bioconductor package \pkg{rhdf5}.
#'
#' @param dir The panel directory (holding the \code{.hdf5} files and the
#'   \code{snpinfo_*} table).
#' @param SNP.info A data.frame of your SNPs on ONE chromosome, in your order:
#'   \code{SNP} (rsID), \code{CHR}, \code{REF}, \code{ALT}.
#' @param format The panel's format. Only \code{"prscs"} for now.
#'
#' @return A list: \code{XtX} (a \code{dgCMatrix} over the SNPs found, dimnames
#'   the rsIDs, in the order of \code{SNP.info}), \code{keep} (the rows of
#'   \code{SNP.info} it covers), \code{missing} (the rows it does not, with the
#'   reason), \code{blocks} (the number of panel blocks used) and
#'   \code{sign.flipped} (how many SNPs had their sign flipped).
#'
#' @seealso \code{\link{convertLD}}, \code{\link{calLD}}, \code{\link{BRIERs}}
#'
#' @examples
#' \dontrun{
#' ld1 = convertLDref("ldblk_ukbb_afr", SNP.info = snps[snps$CHR == 1, ])
#' snps1 = snps[snps$CHR == 1, ][ld1$keep, ]
#' }
#'
#' @export
convertLDref <- function(dir, SNP.info, format = c("prscs")) {
  format <- match.arg(format)
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("convertLDref() needs the Bioconductor package rhdf5.", call. = FALSE)
  }
  need <- c("SNP", "CHR", "REF", "ALT")
  if (!is.data.frame(SNP.info) || !all(need %in% names(SNP.info))) {
    stop("SNP.info must be a data.frame with columns SNP, CHR, REF, ALT.", call. = FALSE)
  }
  chr <- unique(as.numeric(SNP.info$CHR))
  if (length(chr) != 1L || is.na(chr)) {
    stop("SNP.info must hold SNPs of exactly one chromosome.", call. = FALSE)
  }
  if (anyDuplicated(SNP.info$SNP)) { stop("SNP.info$SNP has duplicates.", call. = FALSE) }

  info_file <- list.files(dir, pattern = "^snpinfo_.*", full.names = TRUE)
  if (length(info_file) != 1L) {
    stop("Expected one snpinfo_* table in ", dir, ", found ", length(info_file), ".", call. = FALSE)
  }
  h5 <- list.files(dir, pattern = sprintf("_chr%d\\.hdf5$", as.integer(chr)), full.names = TRUE)
  if (length(h5) != 1L) {
    stop("Expected one *_chr", chr, ".hdf5 in ", dir, ", found ", length(h5), ".", call. = FALSE)
  }
  info <- utils::read.table(info_file, header = TRUE, stringsAsFactors = FALSE)
  info <- info[info$CHR == chr, , drop = FALSE]

  s <- SNP.info
  k <- match(s$SNP, info$SNP)
  a1 <- info$A1[k]
  a2 <- info$A2[k]
  same <- !is.na(k) & a1 == s$ALT & a2 == s$REF
  swap <- !is.na(k) & a1 == s$REF & a2 == s$ALT
  reason <- ifelse(is.na(k), "not in the panel", ifelse(same | swap, NA_character_, "allele mismatch"))
  sg <- ifelse(swap, -1, 1)

  ls5 <- rhdf5::h5ls(h5)
  blks <- unique(sub("^/", "", ls5$group[grepl("^/blk_", ls5$group)]))
  I <- list()
  J <- list()
  V <- list()
  seen <- rep(FALSE, nrow(s))
  nb <- 0L
  for (b in blks) {
    sn <- as.character(rhdf5::h5read(h5, paste0(b, "/snplist")))
    kk <- match(sn, s$SNP)
    use <- which(!is.na(kk) & is.na(reason[kk]))
    if (!length(use)) { next }
    M <- as.matrix(rhdf5::h5read(h5, paste0(b, "/ldblk")))
    if (nrow(M) != length(sn)) {
      stop(b, ": ldblk is ", nrow(M), " x ", ncol(M), " for ", length(sn), " SNPs.", call. = FALSE)
    }
    kk <- kk[use]
    M <- M[use, use, drop = FALSE] * outer(sg[kk], sg[kk])
    ij <- which(M != 0, arr.ind = TRUE)
    I[[b]] <- kk[ij[, 1]]
    J[[b]] <- kk[ij[, 2]]
    V[[b]] <- M[ij]
    if (any(seen[kk])) { stop(b, ": SNPs already seen in an earlier block.", call. = FALSE) }
    seen[kk] <- TRUE
    nb <- nb + 1L
  }
  rhdf5::h5closeAll()
  reason[is.na(reason) & !seen] <- "in no panel block"

  keep <- which(is.na(reason))
  R <- Matrix::sparseMatrix(
    i = unlist(I),
    j = unlist(J),
    x = unlist(V),
    dims = c(nrow(s), nrow(s))
  )
  R <- R[keep, keep, drop = FALSE]
  R <- methods::as(methods::as(R, "generalMatrix"), "CsparseMatrix")
  dimnames(R) <- list(s$SNP[keep], s$SNP[keep])
  list(
    XtX = R,
    keep = keep,
    missing = data.frame(
      row = which(!is.na(reason)),
      SNP = s$SNP[!is.na(reason)],
      reason = reason[!is.na(reason)],
      stringsAsFactors = FALSE
    ),
    blocks = nb,
    sign.flipped = sum(sg[keep] < 0)
  )
}
