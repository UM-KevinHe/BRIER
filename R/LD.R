#' Calculate LD matrix for BRIER.S
#'
#' Compute a sparse LD reference matrix from a genotype reference panel, with
#' optional soft-thresholding and block-structured LD pruning. Used as input
#' to \code{\link{BRIERs}}.
#'
#' @param X An n x p reference genotype matrix (individual-level).
#' @param SNP.info A data.frame describing the columns of \code{X}, in
#'   order. Required when \code{LDB} is supplied; must contain \code{CHR}
#'   and \code{BP} columns for assigning SNPs to LD blocks.
#' @param LDB Optional m x 3 matrix of LD block boundaries (\code{chr},
#'   \code{start}, \code{stop}), e.g. from Berisa et al. If NULL, the full
#'   p x p LD matrix is computed.
#' @param tau A numeric scalar. Soft-threshold parameter applied to LD
#'   entries; entries with |LD| < tau are set to zero. Defaults to 0
#'   (no thresholding).
#'
#' @return An object of class \code{"LD"} containing:
#' \describe{
#'   \item{XtX}{The p x p sparse LD matrix.}
#'   \item{blk}{The block boundary vector (NULL if no \code{LDB} supplied).}
#'   \item{nz}{Integer vector of non-constant column indices retained.}
#' }
#'
#' @seealso \code{\link{BRIERs}}, \code{\link{p2cor}}
#'
#' @examples
#' \dontrun{
#' # full LD matrix
#' ld <- calLD(X, tau = 0.05)
#'
#' # block-structured LD with Berisa et al. blocks
#' ld <- calLD(X, SNP.info = SNP.info, LDB = LDB, tau = 0.05)
#' }
#'
#' @export
calLD <- function(X, SNP.info = NULL, LDB = NULL, tau = 0) {

  if (is.null(LDB)) {

    XX.list <- standardize_X(X)
    std.X  <- XX.list[[1]]
    center <- XX.list[[2]]
    scale  <- XX.list[[3]]
    nz <- which(scale > 1e-6)
    if (length(nz) != ncol(X)) {
      removed <- which(scale <= 1e-6)
      std.X  <- std.X[, nz, drop = FALSE]
      center <- center[nz]
      scale  <- scale[nz]
      warning(
        "Constant columns are removed: ", paste(removed, collapse = ", "),
        call. = FALSE
      )
    }
    XtX <- tcrossprod(t(std.X)) / nrow(std.X)
    XtX[abs(XtX) < tau] <- 0
    out <- list(
      XtX = Matrix::drop0(Matrix::Matrix(XtX, sparse = TRUE)),
      blk = NULL,
      nz  = nz
    )
  } else {
    if (ncol(X) != nrow(SNP.info)) {
      stop(
        "The dimension of reference X and SNP.info does not match.",
        call. = FALSE
      )
    }
    chr <- SNP.info$CHR
    pos <- SNP.info$BP
    if (is.null(chr) || is.null(pos)) {
      stop("SNP.info must contain CHR and BP columns.", call. = FALSE)
    }

    ## Calculate LD matrix from reference data
    blk <- LD_blocks(chr, pos, LDB)
    XX.list <- standardize_X(X)
    std.X  <- XX.list[[1]]
    center <- XX.list[[2]]
    scale  <- XX.list[[3]]
    nz <- which(scale > 1e-6)
    if (length(nz) != ncol(X)) {
      removed <- which(scale <= 1e-6)
      std.X  <- std.X[, nz, drop = FALSE]
      center <- center[nz]
      scale  <- scale[nz]
      warning(
        "Constant columns are removed: ", paste(removed, collapse = ", "),
        call. = FALSE
      )
    }
    XtX <- LD_sigma(X = std.X, blk = blk, tau = tau)
    out <- list(
      XtX = XtX,
      blk = blk,
      nz  = nz
    )
  }

  class(out) <- "LD"
  out
}

#' Path to a bundled Berisa-Pickrell LD block file
#'
#' Returns the file system path to a packaged BED file that describes
#' approximately-independent LD blocks from Berisa and Pickrell (2016).
#' The package ships LD block coordinates for three continental ancestries
#' (AFR, EAS, EUR) in two genome builds (hg19, hg38).
#'
#' These files can be used directly as the \code{blocks} argument to
#' \code{\link{LD_blocks}} or supplied to downstream LD-aware functions
#' such as \code{\link{calLD}}.
#'
#' @param ancestry Character. One of \code{"AFR"}, \code{"EAS"}, or
#'   \code{"EUR"}.
#' @param build Character. One of \code{"hg19"} (GRCh37) or \code{"hg38"}
#'   (GRCh38). Default \code{"hg38"}.
#'
#' @return Absolute path to the requested BED file.
#'
#' @references
#' Berisa, T. and Pickrell, J. K. (2016). Approximately independent
#' linkage disequilibrium blocks in human populations.
#' \emph{Bioinformatics}, 32(2), 283-285.
#'
#' Original hg19 BED files:
#' \url{https://bitbucket.org/nygcresearch/ldetect-data}
#'
#' @examples
#' \dontrun{
#' bed_path <- getLDB("EUR", "hg38")
#' blocks <- read.table(bed_path, header = T)
#' head(blocks)
#' }
#'
#' @seealso \code{\link{LD_blocks}}, \code{\link{calLD}}
#'
#' @export
getLDB <- function(
  ancestry = c("AFR", "EAS", "EUR"),
  build    = c("hg38", "hg19")
  ) {
  ancestry <- match.arg(ancestry)
  build    <- match.arg(build)
  fname    <- sprintf("Berisa.%s.%s.bed", ancestry, build)
  path     <- system.file("extdata", fname, package = "BRIER")
  if (!nzchar(path)) {
    stop(
      "LD block file not found in installed package: ", fname, ". ",
      "Try reinstalling BRIER.", call. = FALSE
    )
  }
  path
}