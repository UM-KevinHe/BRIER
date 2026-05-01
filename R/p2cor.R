#' Convert GWAS p-values to marginal correlations
#'
#' Recover the marginal correlation between genotype and outcome from GWAS
#' summary statistics (p-values, sample sizes, and effect-size signs). Useful
#' for constructing the \code{corr} column of \code{sumstats} input to
#' \code{\link{BRIERs}}.
#'
#' @param p A numeric vector of p-values from marginal associations.
#' @param n A numeric scalar or vector of sample sizes. If scalar, the same
#'   sample size is used for all p-values.
#' @param sign A numeric vector of the signs of the marginal coefficient
#'   estimates (+1 or -1). Defaults to +1 for all.
#'
#' @return A numeric vector of marginal correlations (X'y / n).
#'
#' @seealso \code{\link{BRIERs}}
#'
#' @examples
#' p <- c(0.01, 0.015)
#' n <- c(100, 1000)
#' sign <- c(-1, 1)
#' p2cor(p, n, sign)
#'
#' @export
p2cor <- function(p, n, sign = rep(1, length(p))) {

  stopifnot(length(n) == 1 || length(n) == length(p))
  stopifnot(length(p) == length(sign))

  t <- sign(sign) * qt(p / 2, df = n - 2, lower.tail = FALSE)

  t / sqrt(n - 2 + t^2)
}