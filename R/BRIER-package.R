#' BRIER2: Bregman-Divergence-Based Regularized Integrative Estimator in Regression
#'
#' Penalized regression methods for genetic risk prediction that integrate
#' pretrained external models or external individual-level data with target
#' cohort data. Supports Gaussian, binomial, and Poisson families with LASSO,
#' MCP, and SCAD penalties. Includes individual-level (BRIERi),
#' summary-statistics (BRIERs), and full individual-level external data
#' (BRIERfull) variants.
#'
#' @keywords internal
#' @useDynLib BRIER, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom stats approx cor glm logLik model.matrix optim prcomp predict qt residuals sd var
#' @importFrom rlang .data
"_PACKAGE"