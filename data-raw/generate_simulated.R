# =============================================================================
# generate_simulated.R: regenerate the package's three simulated example datasets
# (Data_BRIERfull, Data_BRIERi, Data_BRIERs). Ray, 2026-11-09.
#
#   Rscript generate_simulated.R [out_dir] [compare_dir]
#
# Recovered from the original scripts (archive/BRIER_test/example_data/scripts/
# generate_example.R, BRIERi_example.R, BRIERs_example.R, written for a pre-release build
# "BRIER2"); same design, same seeds, consolidated into one file without setwd() or unused
# libraries:
#   p = 200 independent N(0, 1) predictors, signal at columns 20, 40, ..., 200, a sparse linear
#   model with Gaussian noise at hsq = 0.2; target n 150 / 150 / 150 (train / validation /
#   testing, seeds 1 / 2 / 3); three externals with n 300 training (seeds 10, 20, 30) and 150
#   validation (seeds 11, 21, 31); outcome noise seed = split seed + 10000. External 1 shares the
#   target's coefficients, externals 2 and 3 differ.
#   Data_BRIERfull  the four cohorts, individual level, plus training summary statistics
#   Data_BRIERi     the target, and beta.external: each external's lasso from BRIERi at eta = 0
#                   (intercept row first), lambda by that external's validation MSPE
#   Data_BRIERs     the target, and beta.external: each external's lasso from BRIERs at eta = 0 on
#                   its training summary statistics and in-sample LD, lambda by its validation MSPE
# With compare_dir, each object is compared with the .rda of the same name there (the package's
# data/); Data_BRIERfull is pure simulation and must be identical; the two beta.external matrices
# come from fits, so they are compared numerically.
# =============================================================================
suppressWarnings(suppressMessages({ library(BRIER); library(Matrix) }))
a = commandArgs(trailingOnly = TRUE)
OUT = if (length(a) >= 1) a[1] else "data"
CMP = if (length(a) >= 2) a[2] else NA_character_
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("BRIER %s\n", as.character(packageVersion("BRIER"))))

# -- 0. design ------------------------------------------------------------------------
p = 200
hsq = 0.2
true_id = (1:10) * 20
target_beta = c(0.5, 1, 1.5, 2, 2.5, -0.5, -1, -1.5, -2, -2.5)
ext1_beta = c(0.5, 1, 1.5, 2, 2.5, -0.5, -1, -1.5, -2, -2.5)
ext2_beta = c(1.2, 1.7, 0.8, 1.3, 1.8, -1.2, -1.7, -0.8, -1.3, -1.8)
ext3_beta = c(1.9, 2.4, 0.1, 0.6, 1.1, -1.9, -2.4, -0.1, -0.6, -1.1)

# -- 1. helpers (as the original) ---------------------------------------------------------
gen_X = function(n, p, mu = 0, sd = 1, seed = NULL) {
  if (!is.null(seed)) { set.seed(seed) }
  X = matrix(rnorm(n * p, mean = mu, sd = sd), nrow = n, ncol = p)
  colnames(X) = paste0("X", 1:p)
  X
}
gen_y_gaussian = function(X, true_id, true_beta, hsq, seed = NULL) {
  if (!is.null(seed)) { set.seed(seed) }
  XB = as.numeric(X[, true_id, drop = FALSE] %*% true_beta)
  var_epsilon = var(XB) * (1 - hsq) / hsq
  y = XB + rnorm(length(XB), mean = 0, sd = sqrt(var_epsilon))
  list(y = y, XB = XB, var_epsilon = var_epsilon)
}
get_sumstats_fast = function(X, y) {
  n = nrow(X)
  corr = apply(X, 2, function(x) cor(x, y))
  stats = corr * sqrt((n - 2) / (1 - corr^2))
  pval = 2 * pt(abs(stats), df = n - 2, lower.tail = FALSE)
  data.frame(variable = colnames(X), corr = as.numeric(corr), stats = as.numeric(stats),
             df = n - 2, pval = as.numeric(pval), n = n)
}
simulate_split = function(n, p, true_id, true_beta, hsq, seed, make_sumstats = FALSE) {
  X = gen_X(n = n, p = p, seed = seed)
  y_obj = gen_y_gaussian(X = X, true_id = true_id, true_beta = true_beta, hsq = hsq, seed = seed + 10000)
  out = list(X = X, y = y_obj$y, XB = y_obj$XB, var_epsilon = y_obj$var_epsilon,
             true_id = true_id, true_beta = true_beta, hsq = hsq)
  if (make_sumstats) { out$sumstats = get_sumstats_fast(X, y_obj$y) }
  out
}
simulate_cohort = function(cohort_name, true_beta, true_id, p, hsq, n_train, seed_train,
                           n_validation = NULL, seed_validation = NULL, n_testing = NULL, seed_testing = NULL) {
  cohort = list()
  cohort$train = simulate_split(n_train, p, true_id, true_beta, hsq, seed_train, make_sumstats = TRUE)
  if (!is.null(n_validation)) { cohort$validation = simulate_split(n_validation, p, true_id, true_beta, hsq, seed_validation) }
  if (!is.null(n_testing)) { cohort$testing = simulate_split(n_testing, p, true_id, true_beta, hsq, seed_testing) }
  cohort$cohort_name = cohort_name
  cohort
}

# -- 2. Data_BRIERfull -----------------------------------------------------------------------
data = list()
data$target = simulate_cohort("target", target_beta, true_id, p, hsq, 150, 1, 150, 2, 150, 3)
data$external1 = simulate_cohort("external1", ext1_beta, true_id, p, hsq, 300, 10, 150, 11)
data$external2 = simulate_cohort("external2", ext2_beta, true_id, p, hsq, 300, 20, 150, 21)
data$external3 = simulate_cohort("external3", ext3_beta, true_id, p, hsq, 300, 30, 150, 31)
Data_BRIERfull = data
save(Data_BRIERfull, file = file.path(OUT, "Data_BRIERfull.rda"))

# -- 3. Data_BRIERi: each external's lasso through BRIERi at eta = 0 --------------------------------
ext_i = function(e) {
  fit = BRIERi(e$train$X, e$train$y, eta.list = 0, beta.external = rep(0, ncol(e$train$X) + 1))
  fit = BRIERi.selection(fit, criteria = "gaussian.mspe", X.val = e$validation$X, y.val = e$validation$y)
  coef(fit)
}
d = Data_BRIERfull
B = cbind(ext_i(d$external1), ext_i(d$external2), ext_i(d$external3))
Data_BRIERi = list(target = d$target, beta.external = B)
save(Data_BRIERi, file = file.path(OUT, "Data_BRIERi.rda"))

# -- 4. Data_BRIERs: each external's lasso through BRIERs at eta = 0 ------------------------------
ext_s = function(e) {
  ld = calLD(e$train$X, e$train$sumstats)
  Xv = standardize_X(e$validation$X)[[1]]
  yv = standardize_X(as.matrix(e$validation$y))[[1]]
  fit = BRIERs(e$train$sumstats, ld$XtX, eta.list = 0, beta.external = rep(0, ncol(e$train$X)))
  fit = BRIERs.selection(fit, criteria = "gaussian.mspe", X.val = Xv, y.val = yv)
  coef(fit)
}
B = cbind(ext_s(d$external1), ext_s(d$external2), ext_s(d$external3))
Data_BRIERs = list(target = d$target, beta.external = B)
save(Data_BRIERs, file = file.path(OUT, "Data_BRIERs.rda"))
cat(sprintf("wrote %s\n", paste(file.path(OUT, c("Data_BRIERfull.rda", "Data_BRIERi.rda", "Data_BRIERs.rda")), collapse = ", ")))

# -- 5. compare with the shipped data ----------------------------------------------------------
if (!is.na(CMP)) {
  for (n in c("Data_BRIERfull", "Data_BRIERi", "Data_BRIERs")) {
    e = new.env(); load(file.path(CMP, paste0(n, ".rda")), envir = e)
    old = get(n, e); new = get(n)
    cat(sprintf("%-15s identical %s", n, identical(old, new)))
    if (!is.null(old$beta.external)) {
      bo = as.matrix(old$beta.external); bn = as.matrix(new$beta.external)
      cat(sprintf(" | target identical %s | beta.external dims %s vs %s, max |diff| %.3g, supports equal %s",
                  identical(old$target, new$target), paste(dim(bo), collapse = "x"), paste(dim(bn), collapse = "x"),
                  if (all(dim(bo) == dim(bn))) max(abs(bo - bn)) else NA,
                  if (all(dim(bo) == dim(bn))) identical(bo != 0, bn != 0) else NA))
    }
    cat("\n")
  }
}
