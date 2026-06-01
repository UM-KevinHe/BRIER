#!/usr/bin/env Rscript
# brier_s.R - fit BRIERs() with summary-statistics target data.
#
# Called by mcp/server.py as:
#   Rscript brier_s.R <input.json> <output.json>
#
# Distinct from brier_i / brier_full:
#   * Target is summary statistics (corr + LD matrix), NOT individual X/y.
#   * beta.external has NO intercept row (asymmetric with BRIERi which
#     requires one); shape is p x M, not (p+1) x M.
#   * Returns coefficients on the STANDARDIZED scale.
#
# input.json: {
#   data_path:           "/path/to/data.rds",     # required
#   sumstats_expr:       "sumstats",              # required (must have $corr column)
#   beta_external_expr:  "beta.external",         # required (p x M)
#   family:              "gaussian" | ...,         # required
#
#   # LD matrix: PREFERRED to pass ld_id (auto-subset by $nz)
#   ld_id:               "ld_xxx",                # OR
#   XtX_expr:            "ld$XtX",                # explicit (caller subsets manually)
#
#   # Other args
#   multi_method:        "stacking" | "PCA" | "ind",
#   eta_list:            [...],
#   trace:               false
# }
#
# output.json: {
#   status: "ok",
#   fit_id, fit_path,
#   family, p, M_external, eta_list_used,
#   timing, ld_id_used (if any),
#   _notice_*, _followup_*
# } or {status: "error", ...}

.script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }
})()
source(file.path(.script_dir, "_common.R"))

suppressPackageStartupMessages({
  library(BRIER)
  library(Matrix)
})


.cache_root_fits <- function() {
  base <- Sys.getenv("XDG_CACHE_HOME", unset = NA)
  if (is.na(base) || !nzchar(base)) {
    base <- if (.Platform$OS.type == "windows") {
      Sys.getenv("LOCALAPPDATA",
                 unset = file.path(Sys.getenv("HOME"), "AppData", "Local"))
    } else {
      file.path(Sys.getenv("HOME"), ".cache")
    }
  }
  d <- file.path(base, "brier-mcp", "fits")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.cache_root_ld <- function() {
  base <- Sys.getenv("XDG_CACHE_HOME", unset = NA)
  if (is.na(base) || !nzchar(base)) {
    base <- if (.Platform$OS.type == "windows") {
      Sys.getenv("LOCALAPPDATA",
                 unset = file.path(Sys.getenv("HOME"), "AppData", "Local"))
    } else {
      file.path(Sys.getenv("HOME"), ".cache")
    }
  }
  file.path(base, "brier-mcp", "ld")
}

.generate_fit_id <- function() {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  suffix <- paste(sample(c(0:9, letters), 6, replace = TRUE), collapse = "")
  paste0("brier_s_", ts, "_", suffix)
}


args <- commandArgs(trailingOnly = TRUE)
io <- read_input(args)

result <- tryCatch({
  inp <- io$input

  if (is.null(inp$data_paths) && is.null(inp$data_path)) {
    stop("either data_paths or data_path is required", call. = FALSE)
  }
  if (is.null(inp$sumstats_expr)) {
    stop("sumstats_expr is required", call. = FALSE)
  }
  if (is.null(inp$beta_external_expr)) {
    stop("beta_external_expr is required", call. = FALSE)
  }

  family_was_supplied <- !is.null(inp$family) && nzchar(inp$family)
  family <- if (family_was_supplied) inp$family else "gaussian"

  multi_method <- if (!is.null(inp$multi_method) && nzchar(inp$multi_method)) {
    inp$multi_method
  } else {
    "stacking"
  }

  # Resolve XtX. Two paths: ld_id (preferred, auto-subset) or XtX_expr (raw).
  # v0.11: multi-file support via load_data_files()
  resolved_paths <- resolve_data_paths_input(inp)
  env <- load_data_files(resolved_paths)
  sumstats <- safe_eval(inp$sumstats_expr, env)
  beta_external <- safe_eval(inp$beta_external_expr, env)

  if (is.null(sumstats)) stop("sumstats_expr resolved to NULL", call. = FALSE)
  if (is.null(beta_external)) {
    stop("beta_external_expr resolved to NULL", call. = FALSE)
  }
  if (!is.matrix(beta_external)) beta_external <- as.matrix(beta_external)
  if (!"corr" %in% colnames(sumstats)) {
    stop(paste(
      "sumstats must have a 'corr' column. Build with p2cor(pval, n,",
      "sign = sign(stats)) if you only have p-values."
    ), call. = FALSE)
  }

  XtX <- NULL
  ld_id_used <- NULL
  used_ld_subset <- FALSE

  if (!is.null(inp$ld_id) && nzchar(inp$ld_id)) {
    ld_path <- file.path(.cache_root_ld(), paste0(inp$ld_id, ".rds"))
    if (!file.exists(ld_path)) {
      stop(sprintf(
        "LD object not found at %s. Re-run cal_ld to regenerate.",
        ld_path
      ), call. = FALSE)
    }
    ld <- readRDS(ld_path)
    XtX <- ld$XtX
    ld_id_used <- inp$ld_id

    # Auto-subset by $nz. This is the silent-failure trap from llms.txt
    # that we systematically remove via the ld_id workflow.
    if (!is.null(ld$nz)) {
      if (nrow(sumstats) > length(ld$nz)) {
        sumstats <- sumstats[ld$nz, , drop = FALSE]
        used_ld_subset <- TRUE
      }
      if (nrow(beta_external) > length(ld$nz)) {
        beta_external <- beta_external[ld$nz, , drop = FALSE]
        used_ld_subset <- TRUE
      }
    }
  } else if (!is.null(inp$XtX_expr) && nzchar(inp$XtX_expr)) {
    XtX <- safe_eval(inp$XtX_expr, env)
    if (is.null(XtX)) stop("XtX_expr resolved to NULL", call. = FALSE)
  } else {
    stop(paste(
      "Either ld_id (from a prior cal_ld call) or XtX_expr is required.",
      "Recommended: use ld_id so sumstats and beta.external get",
      "auto-subset by the LD matrix's $nz indices."
    ), call. = FALSE)
  }

  # Shape sanity: sumstats rows == XtX rows == beta.external rows.
  if (nrow(sumstats) != nrow(XtX) || nrow(beta_external) != nrow(XtX)) {
    stop(sprintf(paste(
      "Row count mismatch after any LD-subset: sumstats has %d rows,",
      "XtX has %d rows, beta.external has %d rows. All three must match."
    ), nrow(sumstats), nrow(XtX), nrow(beta_external)), call. = FALSE)
  }

  # BRIERs takes XtX as Matrix-like; pass through as-is (BRIERs internally
  # coerces to sparse).

  fit_args <- list(
    sumstats = sumstats,
    XtX = XtX,
    family = family,
    beta.external = beta_external,
    multi.method = multi_method,
    trace = isTRUE(inp$trace),
    parallel = FALSE,
    ncores = 1
  )
  if (!is.null(inp$eta_list)) {
    fit_args$eta.list <- as.numeric(unlist(inp$eta_list))
  } else {
    fit_args$eta.list <- c(0, exp(seq(log(0.1), log(10), length.out = 20)))
  }

  # M=1 auto-substitute. With one external, stacking and PCA collapse
  # to ind mathematically. Substitute silently for robustness.
  m_one_auto_ind_applied <- FALSE
  if (ncol(beta_external) == 1L &&
      (identical(multi_method, "stacking") ||
       identical(multi_method, "PCA"))) {
    fit_args$multi.method <- "ind"
    m_one_auto_ind_applied <- TRUE
  }

  t0 <- Sys.time()
  fit <- do.call(BRIER::BRIERs, fit_args)
  t1 <- Sys.time()
  fit_seconds <- as.numeric(difftime(t1, t0, units = "secs"))

  cache_dir <- .cache_root_fits()
  fit_id <- .generate_fit_id()
  fit_path <- file.path(cache_dir, paste0(fit_id, ".rds"))

  saveRDS(
    list(
      fit = fit,
      meta = list(
        family = family,
        multi_method = multi_method,
        data_path = inp$data_path,
        data_paths = resolved_paths,
        sumstats_expr = inp$sumstats_expr,
        beta_external_expr = inp$beta_external_expr,
        XtX_expr = inp$XtX_expr,
        ld_id_used = ld_id_used,
        eta_list = inp$eta_list,
        tool = "brier_s",
        prep_session_ids = extract_prep_session_ids(resolved_paths)
      )
    ),
    file = fit_path
  )

  eta_used <- tryCatch(fit$eta.list, error = function(e) NULL)
  eta_serialized <- if (is.null(eta_used)) {
    NULL
  } else if (is.list(eta_used) && length(eta_used) == 1) {
    as.numeric(eta_used[[1]])
  } else if (is.list(eta_used)) {
    lapply(eta_used, as.numeric)
  } else {
    as.numeric(eta_used)
  }

  M_external <- ncol(beta_external)

  out <- list(
    status = "ok",
    fit_id = fit_id,
    fit_path = fit_path,
    family = family,
    p = nrow(XtX),
    M_external = M_external,
    eta_list_used = eta_serialized,
    multi_method_used = multi_method,
    ld_id_used = ld_id_used,
    timing = list(fit_seconds = round(fit_seconds, 3))
  )

  if (!family_was_supplied) {
    out$`_notice_family_default` <- paste(
      "Family was not explicitly supplied; BRIERs used the gaussian default.",
      "If the outcome is binary or count, refit with family='binomial' or",
      "'poisson'."
    )
  }

  # Always-on standardization warning. This is THE big BRIERs pitfall.
  out$`_notice_brier_s_standardize` <- paste(
    "BRIERs returns coefficients on the STANDARDIZED scale. When you call",
    "brier_s_selection with a validation-set criterion (gaussian.mspe,",
    "binomial.dev, etc.), the X.val matrix you pass MUST be column-",
    "standardized (e.g., via standardize_X). Pass y.val standardized for",
    "family='gaussian' only; for family='binomial' or 'poisson' pass raw",
    "y.val. See llms.txt 'BRIERs() returns standardized coefficients'."
  )

  if (used_ld_subset) {
    out$`_notice_ld_subset_applied` <- paste(
      "ld_id was used; sumstats and beta.external were automatically",
      "subset by the LD matrix's $nz indices before fitting. The fit's",
      "p =", nrow(XtX), "reflects the retained variants, which may be",
      "smaller than the input sumstats / beta.external row counts."
    )
  }

  if (identical(multi_method, "ind") && M_external >= 5) {
    out$`_notice_multi_ind_slow` <- sprintf(
      paste(
        "multi.method='ind' tunes an independent eta per external model.",
        "With M=%d this grid is multiplicative; expect very slow fits.",
        "Consider multi.method='stacking'."
      ), M_external)
  }

  if (m_one_auto_ind_applied) {
    out$`_notice_m_one_auto_ind` <- paste(
      "Detected M=1 (a single external) with multi_method=",
      sprintf("'%s'.", multi_method),
      "With one external, stacking and PCA collapse mathematically to",
      "ind. Silently switched to multi_method='ind'; the result is",
      "identical. Pass multi_method='ind' explicitly to skip this notice."
    )
    out$multi_method_used <- "ind"
  }

  out$`_followup_offer_selection` <- paste(
    "This is the raw fit across the (eta, lambda) grid. To select the",
    sprintf("optimal hyperparameters, call brier_s_selection with fit_id='%s'", fit_id),
    "and one of: (a) an IC criterion like 'Cp', 'GIC', 'pseu.val';",
    "(b) a family-specific validation metric ('gaussian.mspe', etc.)",
    "plus X_val_expr / y_val_expr / data_path with STANDARDIZED X.val and",
    "(for gaussian only) standardized y.val."
  )

  out
}, error = function(err) {
  make_error(
    msg = conditionMessage(err),
    where = "brier_s.R",
    class = class(err)[1]
  )
})

write_output(result, io$output_path)
