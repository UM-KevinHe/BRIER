#!/usr/bin/env Rscript
# brier_s_selection.R - select optimal (eta, lambda) from a cached BRIERs fit.
#
# Called by mcp/server.py as:
#   Rscript brier_s_selection.R <input.json> <output.json>
#
# UNLIKE brier_full_selection, BRIERs.selection accepts BOTH IC criteria
# (Cp, GIC, pseu.val) AND validation-set criteria. UNLIKE BRIERi.selection,
# the validation-set path requires the caller to PRE-STANDARDIZE X.val
# (and y.val for gaussian only). This is documented in llms.txt as a
# silent-failure trap.
#
# input.json: {
#   fit_id:    "brier_s_xxx",   # required
#   criteria:  "Cp" | "GIC" | "pseu.val"     # IC, no val data needed
#              | "gaussian.mspe" | "gaussian.rsq"
#              | "binomial.dev" | "binomial.mcfrsq" | "binomial.tjursq" | "binomial.auc"
#              | "poisson.dev",
#
#   # For validation-set criteria, ONE of:
#   X_val_expr:    "X.val",        # individual-level val data
#   y_val_expr:    "y.val",
#       OR
#   XtX_val_expr:  "ld_val$XtX",   # summary-level val data
#   sumstats_val_expr: "sumstats.val",
#   TN:            integer,         # total N at validation
#   h2:            numeric,         # heritability (can be 0 as fallback)
#
#   data_path:     ...              # required for either val mode
# }
#
# Caches the selection object and returns selection_id so brier_predict /
# brier_evaluate can use it.

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

.cache_root <- function() {
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

.IC_CRITERIA <- c("Cp", "GIC", "pseu.val")
.VAL_CRITERIA <- c("gaussian.mspe", "gaussian.rsq",
                    "binomial.dev", "binomial.mcfrsq",
                    "binomial.tjursq", "binomial.auc",
                    "poisson.dev")


.x_looks_standardized <- function(X) {
  if (!is.matrix(X) || ncol(X) < 5 || nrow(X) < 10) return(FALSE)
  k <- min(20, ncol(X))
  cols <- sample.int(ncol(X), k)
  means <- apply(X[, cols, drop = FALSE], 2, mean)
  sds <- apply(X[, cols, drop = FALSE], 2, sd)
  mean(abs(means) < 0.05) > 0.8 && mean(abs(sds - 1) < 0.1) > 0.8
}


args <- commandArgs(trailingOnly = TRUE)
io <- read_input(args)

result <- tryCatch({
  inp <- io$input

  if (is.null(inp$criteria) || !nzchar(inp$criteria)) {
    stop("criteria is required", call. = FALSE)
  }
  if (is.null(inp$fit_id) || !nzchar(inp$fit_id)) {
    stop("fit_id is required", call. = FALSE)
  }

  is_ic <- inp$criteria %in% .IC_CRITERIA
  is_val <- inp$criteria %in% .VAL_CRITERIA
  if (!is_ic && !is_val) {
    stop(sprintf(paste(
      "Unknown criteria='%s'. Valid IC criteria: %s. Valid",
      "validation-set criteria: %s."),
      inp$criteria,
      paste(.IC_CRITERIA, collapse = ", "),
      paste(.VAL_CRITERIA, collapse = ", ")
    ), call. = FALSE)
  }

  # IC criteria require TN (training sample size). Surface this before
  # BRIERs.selection does, with a clearer error.
  if (is_ic && is.null(inp$TN)) {
    stop(sprintf(paste(
      "criteria='%s' is an IC-based criterion and requires TN (training",
      "sample size, the n used to fit the model). Pass TN as an integer."
    ), inp$criteria), call. = FALSE)
  }

  fit_path <- file.path(.cache_root(), paste0(inp$fit_id, ".rds"))
  if (!file.exists(fit_path)) {
    stop(sprintf(
      "Fitted object not found at %s. Refit with brier_s.",
      fit_path
    ), call. = FALSE)
  }
  cached <- readRDS(fit_path)
  fit <- cached$fit
  meta <- cached$meta

  selection_args <- list(object = fit, criteria = inp$criteria)
  std_warning <- NULL

  # TN (training sample size) is required by IC criteria like Cp.
  # It's also accepted for some validation paths.
  if (!is.null(inp$TN)) {
    selection_args$TN <- as.numeric(inp$TN)
  }
  if (!is.null(inp$h2)) {
    selection_args$h2 <- as.numeric(inp$h2)
  }

  if (is_val) {
    has_individual <- !is.null(inp$X_val_expr) && !is.null(inp$y_val_expr)
    has_summary <- !is.null(inp$XtX_val_expr) && !is.null(inp$sumstats_val_expr)

    if (!has_individual && !has_summary) {
      stop(sprintf(paste(
        "criteria='%s' requires validation data. Provide EITHER",
        "X_val_expr + y_val_expr (individual-level), OR",
        "XtX_val_expr + sumstats_val_expr + TN + h2 (summary-level)."
      ), inp$criteria), call. = FALSE)
    }

    if (is.null(inp$data_paths) && is.null(inp$data_path)) {
      stop(paste("data_paths (or data_path) is required for",
                  "validation-set criteria"),
           call. = FALSE)
    }

    resolved_paths <- resolve_data_paths_input(inp)
    env <- load_data_files(resolved_paths)

    if (has_individual) {
      X_val <- safe_eval(inp$X_val_expr, env)
      y_val <- safe_eval(inp$y_val_expr, env)
      if (!is.matrix(X_val)) X_val <- as.matrix(X_val)
      if (is.matrix(y_val) && ncol(y_val) == 1) y_val <- as.vector(y_val)
      selection_args$X.val <- X_val
      selection_args$y.val <- y_val

      if (!.x_looks_standardized(X_val)) {
        std_warning <- paste(
          "Heuristic: X.val does NOT appear to be column-standardized",
          "(column means not near 0, sds not near 1). BRIERs returned",
          "coefficients on the standardized scale; passing un-standardized",
          "X.val to selection produces meaningless predictions. Pass",
          "X.val <- standardize_X(X)[['standardized']] instead."
        )
      }
    } else {
      XtX_val <- safe_eval(inp$XtX_val_expr, env)
      sumstats_val <- safe_eval(inp$sumstats_val_expr, env)
      selection_args$XtX <- XtX_val
      selection_args$sumstats <- sumstats_val
      # TN / h2 already set above if provided.
    }
  }

  sel <- do.call(BRIER::BRIERs.selection, selection_args)

  # Same return-shape conventions as BRIERi.selection / BRIERfull.selection.
  selected_eta <- if (is.list(sel$eta.min)) {
    lapply(sel$eta.min, function(x) as.numeric(unname(x)))
  } else {
    as.numeric(unname(sel$eta.min))
  }
  selected_lambda <- as.numeric(unname(sel$lambda.min))
  selected_metric <- tryCatch({
    idx <- sel$eta.min.index
    if (!is.null(idx) && !is.null(sel$eta.lambda$measure.min)) {
      as.numeric(sel$eta.lambda$measure.min[idx])
    } else {
      NULL
    }
  }, error = function(e) NULL)

  selection_id <- paste0(
    "brier_s_sel_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    paste(sample(c(0:9, letters), 6, replace = TRUE), collapse = "")
  )
  selection_path <- file.path(.cache_root(), paste0(selection_id, ".rds"))
  saveRDS(
    list(
      selection = sel,
      source_fit_id = inp$fit_id,
      criteria = inp$criteria
    ),
    file = selection_path
  )

  eta_grid_values <- tryCatch({
    el <- sel$eta.lambda
    eta_cols <- grep("^eta_\\d+$", colnames(el), value = TRUE)
    if (length(eta_cols) > 0) {
      as.numeric(unlist(el[, eta_cols, drop = FALSE]))
    } else {
      NULL
    }
  }, error = function(e) NULL)

  out <- list(
    status = "ok",
    fit_id = inp$fit_id,
    selection_id = selection_id,
    selection_path = selection_path,
    criteria = inp$criteria,
    criteria_mode = if (is_ic) "ic" else "validation",
    selected_eta = selected_eta,
    selected_lambda = selected_lambda,
    selected_metric = selected_metric,
    eta_grid_values = eta_grid_values
  )

  if (!is.null(std_warning)) {
    out$`_notice_x_val_not_standardized` <- std_warning
  }

  # Family-criterion mismatch (same pattern as the other selection tools).
  if (!is.null(meta$family) && meta$family == "binomial" &&
      grepl("^gaussian\\.", inp$criteria)) {
    out$`_notice_family_criterion_mismatch` <- paste(
      "Fit used family='binomial' but selection criterion is",
      sprintf("'%s' (a gaussian metric).", inp$criteria),
      "Use 'binomial.dev', 'binomial.mcfrsq', 'binomial.tjursq', or",
      "'binomial.auc' instead."
    )
  }

  # Boundary-saturation warning. BRIERs often shows high-eta optima
  # (pretrained external coefficients carry signal worth borrowing
  # heavily). If the selected eta is at or near the top of the fitted
  # grid, the true optimum may lie beyond.
  eta_grid_used <- tryCatch({
    if (!is.null(fit$eta.list) && is.list(fit$eta.list) &&
        length(fit$eta.list) >= 1L) {
      as.numeric(fit$eta.list[[1]])
    } else if (!is.null(fit$eta.list)) {
      as.numeric(fit$eta.list)
    } else NULL
  }, error = function(e) NULL)
  if (!is.null(eta_grid_used) && length(eta_grid_used) >= 3L) {
    max_eta <- max(eta_grid_used, na.rm = TRUE)
    sel_eta_scalar <- if (is.list(selected_eta)) {
      max(unlist(selected_eta), na.rm = TRUE)
    } else if (length(selected_eta) > 0L) {
      max(as.numeric(selected_eta), na.rm = TRUE)
    } else NA_real_
    # Flag if the selected eta is at the maximum OR within the top
    # quintile of the grid (heuristic).
    sorted_eta <- sort(eta_grid_used)
    top_quintile_threshold <- sorted_eta[
      max(1L, ceiling(0.8 * length(sorted_eta)))
    ]
    if (!is.na(sel_eta_scalar) && sel_eta_scalar >= top_quintile_threshold) {
      out$`_notice_eta_at_boundary` <- sprintf(paste(
        "Selected eta = %.4f is at or near the top of the fitted grid",
        "(max = %g). The true optimum may lie beyond. BRIERs often",
        "shows high-eta optima because pretrained external coefficients",
        "encode useful signal. To check, refit with an extended grid,",
        "e.g. eta_list = [0, 1, 5, 10, 25, 50, 100]."
      ), sel_eta_scalar, max_eta)
    }
  }

  out
}, error = function(err) {
  make_error(
    msg = conditionMessage(err),
    where = "brier_s_selection.R",
    class = class(err)[1]
  )
})

write_output(result, io$output_path)
