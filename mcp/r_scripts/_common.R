#!/usr/bin/env Rscript
# _common.R - shared utilities sourced by every BRIER MCP dispatcher script.
#
# What lives here:
#   * read_input(args)        - parse argv positions [1] (input.json path)
#                               and [2] (output.json path), load input JSON.
#   * write_output(result, p) - serialize a result list to output JSON
#                               with consistent options across all scripts.
#   * load_data_file(path)    - load .rda / .RData / .rds into a fresh env;
#                               returns the env (use `as.list(env)` to peek).
#   * safe_eval(expr_str, env)- eval an R expression string from MCP input,
#                               with belt-and-suspenders deny-list check
#                               matching server.py's pre-flight filter.
#   * make_error(msg, where, class)
#                             - construct a standard error payload.
#
# What does NOT live here:
#   * Statistical logic. Each dispatcher calls BRIER:: directly.
#   * Tool-specific shaping. Each dispatcher returns its own result list.

suppressPackageStartupMessages({
  library(jsonlite)
})


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------

read_input <- function(args) {
  if (length(args) < 2) {
    stop(
      "Usage: Rscript <script>.R <input.json> <output.json>",
      call. = FALSE
    )
  }
  list(
    input  = fromJSON(args[1], simplifyVector = FALSE),
    output_path = args[2]
  )
}

write_output <- function(result, path) {
  # Consistent JSON output across all scripts:
  #   * auto_unbox: scalars come out as JSON scalars, not 1-element arrays.
  #   * matrix = "rowmajor": match jsonlite default for matrices.
  #   * na, null: explicit "null" for missing values.
  #   * pretty: readable in audit logs and easier to diff in tests.
  writeLines(
    toJSON(
      result,
      auto_unbox = TRUE,
      matrix = "rowmajor",
      na = "null",
      null = "null",
      pretty = TRUE
    ),
    con = path
  )
}


# --------------------------------------------------------------------------
# Data file loading
# --------------------------------------------------------------------------

load_data_file <- function(path) {
  # Returns a fresh environment containing the loaded object(s).
  # For .rda / .RData: multi-object load() into env.
  # For .rds: single-object readRDS, named after the file basename.
  if (is.null(path) || !nzchar(path)) {
    stop("data_path is required", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  e <- new.env()

  if (ext %in% c("rda", "rdata")) {
    load(path, envir = e)
  } else if (ext == "rds") {
    obj_name <- tools::file_path_sans_ext(basename(path))
    assign(obj_name, readRDS(path), envir = e)
  } else {
    stop(
      sprintf(
        "Unsupported file extension: .%s (supported: .rda .RData .rds)",
        ext
      ),
      call. = FALSE
    )
  }

  e
}


load_data_files <- function(paths) {
  # v0.11: multi-file loading with basename wrapping.
  #
  # Each file's contents are wrapped under a top-level variable named
  # after the file basename, regardless of the file format. This makes
  # expressions like "height_AFR$sumstats" resolve consistently whether
  # the source was .rds or .RData.
  #
  # Behavior:
  # - For .rds: the readRDS() result is assigned to a variable named
  #   after the file basename.
  # - For .RData / .rda with one internal variable: that variable's
  #   value is assigned to a variable named after the file basename
  #   (whatever the internal name was, it gets renamed).
  # - For .RData / .rda with multiple internal variables: a list of
  #   all variables is assigned to a variable named after the file
  #   basename. Each can be accessed as basename$internal_var_name.
  #
  # This is a backward-incompatible change for .RData files: previously
  # load() exposed internal variable names directly. Now they live
  # under the file basename. This was a deliberate trade for predictable
  # naming across file formats. The single-file load_data_file() is
  # retained for use by code that explicitly wants the legacy behavior.
  if (length(paths) == 0L) {
    stop("data_paths must contain at least one path", call. = FALSE)
  }

  e <- new.env()

  # Determine if this is single-file legacy mode (one .RData/.rda path).
  # In that case we ALSO expose the internal variables at the top level
  # so pre-v0.11 expressions like "X" / "y" / "beta.external" still
  # resolve. The basename-wrapped name is also assigned, so new code
  # that uses "basename$X" continues to work too.
  single_legacy_rdata <- (length(paths) == 1L &&
                          tolower(tools::file_ext(paths[1])) %in%
                            c("rda", "rdata"))

  for (path in paths) {
    if (is.null(path) || !nzchar(path)) {
      stop("encountered an empty path in data_paths", call. = FALSE)
    }
    if (!file.exists(path)) {
      stop(sprintf("File not found: %s", path), call. = FALSE)
    }
    ext <- tolower(tools::file_ext(path))
    basename_var <- tools::file_path_sans_ext(basename(path))

    if (ext == "rds") {
      assign(basename_var, readRDS(path), envir = e)
    } else if (ext %in% c("rda", "rdata")) {
      tmp <- new.env()
      load(path, envir = tmp)
      vars <- ls(tmp)
      if (length(vars) == 1L) {
        # Single internal variable -> alias under basename
        assign(basename_var, get(vars[1], envir = tmp), envir = e)
        if (single_legacy_rdata) {
          # Backward compat: also expose the internal variable at the
          # top level under its original name.
          assign(vars[1], get(vars[1], envir = tmp), envir = e)
        }
      } else {
        # Multiple internal variables -> wrap as list under basename
        assign(basename_var, mget(vars, envir = tmp), envir = e)
        if (single_legacy_rdata) {
          # Backward compat: ALSO expose each variable at the top level.
          for (v in vars) {
            assign(v, get(v, envir = tmp), envir = e)
          }
        }
      }
    } else {
      stop(sprintf("Unsupported file extension: .%s", ext), call. = FALSE)
    }
  }

  e
}


resolve_data_paths_input <- function(inp) {
  # Backward-compat shim. Accept either `data_paths` (preferred, list)
  # or `data_path` (legacy, single string) from the input JSON.
  # Returns a character vector of paths.
  if (!is.null(inp$data_paths) && length(inp$data_paths) > 0L) {
    return(as.character(inp$data_paths))
  }
  if (!is.null(inp$data_path) && nzchar(inp$data_path)) {
    return(as.character(inp$data_path))
  }
  stop("either data_paths or data_path is required", call. = FALSE)
}


extract_prep_session_ids <- function(paths) {
  # v0.13: when a fit reads a .rds that was written by prep_data's
  # `persist` op, the file contains a `.prep_meta` list with the prep
  # session id. Scan every input path and collect any prep_session_ids
  # found. Returns a character vector (possibly empty).
  ids <- character(0)
  for (p in paths) {
    if (is.null(p) || !nzchar(p)) next
    ext <- tolower(tools::file_ext(p))
    if (ext != "rds") next  # only persisted-prep files use .rds
    if (!file.exists(p)) next
    obj <- tryCatch(readRDS(p), error = function(e) NULL)
    if (is.list(obj) && !is.data.frame(obj) &&
        !is.null(obj$.prep_meta) &&
        is.list(obj$.prep_meta) &&
        !is.null(obj$.prep_meta$prep_session_id)) {
      ids <- c(ids, as.character(obj$.prep_meta$prep_session_id))
    }
  }
  unique(ids)
}


# --------------------------------------------------------------------------
# Expression evaluation (with deny-list)
# --------------------------------------------------------------------------

# Mirror of server.py:DENY_PATTERNS. Belt-and-suspenders: server.py runs
# this check before writing JSON, but R-side enforcement protects against
# any path that bypasses the Python pre-flight (e.g. direct Rscript
# invocation during development).
# Mirror of server.py:DENY_PATTERNS. Belt-and-suspenders: server.py runs
# this check before writing JSON, but R-side enforcement protects against
# any path that bypasses the Python pre-flight (e.g. direct Rscript
# invocation during development).
#
# `::` is not on the deny list; instead it is whitelist-checked below so
# `BRIER::standardize_X(X)` and similar safe namespace calls are allowed.
# `:::` (non-exported access) stays denied.
.DENY_PATTERNS <- c(
  "system(", "system2(", "shell(", "shell.exec(",
  "unlink(", "file.remove(", "file.rename(",
  "file.create(", "file.copy(",
  "eval(", "parse(", "source(",
  "Sys.setenv(", "Sys.unsetenv(",
  "do.call(",
  ":::",
  "`", ";"
)

# Safe namespace prefixes mirroring server.py:SAFE_NAMESPACE_PREFIXES.
# An expression containing `::` is rejected unless every occurrence is
# preceded by one of these prefixes.
.SAFE_NAMESPACE_PREFIXES <- c(
  "BRIER::", "base::", "stats::", "utils::", "Matrix::"
)

.expr_uses_only_safe_namespaces <- function(expr_str) {
  if (grepl(":::", expr_str, fixed = TRUE)) return(FALSE)
  if (!grepl("::", expr_str, fixed = TRUE)) return(TRUE)
  remaining <- expr_str
  while (grepl("::", remaining, fixed = TRUE)) {
    idx <- regexpr("::", remaining, fixed = TRUE)
    pos <- as.integer(idx)
    # Walk back from `pos` to find the start of the identifier
    i <- pos - 1L
    while (i >= 1L) {
      ch <- substr(remaining, i, i)
      if (grepl("[A-Za-z0-9_.]", ch)) {
        i <- i - 1L
      } else {
        break
      }
    }
    ns <- substr(remaining, i + 1L, pos + 1L)  # includes the "::"
    if (!(ns %in% .SAFE_NAMESPACE_PREFIXES)) return(FALSE)
    remaining <- substr(remaining, pos + 2L, nchar(remaining))
  }
  TRUE
}

safe_eval <- function(expr_str, env) {
  # Returns NULL if expr_str is NULL / empty; otherwise evaluates inside `env`.
  # Throws if the expression matches any deny-list pattern OR uses an
  # unwhitelisted namespace.
  if (is.null(expr_str) || !is.character(expr_str) || !nzchar(expr_str)) {
    return(NULL)
  }
  for (pat in .DENY_PATTERNS) {
    if (grepl(pat, expr_str, fixed = TRUE)) {
      stop(
        sprintf(
          "Refusing to evaluate expression: contains disallowed pattern %s",
          shQuote(pat)
        ),
        call. = FALSE
      )
    }
  }
  if (!.expr_uses_only_safe_namespaces(expr_str)) {
    stop(
      sprintf(
        "Refusing to evaluate expression: contains '::' but not from an allowed namespace (%s)",
        paste(.SAFE_NAMESPACE_PREFIXES, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  eval(parse(text = expr_str), envir = env)
}


# --------------------------------------------------------------------------
# Error construction
# --------------------------------------------------------------------------

make_error <- function(msg, where, class = "Error") {
  list(
    status = "error",
    message = msg,
    class = class,
    where = where
  )
}
