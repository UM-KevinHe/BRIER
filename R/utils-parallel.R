# Errors raised inside parallel::mclapply workers do NOT propagate: mclapply returns a
# "try-error" object in that slot (with only a warning), and a worker killed by the system
# (for example out of memory) returns NULL. Before this check, a fit whose every eta failed
# (an invalid penalty, a singular input) "succeeded", and the failure surfaced later as an
# unrelated error in *.selection ("$ operator is invalid for atomic vectors"). The serial
# path (lapply) already stops on the first error; this makes the parallel path do the same.
.stop_on_worker_errors <- function(res, eta.grid.list) {
  bad <- which(vapply(
    res,
    function(r) is.null(r) || inherits(r, "try-error"),
    logical(1)
  ))
  if (length(bad) == 0L) {
    return(invisible(res))
  }
  first <- res[[bad[1L]]]
  why <- if (is.null(first)) {
    "the worker returned no result (it may have been killed, for example out of memory)"
  } else {
    trimws(sub("^Error[^:]*:\\s*", "", as.character(first)))
  }
  etas <- vapply(
    eta.grid.list[bad],
    function(e) paste(format(signif(unlist(e), 4)), collapse = ", "),
    character(1)
  )
  stop(
    sprintf(
      paste0(
        "fitting failed in the parallel workers (%d of %d eta value(s) returned an error, ",
        "eta = %s; with mclapply prescheduling one failure marks every eta on the same core): %s"
      ),
      length(bad),
      length(res),
      paste(head(etas, 5L), collapse = "; "),
      why
    ),
    call. = FALSE
  )
}
