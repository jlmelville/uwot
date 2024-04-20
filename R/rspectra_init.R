rspectra_is_installed <- function() {
  is_installed("RSpectra")
}

rspectra_eigs_asym <- function(L, ndim) {
  res <- NULL
  suppressWarnings(res <- tryCatch(
    RSpectra::eigs(
      L,
      k = ndim + 1,
      which = "LR",
      opt = list(tol = 1e-4)
    ),
    error = function(c) {
      NULL
    }
  ))
  res
}

rspectra_eigs_sym <- function(L, ndim, verbose = FALSE) {
  k <- ndim + 1
  opt <- list(tol = 1e-4)
  suppressWarnings(res <-
    tryCatch(
      RSpectra::eigs_sym(L, k = k, which = "SM", opt = opt),
      error = function(c) {
        tsmessage("RSpectra calculation failed, retrying with shifted")
      }
    ))
  if (is.null(res) ||
    !is.list(res) ||
    !"vectors" %in% names(res) ||
    is.null(res$vectors) ||
    tryCatch(
      is.na(ncol(res$vectors)),
      error = function(e) {
        TRUE
      }
    ) ||
    ncol(res$vectors) < ndim) {
    suppressWarnings(res <- tryCatch(
      RSpectra::eigs_sym(
        L,
        k = k,
        which = "LM",
        sigma = 0,
        opt = opt
      ),
      error = function(c) {
        NULL
      }
    ))
  }
  res
}
