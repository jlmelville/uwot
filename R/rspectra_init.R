rspectra_eigs_asym <- function(L, ndim) {
  res <- NULL
  suppressWarnings(res <- tryCatch(
    RSpectra::eigs(
      L,
      k = ndim + 1,
      which = "LM",
      opt = list(tol = 1e-4)
    ),
    error = function(c) {
      NULL
    }
  ))
}

rspectra_eigs_sym <- function(L, ndim) {
  k <- ndim + 1
  opt <- list(tol = 1e-4)
  suppressWarnings(res <-
    tryCatch(
      RSpectra::eigs_sym(L, k = k, which = "SM", opt = opt),
      error = function(c) {
        NULL
      }
    ))
  if (is.null(res) || ncol(res$vectors) < ndim) {
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
