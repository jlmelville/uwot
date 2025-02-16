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
      opts = list(tol = 1e-4)
    ),
    error = function(c) {
      NULL
    }
  ))
  res
}

rspectra_eigs_sym <- function(L, ndim, verbose = FALSE, ...) {
  k <- ndim + 1
  opt <- lmerge(list(tol = 1e-4), list(...))

  suppressWarnings(res <-
    tryCatch(
      RSpectra::eigs_sym(L, k = k, which = "SM", opts = opt),
      error = function(c) {
       NULL
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
    tsmessage("RSpectra calculation failed, retrying with shifted")
    if ("initvec" %in% names(opt)) {
      opt$initvec <- 1 / opt$initvec
    }
    suppressWarnings(res <- tryCatch(
      RSpectra::eigs_sym(
        L,
        k = k,
        which = "LM",
        sigma = 1e-9,
        opts = opt
      ),
      error = function(c) {
        tsmessage("RSpectra shifted calculation also failed")
        NULL
      }
    ))
  }
  res
}

rspectra_eigs_shift_sym <- function(L, ndim, verbose = FALSE, ...) {
  k <- ndim + 1
  opts <- lmerge(list(tol = 1e-4), list(...))
  suppressWarnings(res <- tryCatch(
    RSpectra::eigs_sym(
      L,
      k = k,
      which = "LM",
      opts = opts
    ),
    error = function(c) {
      NULL
    }
  ))
  res
}
