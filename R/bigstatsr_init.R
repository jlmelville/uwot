bigstatsr_is_installed <- function() {
  is_installed("bigstatsr")
}

bigstatsr_scores <- function(X,
                             ncol,
                             center = TRUE,
                             ret_extra = FALSE,
                             ncores = 1,
                             verbose = FALSE) {
  res <- bigstatsr::big_randomSVD(
    X = bigstatsr::as_FBM(X),
    fun.scaling = bigstatsr::big_scale(center = center, scale = FALSE),
    k = ncol,
    ncores = ncores
  )
  if (verbose) {
    totalvar <- sum(apply(X, 2, stats::var))
    lambda <- sum((res$d^2) / (nrow(X) - 1))
    varex <- lambda / totalvar
    tsmessage(
      "PCA: ",
      ncol,
      " components explained ",
      formatC(varex * 100),
      "% variance"
    )
  }
  scores <- stats::predict(res)
  if (ret_extra) {
    list(
      scores = scores,
      rotation = res$v,
      center = res$center
    )
  } else {
    scores
  }
}
