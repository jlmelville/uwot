hnsw_nn <- function(X,
                    k = 10,
                    metric = "euclidean",
                    n_threads = NULL,
                    ret_index = FALSE,
                    verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  ann <-
    hnsw_build(X,
      metric = metric,
      n_threads = n_threads,
      verbose = verbose
    )
  res <-
    hnsw_search(X, k, ann, n_threads = n_threads, verbose = verbose)


  res <- list(idx = res$idx, dist = res$dist)
  if (ret_index) {
    res$index <- ann
  }
  res
}

hnsw_build <- function(X, metric, n_threads, verbose) {

  hnsw_distance <- metric
  if (metric == "correlation") {
    tsmessage("Annoy build: subtracting row means for correlation")
    X <- sweep(X, 1, rowMeans(X))
    hnsw_distance <- "cosine"
  }

  # FIXME: allow for different M and ef values
  index <-
    RcppHNSW::hnsw_build(
      X,
      distance = hnsw_distance,
      M = 16,
      ef = 200,
      verbose = verbose,
      n_threads = n_threads
    )
  list(
    ann = index,
    type = "hnswv1",
    metric = metric,
    ndim = ncol(X)
  )
}

hnsw_search <-
  function(X,
           k,
           ann,
           n_threads = NULL,
           verbose = FALSE) {
    if (is.null(n_threads)) {
      n_threads <- default_num_threads()
    }
    if (ann$metric == "correlation") {
      tsmessage("HNSW search: subtracting row means for correlation")
      X <- sweep(X, 1, rowMeans(X))
    }

    # FIXME: allow for different ef values
    res <- RcppHNSW::hnsw_search(
      X = X,
      k = k,
      ann = ann$ann,
      ef = 10,
      n_threads = n_threads,
      verbose = verbose
    )

    res
  }

hnsw_load <- function(name, ndim, filename) {
  class_name <- switch(
    name,
    cosine =RcppHNSW::HnswCosine,
    euclidean = RcppHNSW::HnswL2,
    correlation = RcppHNSW::HnswCosine,
    stop("BUG: unknown HNSW metric '", name, "'")
  )
  methods::new(class_name, ndim, filename)
}
