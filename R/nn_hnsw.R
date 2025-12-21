hnsw_nn <- function(X,
                    k = 10,
                    metric = "euclidean",
                    M = 16,
                    ef_construction = 200,
                    ef = 10,
                    n_threads = NULL,
                    n_build_threads = NULL,
                    ret_index = FALSE,
                    verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (is.null(n_build_threads)) {
    n_build_threads <- n_threads
  }
  ann <-
    hnsw_build(X,
      metric = metric,
      M = M,
      ef_construction = ef_construction,
      n_build_threads = n_build_threads,
      verbose = verbose
    )
  res <-
    hnsw_search(X, k, ann, ef = ef, n_threads = n_threads, verbose = verbose)

  # We actually use the L2 HNSW metric so we need to convert here
  # (also umap_transform must do this)
  if (metric == "euclidean") {
    res$dist <- sqrt(res$dist)
  }

  res <- list(idx = res$idx, dist = res$dist)
  if (ret_index) {
    res$index <- ann
  }
  res
}

hnsw_build <- function(X, metric, M, ef_construction, n_build_threads, verbose) {
  hnsw_distance <- metric
  if (metric == "correlation") {
    tsmessage("Annoy build: subtracting row means for correlation")
    X <- sweep(X, 1, rowMeans(X))
    hnsw_distance <- "cosine"
  }
  # To avoid issues with whether a dedicated Euclidean class exists in RcppHNSW
  # we will always use L2 and manually process the distances when we are done
  if (metric == "euclidean") {
    hnsw_distance <- "l2"
  }

  index <-
    RcppHNSW::hnsw_build(
      X,
      distance = hnsw_distance,
      M = M,
      ef = ef_construction,
      verbose = verbose,
      n_threads = n_build_threads
    )
  list(
    ann = index,
    type = "hnswv1",
    metric = metric,
    ndim = ncol(X)
  )
}

# called by hnsw_nn when building a model, and by umap_transform directly
hnsw_search <-
  function(X,
           k,
           ann,
           ef,
           n_threads = NULL,
           verbose = FALSE) {
    if (is.null(n_threads)) {
      n_threads <- default_num_threads()
    }
    if (ann$metric == "correlation") {
      tsmessage("HNSW search: subtracting row means for correlation")
      X <- sweep(X, 1, rowMeans(X))
    }

    res <- RcppHNSW::hnsw_search(
      X = X,
      k = k,
      ann = ann$ann,
      ef = ef,
      n_threads = n_threads,
      verbose = verbose
    )

    res
  }

hnsw_load <- function(name, ndim, filename) {
  class_name <- switch(
    name,
    cosine = RcppHNSW::HnswCosine,
    euclidean = RcppHNSW::HnswL2,
    correlation = RcppHNSW::HnswCosine,
    stop("BUG: unknown HNSW metric '", name, "'")
  )
  methods::new(class_name, ndim, filename)
}

is_ok_hnsw_metric <- function(metric) {
  hnsw_metrics <- c("euclidean", "cosine", "correlation")
  metric %in% hnsw_metrics
}
