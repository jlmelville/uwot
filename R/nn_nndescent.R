nndescent_nn <- function(X,
                         k = 10,
                         metric = "euclidean",
                         nn_args = list(),
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

  if (!ret_index) {
    nn_knn_args <- get_nndescent_knn_args(nn_args)
    nn_knn_args <- lmerge(
      nn_knn_args,
      list(
        data = X,
        k = k,
        metric = metric,
        n_threads = n_build_threads,
        verbose = verbose
      )
    )
    return(do.call(rnndescent::rnnd_knn, nn_knn_args))
  }

  ann <- nndescent_build(
    X,
    k,
    metric,
    nn_args = nn_args,
    n_build_threads = n_build_threads,
    verbose = verbose
  )
  res <-
    list(
      idx = ann$ann$graph$idx,
      dist = ann$ann$graph$dist,
      index = ann
    )
  res$index$ann$ann$graph <- NULL
  res
}

nndescent_build <- function(X,
                            k,
                            metric,
                            nn_args = list(),
                            n_build_threads = NULL,
                            verbose = FALSE) {
  nn_build_args <- get_nndescent_build_args(nn_args)
  nn_build_args <- lmerge(
    nn_build_args,
    list(
      data = X,
      k = k,
      metric = metric,
      n_threads = n_build_threads,
      verbose = verbose
    )
  )

  index <- do.call(rnndescent::rnnd_build, nn_build_args)
  list(
    ann = index,
    type = "nndescentv1",
    metric = metric,
    ndim = ncol(X)
  )
}


nndescent_search <- function(X,
                             k,
                             ann,
                             nn_args = list(),
                             n_threads = NULL,
                             verbose = FALSE) {
  nn_query_args <- get_nndescent_query_args(nn_args)
  nn_query_args <- lmerge(
    nn_query_args,
    list(
      index = ann$ann,
      query = X,
      k = k,
      n_threads = n_threads,
      verbose = verbose
    )
  )

  do.call(rnndescent::rnnd_query, nn_query_args)
}

get_nndescent_knn_args <- function(nn_args) {
  nn_knn_args <- list()
  nnd_knn_names <- c(
    "use_alt_metric",
    "init",
    "n_trees",
    "leaf_size",
    "max_tree_depth",
    "margin",
    "n_iters",
    "delta",
    "max_candidates",
    "weight_by_degree",
    "low_memory"
  )
  for (name in nnd_knn_names) {
    if (name %in% names(nn_args)) {
      nn_knn_args[[name]] <- nn_args[[name]]
    }
  }
  nn_knn_args
}

get_nndescent_build_args <- function(nn_args) {
  # prune_reverse should probably always be TRUE
  nn_build_args <- list(prune_reverse = TRUE)
  nnd_build_names <- c(
    "use_alt_metric",
    "init",
    "n_trees",
    "leaf_size",
    "max_tree_depth",
    "margin",
    "n_iters",
    "delta",
    "max_candidates",
    "weight_by_degree",
    "low_memory",
    "n_search_trees",
    "pruning_degree_multiplier",
    "diversify_prob",
    "prune_reverse"
  )
  for (name in nnd_build_names) {
    if (name %in% names(nn_args)) {
      nn_build_args[[name]] <- nn_args[[name]]
    }
  }
  nn_build_args
}

get_nndescent_query_args <- function(nn_args) {
  nn_query_args <- list()
  nnd_query_names <- c(
    "epsilon",
    "max_search_fraction"
  )
  for (name in nnd_query_names) {
    if (name %in% names(nn_args)) {
      nn_query_args[[name]] <- nn_args[[name]]
    }
  }
  nn_query_args
}
