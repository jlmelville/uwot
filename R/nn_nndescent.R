nndescent_nn <- function(X,
                         k = 10,
                         metric = "euclidean",
                         nn_args = list(),
                         n_threads = NULL,
                         ret_index = FALSE,
                         verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }

  if (!ret_index) {
    nn_args$data <- X
    nn_args$k <- k
    nn_args$metric <- metric
    nn_args$n_threads <- n_threads
    nn_args$verbose <- verbose
    return(do.call(rnndescent::rnnd_knn, nn_args))
  }

  ann <- nndescent_build(
    X,
    k,
    metric,
    nn_args = nn_args,
    n_threads = n_threads,
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
                            n_threads = NULL,
                            verbose = FALSE) {
  nn_args$data <- X
  nn_args$k <- k
  nn_args$metric <- metric
  nn_args$n_threads <- n_threads
  nn_args$verbose <- verbose
  index <- do.call(rnndescent::rnnd_build, nn_args)
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
  nn_args$index <- ann$ann
  nn_args$query <- X
  nn_args$k <- k
  nn_args$n_threads <- n_threads
  nn_args$verbose <- verbose
  do.call(rnndescent::rnnd_query, nn_args)
}
