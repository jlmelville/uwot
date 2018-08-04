umap_transform <- function(X, model,
                           weighted = TRUE,
                           n_threads = max(1, RcppParallel::defaultNumThreads() / 2),
                           grain_size = 1, verbose = FALSE) {
  nn_index <- model$nn_index
  n_neighbors <- model$n_neighbors
  search_k <- model$search_k
  local_connectivity <- model$local_connectivity
  train_embedding <- model$embedding
  
  # TODO return scaling info from uwot and apply it
  scale <- FALSE
  # TODO: establish what kind of input is acceptable
  if (methods::is(X, "dist")) {
    n_vertices <- attr(X, "Size")
    tsmessage("Read ", n_vertices, " rows")
  }
  else if (methods::is(X, "sparseMatrix")) {
    n_vertices <- nrow(X)
    if (ncol(X) != n_vertices) {
      stop("Sparse matrices are only supported as distance matrices")
    }
    tsmessage("Read ", n_vertices, " rows of sparse distance matrix")
  }
  else {
    if (methods::is(X, "data.frame")) {
      indexes <- which(vapply(X, is.numeric, logical(1)))
      if (length(indexes) == 0) {
        stop("No numeric columns found")
      }
      X <- as.matrix(X[, indexes])
    }
    n_vertices <- nrow(X)
    tsmessage("Read ", n_vertices, " rows and found ", ncol(X),
              " numeric columns")
    X <- scale_input(X, scale_type = scale, verbose = verbose)
  }
  
  
  if (n_threads > 0) {
    RcppParallel::setThreadOptions(numThreads = n_threads)
  }
  
  nn <- annoy_search(X, k = n_neighbors, nn_index, search_k = search_k,
                     n_threads = n_threads, grain_size = grain_size,
                     verbose = verbose)
  if (weighted) {
    adjusted_local_connectivity <- max(0, local_connectivity - 1.0)
    graph <- smooth_knn(nn,
                        local_connectivity = adjusted_local_connectivity,
                        n_threads = n_threads,
                        grain_size = grain_size,
                        self_nbr = FALSE,
                        n_reference_vertices = nrow(train_embedding),
                        verbose = verbose)
  }
  if (weighted) {
    tsmessage("Initializing by weighted average of neighbor coordinates")
    embedding <- init_transform_cpp(train_embedding, nn$idx, graph)
    # embedding <- init_transform(train_embedding, nn$idx, graph)
  }
  else {
    tsmessage("Initializing by average of neighbor coordinates")
    embedding <- init_transform_av_cpp(train_embedding, nn$idx)
    # embedding <- init_transform(train_embedding, nn$idx)
  }
  
  tsmessage("Finished")
  embedding
}

init_transform <- function(train_embedding, nn_index, weights = NULL) {
  nr <- nrow(nn_index)
  nc <- ncol(train_embedding)
  nnbrs <- ncol(nn_index)
  
  embedding <- matrix(nrow = nr, ncol = nc)
  if (is.null(weights)) {
    for (i in 1:nr) {
      nbr_embedding <- train_embedding[nn_index[i, ], ]
      embedding[i, ] <- apply(nbr_embedding, 2, mean)
    }
  }
  else {
    offset <- nrow(train_embedding)
    for (i in 1:nr) {
      nbr_embedding <- train_embedding[nn_index[i, ], ]
      nbr_weights <- weights[nn_index[i, ], i]
      embedding[i, ] <- apply(nbr_embedding, 2,
                              function(x) { stats::weighted.mean(x, nbr_weights) })
    }
  }
  
  embedding
}
