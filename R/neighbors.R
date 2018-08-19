find_nn <- function(X, k, include_self = TRUE, method = "fnn",
                    metric = "euclidean",
                    n_trees = 50, search_k = 2 * k * n_trees,
                    n_threads = max(1, RcppParallel::defaultNumThreads() / 2),
                    grain_size = 1,
                    ret_index = FALSE,
                    verbose = FALSE) {
  if (methods::is(X, "dist")) {
    res <- dist_nn(X, k, include_self = include_self)
  }
  else if (methods::is(X, "sparseMatrix")) {
    # sparse distance matrix
    res <- sparse_nn(X, k, include_self = include_self)
  }
  else {
    # normal matrix
    if (method == "fnn") {
      res <- FNN_nn(X, k = k, include_self = include_self)
    }
    else {
      res <- annoy_nn(X,
        k = k, include_self = include_self,
        metric = metric,
        n_trees = n_trees, search_k = search_k,
        n_threads = n_threads,
        ret_index = ret_index,
        verbose = verbose
      )
    }
  }
  res
}

# n_trees - number of trees to build when constructing the index. The more trees
# specified, the larger the index, but the better the results. largeVis uses 10
# trees for datasets with N = 10,000 observations, 20 trees for datasets up to N
# = 1,000,000, 50 trees for N up to 5,000,000 and 100 trees otherwise
# search_k - the number of nodes to search during the neighbor retrieval. The
# larger k, the more accurate results, but the longer the search takes. Default
# is k * n_trees.
#' @importFrom methods new
annoy_nn <- function(X, k = 10, include_self = TRUE,
                     metric = "euclidean",
                     n_trees = 50, search_k = 2 * k * n_trees,
                     n_threads = max(1, RcppParallel::defaultNumThreads() / 2),
                     grain_size = 1,
                     ret_index = FALSE,
                     verbose = FALSE) {
  ann <- annoy_build(X,
    metric = metric, n_trees = n_trees,
    n_threads = n_threads, grain_size = grain_size,
    verbose = verbose
  )

  # Search index
  if (!include_self) {
    k <- k + 1
  }
  res <- annoy_search(X,
    k = k, ann = ann, search_k = search_k,
    n_threads = n_threads,
    grain_size = grain_size, verbose = verbose
  )
  idx <- res$idx
  dist <- res$dist
  if (!include_self) {
    idx <- idx[, -1]
    dist <- dist[, -1]
    k <- k - 1
  }

  res <- list(idx = idx, dist = dist)
  if (ret_index) {
    res$index <- ann
  }
  res
}

annoy_build <- function(X, metric = "euclidean", n_trees = 50,
                        n_threads =
                          max(1, RcppParallel::defaultNumThreads() / 2),
                        grain_size = 1, verbose = FALSE) {
  nr <- nrow(X)
  nc <- ncol(X)

  if (metric == "cosine") {
    ann <- methods::new(RcppAnnoy::AnnoyAngular, nc)
  }
  else if (metric == "manhattan") {
    ann <- methods::new(RcppAnnoy::AnnoyManhattan, nc)
  }
  else {
    ann <- methods::new(RcppAnnoy::AnnoyEuclidean, nc)
  }

  tsmessage("Building Annoy index with metric = ", metric)
  progress <- Progress$new(max = nr, display = verbose)

  # Add items
  for (i in 1:nr) {
    ann$addItem(i - 1, X[i, ])
    progress$increment()
  }

  # Build index
  ann$build(n_trees)

  ann
}

# Search a pre-built Annoy index for neighbors of X
annoy_search <- function(X, k = 10, ann,
                         search_k = 100 * k,
                         n_threads =
                           max(1, RcppParallel::defaultNumThreads() / 2),
                         grain_size = 1,
                         verbose = FALSE) {
  if (n_threads > 0) {
    index_file <- tempfile()
    ann$save(index_file)

    tsmessage("Searching Annoy index using ", pluralize("thread", n_threads))

    ann_class <- class(ann)
    if (endsWith(ann_class, "Angular")) {
      search_nn_func <- annoy_cosine_nns
    }
    else if (endsWith(ann_class, "Manhattan")) {
      search_nn_func <- annoy_manhattan_nns
    }
    else if (endsWith(ann_class, "Euclidean")) {
      search_nn_func <- annoy_euclidean_nns
    }
    else {
      stop("BUG: Unknown ANN class '", ann_class, "'")
    }

    res <- search_nn_func(index_file,
      X,
      k, search_k,
      grain_size = grain_size,
      verbose = verbose
    )
    idx <- res$item
    dist <- res$distance

    unlink(index_file)
  }
  else {
    tsmessage("Searching Annoy index")
    nr <- nrow(X)
    search_progress <- Progress$new(max = nr, display = verbose)
    idx <- matrix(nrow = nr, ncol = k)
    dist <- matrix(nrow = nr, ncol = k)
    for (i in 1:nr) {
      res <- ann$getNNsByVectorList(X[i, ], k, search_k, TRUE)
      if (length(res$item) != k) {
        stop(
          "search_k/n_trees settings were unable to find ", k,
          " neighbors for item ", i
        )
      }
      idx[i, ] <- res$item
      dist[i, ] <- res$distance
      search_progress$increment()
    }
  }

  list(idx = idx + 1, dist = dist)
}


FNN_nn <- function(X, k = 10, include_self = TRUE) {
  if (include_self) {
    k <- k - 1
  }

  fnn <- FNN::get.knn(X, k)
  idx <- fnn$nn.index
  dist <- fnn$nn.dist

  if (include_self) {
    idx <- cbind(seq_len(nrow(X)), idx)
    dist <- cbind(rep(0, nrow(X)), dist)
  }

  list(idx = idx, dist = dist)
}

dist_nn <- function(X, k, include_self = TRUE) {
  X <- as.matrix(X)

  if (!include_self) {
    k <- k + 1
  }

  nn_idx <- t(apply(X, 2, order))[, 1:k]
  nn_dist <- matrix(0, nrow = nrow(X), ncol = k)
  for (i in seq_len(nrow(nn_idx))) {
    nn_dist[i, ] <- X[i, nn_idx[i, ]]
  }

  if (!include_self) {
    nn_idx <- nn_idx[, 2:ncol(nn_idx)]
    nn_dist <- nn_dist[, 2:ncol(nn_dist)]
  }

  attr(nn_idx, "dimnames") <- NULL
  attr(nn_dist, "dimnames") <- NULL

  list(idx = nn_idx, dist = nn_dist)
}

sparse_nn <- function(X, k, include_self = TRUE) {
  if (include_self) {
    k <- k - 1
  }

  n <- nrow(X)
  nn_idx <- matrix(0, nrow = n, ncol = k)
  nn_dist <- matrix(0, nrow = n, ncol = k)

  for (i in 1:n) {
    dists <- X[, i]
    is_nonzero <- dists != 0
    dist_nonzero <- dists[is_nonzero]
    if (length(dist_nonzero) < k) {
      stop(
        "Row ", i, " of distance matrix has only ", length(dist_nonzero),
        " defined distances"
      )
    }

    k_order <- order(dist_nonzero)[1:k]

    idx_nonzero <- which(is_nonzero, arr.ind = TRUE)

    nn_idx[i, ] <- idx_nonzero[k_order]
    nn_dist[i, ] <- dist_nonzero[k_order]
  }

  if (include_self) {
    nn_idx <- cbind(1:n, nn_idx)
    nn_dist <- cbind(rep(0, n), nn_dist)
  }

  list(idx = nn_idx, dist = nn_dist)
}
