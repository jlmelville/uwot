find_nn <- function(X, k, include_self = TRUE, method = "fnn",
                    n_trees = 50, search_k = 2 * k * n_trees,
                    n_threads = RcppParallel::defaultNumThreads() / 2,
                    grain_size = 1000,
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
      res <- annoy_nn(X, k = k, include_self = include_self,
                      n_trees = n_trees, search_k = search_k,
                      n_threads = n_threads,
                      verbose = verbose)
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
                     n_trees = 50, search_k = 2 * k * n_trees,
                     n_threads = RcppParallel::defaultNumThreads() / 2,
                     grain_size = 1000,
                     verbose = FALSE) {
  nr <- nrow(X)
  nc <- ncol(X)
  ann <- methods::new(RcppAnnoy::AnnoyEuclidean, nc)

  tsmessage("Building Annoy index")
  progress <- Progress$new(max = nr, display = verbose)

  for (i in 1:nr) {
    ann$addItem(i - 1, X[i, ])
    progress$increment()
  }

  ann$build(n_trees)
  index_file = tempfile()
  ann$save(index_file)

  if (!include_self) {
    k <- k + 1
  }

  if (n_threads > 0) {
    RcppParallel::setThreadOptions(numThreads = n_threads)
    tsmessage("Searching Annoy index with ", n_threads, " thread", ifelse(n_threads > 1, "s", ""))
    res <- annoy_euclidean_nns(index_file,
                               X,
                               k, search_k,
                               grain_size = grain_size,
                               verbose = verbose)
    idx <- res$item
    dist <- res$distance

    unlink(index_file)
  }
  else {
    tsmessage("Searching Annoy index")
    search_progress <- Progress$new(max = nr, display = verbose)
    idx <- matrix(nrow = nr, ncol = k)
    dist <- matrix(nrow = nr, ncol = k)
    for (i in 1:nr) {
      res <- ann$getNNsByItemList(i - 1, k, search_k, TRUE)
      idx[i, ] <- res$item
      dist[i, ] <- res$distance
      search_progress$increment()
    }
  }

  if (!include_self) {
    idx <- idx[, -1]
    dist <- dist[, -1]
    k <- k - 1
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
    idx <- cbind(1:nrow(X), idx)
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
  for (i in 1:nrow(nn_idx)) {
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
      stop("Row ", i, " of distance matrix has only ", length(dist_nonzero), " defined distances")
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
