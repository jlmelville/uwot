find_nn <- function(
  X,
  k,
  include_self = TRUE,
  method = "fnn",
  metric = "euclidean",
  n_trees = 50,
  search_k = 2 * k * n_trees,
  nn_args = nn_args,
  tmpdir = tempdir(),
  n_threads = NULL,
  n_build_threads = NULL,
  grain_size = 1,
  ret_index = FALSE,
  sparse_is_distance = TRUE,
  verbose = FALSE
) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (is.null(n_build_threads)) {
    n_build_threads <- n_threads
  }

  if (inherits(X, "dist")) {
    res <- dist_nn(X, k, include_self = include_self, verbose = verbose)
  } else if (sparse_is_distance && is_sparse_matrix(X)) {
    # sparse distance matrix
    if (Matrix::isTriangular(X)) {
      res <- sparse_tri_nn(X, k, include_self = include_self, verbose = verbose)
    } else {
      res <- sparse_nn(X, k, include_self = include_self, verbose = verbose)
    }
  } else {
    if (is_sparse_matrix(X) && method != "nndescent") {
      stop("Sparse matrix input only supported for nndescent method.")
    }
    # normal matrix
    switch(
      method,
      "fnn" = {
        res <- FNN_nn(X, k = k, include_self = include_self)
      },
      "annoy" = {
        nn_args_names <- names(nn_args)

        if ("n_trees" %in% nn_args_names) {
          n_trees <- nn_args$n_trees
        }

        res <- annoy_nn(
          X,
          k = k,
          metric = metric,
          n_trees = n_trees,
          search_k = annoy_search_k(nn_args, search_k),
          tmpdir = tmpdir,
          n_threads = n_threads,
          ret_index = ret_index,
          verbose = verbose
        )
      },
      "hnsw" = {
        nn_args$X <- X
        nn_args$k <- k
        nn_args$metric <- metric
        nn_args$n_threads <- n_threads
        nn_args$n_build_threads <- n_build_threads
        nn_args$verbose <- verbose
        nn_args$ret_index <- ret_index

        res <- do.call(hnsw_nn, nn_args)
      },
      "nndescent" = {
        res <- nndescent_nn(
          X,
          k = k,
          metric = metric,
          nn_args = nn_args,
          n_threads = n_threads,
          n_build_threads = n_build_threads,
          ret_index = ret_index,
          verbose = verbose
        )
      },
      stop("Unknown method: ", method)
    )
    if (include_self && method != "fnn") {
      res <- self_first(res)
    }
  }
  res
}

annoy_search_k <- function(nn_args, search_k) {
  if (is.null(nn_args[["search_k"]])) search_k else nn_args[["search_k"]]
}

# allow (idx, dist) or (index, distance): convert the latter to the former
normalize_nn_graph <- function(graph) {
  if (is.null(graph[["idx"]]) && is.matrix(graph[["index"]])) {
    graph[["idx"]] <- graph[["index"]]
    graph[["index"]] <- NULL
  }
  if (is.null(graph[["dist"]]) && is.matrix(graph[["distance"]])) {
    graph[["dist"]] <- graph[["distance"]]
    graph[["distance"]] <- NULL
  }
  graph
}

# in the knn graph case we will allow (idx, dist) or (index, distance)
# for compatibility with BiocNeighbors
normalize_nn_method <- function(nn_method) {
  if (is.null(nn_method) || !is.list(nn_method) || !is.null(nn_method$type)) {
    return(nn_method)
  }

  # the case where we have a single list with (idx, dist) or (index, distance)
  if (
    !is.null(nn_method[["idx"]]) ||
      !is.null(nn_method[["dist"]]) ||
      !is.null(nn_method[["index"]]) ||
      !is.null(nn_method[["distance"]])
  ) {
    nn_method <- normalize_nn_graph(nn_method)
    return(nn_method)
  }

  # the case where we have a list of knn graphs
  for (i in seq_along(nn_method)) {
    graph <- nn_method[[i]]
    if (is.list(graph) && is.null(graph$type)) {
      nn_method[[i]] <- normalize_nn_graph(graph)
    }
  }
  nn_method
}

# an nn graph not in a list
nn_is_single <- function(nn) {
  (is.list(nn) && !is.null(nn$idx)) || is_sparse_matrix(nn)
}

# TRUE if nn is a sparse matrix or an untagged list. This covers passing in
# a single nn graph, sparse distance matrix or list thereof, but excludes a
# tagged annoy index or a string like "euclidean"
nn_is_precomputed <- function(nn) {
  (is.list(nn) && is.null(nn$type)) || is_sparse_matrix(nn)
}

# TRUE if we are using an annoy index
nn_is_annoy <- function(ann) {
  is.list(ann) && !is.null(ann$type) && startsWith(ann$type, "annoy")
}

nn_is_hnsw <- function(ann) {
  is.list(ann) && !is.null(ann$type) && startsWith(ann$type, "hnsw")
}

# n_trees - number of trees to build when constructing the index. The more trees
# specified, the larger the index, but the better the results. largeVis uses 10
# trees for datasets with N = 10,000 observations, 20 trees for datasets up to N
# = 1,000,000, 50 trees for N up to 5,000,000 and 100 trees otherwise
# search_k - the number of nodes to search during the neighbor retrieval. The
# larger k, the more accurate results, but the longer the search takes. Default
# is k * n_trees.
#' @importFrom methods new
annoy_nn <- function(
  X,
  k = 10,
  metric = "euclidean",
  n_trees = 50,
  search_k = 2 * k * n_trees,
  tmpdir = tempdir(),
  n_threads = NULL,
  grain_size = 1,
  ret_index = FALSE,
  verbose = FALSE
) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  ann <- annoy_build(X, metric = metric, n_trees = n_trees, verbose = verbose)
  if (!ret_index) {
    on.exit(release_annoy_backing(annoy_backing(ann$ann)))
  }

  res <- annoy_search(
    X,
    k = k,
    ann = ann,
    search_k = search_k,
    tmpdir = tmpdir,
    n_threads = n_threads,
    prep_data = TRUE,
    owns_index = TRUE,
    grain_size = grain_size,
    verbose = verbose
  )

  nn_acc <- sum(res$idx == 1:nrow(X)) / nrow(X)
  tsmessage("Annoy recall = ", formatC(nn_acc * 100.0), "%")

  res <- list(idx = res$idx, dist = res$dist, recall = nn_acc)
  if (ret_index) {
    res$index <- ann
  }
  res
}

annoy_create <- function(metric, ndim) {
  rcppannoy <- create_ann(metric, ndim)

  list(
    ann = rcppannoy,
    type = "annoyv1",
    metric = metric,
    ndim = ndim
  )
}

annoy_build <- function(
  X,
  metric = "euclidean",
  n_trees = 50,
  verbose = FALSE
) {
  nr <- nrow(X)
  nc <- ncol(X)

  annoy <- annoy_create(metric, nc)

  if (metric == "correlation") {
    tsmessage("Annoy build: subtracting row means for correlation")
    X <- sweep(X, 1, rowMeans(X))
  }

  tsmessage(
    "Building Annoy index with metric = ",
    metric,
    ", n_trees = ",
    n_trees
  )
  ann <- annoy$ann
  nstars <- 50
  if (verbose && nr > nstars) {
    progress_for(
      nr,
      nstars,
      function(chunk_start, chunk_end) {
        for (i in chunk_start:chunk_end) {
          ann$addItem(i - 1, X[i, , drop = FALSE])
        }
      }
    )
  } else {
    for (i in 1:nr) {
      ann$addItem(i - 1, X[i, ])
    }
  }

  # Build index
  ann$build(n_trees)

  annoy
}

# create RcppAnnoy class from metric name with ndim dimensions
# Correlation uses AnnoyAngular, input data needs to be centered first
create_ann <- function(name, ndim) {
  ann <- switch(
    name,
    cosine = methods::new(RcppAnnoy::AnnoyAngular, ndim),
    manhattan = methods::new(RcppAnnoy::AnnoyManhattan, ndim),
    euclidean = methods::new(RcppAnnoy::AnnoyEuclidean, ndim),
    hamming = methods::new(RcppAnnoy::AnnoyHamming, ndim),
    correlation = methods::new(RcppAnnoy::AnnoyAngular, ndim),
    stop("BUG: unknown Annoy metric '", name, "'")
  )
}

# fetch the underlying RcppAnnoy class from inside an index
get_rcppannoy <- function(nni) {
  if (startsWith(class(nni), "Rcpp_Annoy")) {
    rcppannoy <- nni
  } else if (nn_is_annoy(nni)) {
    rcppannoy <- nni$ann
  } else if (nn_is_hnsw(nni)) {
    rcppannoy <- nni$ann
  } else {
    stop(
      "BUG: Found an unknown ann implementation of class: '",
      class(nni),
      "'"
    )
  }
  rcppannoy
}

# Search a pre-built Annoy index for neighbors of X
annoy_search <- function(
  X,
  k,
  ann,
  search_k = 100 * k,
  prep_data = FALSE,
  tmpdir = tempdir(),
  n_threads = NULL,
  grain_size = 1,
  verbose = FALSE,
  owns_index = FALSE,
  index_file = NULL
) {
  # newer NN structures hide impl in a tagged list
  if (nn_is_annoy(ann)) {
    lann <- ann
    ann <- lann$ann
    if (prep_data && lann$metric == "correlation") {
      tsmessage("Annoy search: subtracting row means for correlation")
      X <- sweep(X, 1, rowMeans(X))
    }
  }

  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (n_threads > 0) {
    annoy_res <- annoy_search_parallel(
      X = X,
      k = k,
      ann = ann,
      search_k = search_k,
      tmpdir = tmpdir,
      n_threads = n_threads,
      grain_size = grain_size,
      verbose = verbose,
      owns_index = owns_index,
      index_file = index_file
    )
    res <- list(idx = annoy_res$item + 1, dist = annoy_res$distance)
  } else {
    res <- annoy_search_serial(
      X = X,
      k = k,
      ann = ann,
      search_k = search_k,
      verbose = verbose
    )
  }
  # Convert from angular distance to the UMAP/sklearn definition of cosine
  # distance
  # Current Annoy README defines cosine distance as sqrt(2 - 2 cos(u,v))
  # where cos(u, v) is the cosine of the angle between two unit-scaled vectors
  # u and v (i.e. the cosine similarity). That expression is known to be
  # equivalent to the euclidean distance between u and v.
  # We shall convert back to 1 - cos(u, v) which is the definition of cosine
  # distance used by UMAP.
  if (methods::is(ann, "Rcpp_AnnoyAngular")) {
    res$dist <- 0.5 * res$dist * res$dist
  }

  res
}

annoy_search_serial <- function(
  X,
  k,
  ann,
  search_k = 100 * k,
  verbose = FALSE
) {
  tsmessage("Searching Annoy index, search_k = ", search_k)
  nr <- nrow(X)
  idx <- matrix(nrow = nr, ncol = k)
  dist <- matrix(nrow = nr, ncol = k)
  nstars <- 50
  if (verbose && nr > nstars) {
    progress_for(
      nr,
      nstars,
      function(chunk_start, chunk_end) {
        for (i in chunk_start:chunk_end) {
          res <- ann$getNNsByVectorList(X[i, ], k, search_k, TRUE)
          if (length(res$item) != k) {
            stop(
              "search_k/n_trees settings were unable to find ",
              k,
              " neighbors for item ",
              i
            )
          }
          idx[i, ] <<- res$item
          dist[i, ] <<- res$distance
        }
      }
    )
  } else {
    for (i in 1:nr) {
      res <- ann$getNNsByVectorList(X[i, ], k, search_k, TRUE)
      if (length(res$item) != k) {
        stop(
          "search_k/n_trees settings were unable to find ",
          k,
          " neighbors for item ",
          i
        )
      }
      idx[i, ] <- res$item
      dist[i, ] <- res$distance
    }
  }
  list(idx = idx + 1, dist = dist)
}

# A saved Annoy index remains mapped on Windows. Keep its backing file with
# the live index, including when several R models share that index.
annoy_backing <- function(ann) {
  get0(".uwot_backing", envir = ann, inherits = FALSE)
}

remember_annoy_backing <- function(ann, path, cleanup_dir = NULL) {
  backing <- annoy_backing(ann)
  if (is.null(backing)) {
    backing <- new.env(parent = emptyenv())
    backing$unload <- ann$unload
    reg.finalizer(backing, release_annoy_backing, onexit = TRUE)
    assign(".uwot_backing", backing, envir = ann)
  }
  backing$path <- normalizePath(path, mustWork = TRUE)
  backing$cleanup_dir <- cleanup_dir
}

release_annoy_backing <- function(backing, unload = TRUE) {
  if (!is.null(backing$path)) {
    if (unload) {
      backing$unload()
    }
    unlink(backing$path)
    if (!is.null(backing$cleanup_dir)) {
      unlink(backing$cleanup_dir, recursive = TRUE)
      backing$cleanup_dir <- NULL
    }
    backing$path <- NULL
  }
}

# RcppAnnoy returns NULL even when fopen, fwrite or fclose fails. Its native
# diagnostic uses REprintf. Preserve the caller's message connection, including
# when this save runs inside another capture.output().
annoy_save <- function(ann, path) {
  previous <- getConnection(sink.number(type = "message"))
  on.exit(sink(previous, type = "message"))
  diagnostics <- utils::capture.output(
    status <- ann$save(path),
    type = "message"
  )
  errors <- startsWith(diagnostics, "Unable to ")
  if (
    identical(status, FALSE) ||
      any(errors) ||
      !utils::file_test("-f", path) ||
      file.size(path) == 0
  ) {
    stop(
      "Could not save nearest neighbor index: ",
      path,
      if (any(errors)) {
        paste0(" (", paste(diagnostics[errors], collapse = "; "), ")")
      }
    )
  }
  if (length(diagnostics)) writeLines(diagnostics, previous)
}

annoy_search_parallel <- function(
  X,
  k,
  ann,
  search_k = 100 * k,
  tmpdir = tempdir(),
  n_threads = NULL,
  grain_size = 1,
  verbose = FALSE,
  owns_index = FALSE,
  index_file = NULL
) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  # An index loaded from a saved model already has a searchable backing file.
  if (is.null(index_file) || !file.exists(index_file)) {
    index_file <- annoy_backing(ann)$path
  }
  if (is.null(index_file)) {
    index_file <- tempfile(tmpdir = tmpdir)
    saved <- FALSE
    finished <- FALSE
    on.exit({
      if (owns_index && !finished) {
        ann$unload()
      }
      removed <- unlink(index_file) == 0L
      if (!removed && saved && .Platform$OS.type == "windows") {
        remember_annoy_backing(ann, index_file)
      }
    })
    tsmessage("Writing NN index file to temp file ", index_file)
    annoy_save(ann, index_file)
    saved <- TRUE
  }
  fsize <- file.size(index_file)
  tsmessage(
    "Searching Annoy index using ",
    pluralize("thread", n_threads),
    ", search_k = ",
    search_k
  )

  ann_class <- class(ann)
  metric <- switch(
    ann_class,
    Rcpp_AnnoyAngular = "cosine",
    Rcpp_AnnoyManhattan = "manhattan",
    Rcpp_AnnoyEuclidean = "euclidean",
    Rcpp_AnnoyHamming = "hamming",
    stop("BUG: unknown Annoy class '", ann_class, "'")
  )

  res <- annoy_search_parallel_cpp(
    index_file,
    X,
    k,
    search_k,
    metric = metric,
    n_threads = n_threads,
    grain_size = grain_size
  )
  if (any(res$item == -1)) {
    msg <- paste0(
      "search_k/n_trees settings were unable to find ",
      k,
      " neighbors for all items."
    )
    if (fsize > 2147483647) {
      msg <- paste0(
        msg,
        " Index file may have been too large to process.",
        " Try repeating with n_threads = 0, reducing n_trees,",
        " or reducing to a smaller dimensionality, e.g. pca = 50"
      )
    }
    stop(msg)
  }

  finished <- TRUE
  res
}

# Training searches can return a duplicate observation before the query itself,
# or omit self at the cutoff. Keep the backend's ordering among other ties.
self_first <- function(nn) {
  ids <- seq_len(nrow(nn$idx))
  k <- ncol(nn$idx)
  misplaced <- which(nn$idx[, 1] != ids)
  if (!length(misplaced) && all(nn$dist[, 1] == 0)) {
    return(nn)
  }
  for (i in misplaced) {
    other <- which(nn$idx[i, ] != i)[seq_len(k - 1L)]
    nn$idx[i, ] <- c(i, nn$idx[i, other])
    nn$dist[i, ] <- c(0, nn$dist[i, other])
  }
  nn$dist[, 1] <- 0
  nn
}

FNN_nn <- function(X, k = 10, include_self = TRUE) {
  # An explicit query returns k distinct candidates even for tied duplicates.
  # get.knn's implicit removal of self can discard a different observation.
  k_search <- if (include_self) k else k + 1L
  fnn <- FNN::get.knnx(X, X, k = k_search)
  nn <- self_first(list(idx = fnn$nn.index, dist = fnn$nn.dist))
  if (!include_self) {
    nn <- lapply(nn, function(x) x[, -1, drop = FALSE])
  }
  nn
}

dist_nn <- function(X, k, include_self = TRUE, verbose = FALSE) {
  tsmessage("Finding nearest neighbors from distance matrix")
  X <- as.matrix(X)
  # Sort self before every other observation, including exact duplicates.
  diag(X) <- -Inf
  k_search <- if (include_self) k else k + 1L
  idx <- t(apply(X, 2, order))[, seq_len(k_search), drop = FALSE]
  diag(X) <- 0
  if (!include_self) {
    idx <- idx[, -1, drop = FALSE]
  }
  dist <- matrix(0, nrow(X), ncol(idx))
  for (i in seq_len(nrow(X))) {
    dist[i, ] <- X[i, idx[i, ]]
  }
  dimnames(idx) <- NULL
  list(idx = idx, dist = dist)
}

sparse_nn <- function(X, k, include_self = TRUE, verbose = FALSE) {
  tsmessage("Finding nearest neighbors from sparse matrix")
  X <- as_neighbor_graph(X)
  n <- nrow(X)
  k_other <- if (include_self) k - 1L else k
  idx <- matrix(0L, n, k_other)
  dist <- matrix(0, n, k_other)
  for (i in seq_len(n)) {
    start <- X@p[i]
    entries <- start + seq_len(X@p[i + 1L] - start)
    # Distances on the diagonal cannot supply one of the k nonself neighbors.
    entries <- entries[X@i[entries] != i - 1L]
    degree <- length(entries)
    if (degree < k_other) {
      stop(
        "Row ",
        i,
        " of distance matrix has only ",
        degree,
        " defined distances"
      )
    }
    selected <- entries[order(X@x[entries])[seq_len(k_other)]]
    idx[i, ] <- X@i[selected] + 1L
    dist[i, ] <- X@x[selected]
  }
  if (include_self) {
    idx <- cbind(seq_len(n), idx)
    dist <- cbind(rep(0, n), dist)
  }
  list(idx = idx, dist = dist)
}

sparse_tri_nn <- function(X, k, include_self = TRUE, verbose = FALSE) {
  uplo <- if (Matrix::isTriangular(X, upper = TRUE)) "U" else "L"
  sparse_nn(Matrix::forceSymmetric(X, uplo), k, include_self, verbose)
}

is_binary_metric <- function(metric) {
  metric %in%
    c(
      "dice",
      "hamming",
      "jaccard",
      "kulsinski",
      "matching",
      "rogerstanimoto",
      "russellrao",
      "sokalmichener",
      "sokalsneath",
      "yule"
    )
}
