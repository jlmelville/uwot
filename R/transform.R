#' Add New Points to an Existing Embedding
#'
#' Carry out an embedding of new data using an existing embedding. Requires
#' using the result of calling \code{\link{umap}} or \code{\link{tumap}} with
#' \code{ret_model = TRUE}.
#'
#' Note that some settings are incompatible with the production of a UMAP model
#' via \code{\link{umap}}: external neighbor data (passed via a list to the
#' argument of the \code{nn_method} parameter), and factor columns that were
#' included in the UMAP calculation via the \code{metric} parameter. In the
#' latter case, the model produced is based only on the numeric data.
#' A transformation is possible, but factor columns in the new data are ignored.
#'
#' @param X The new data to be transformed, either a matrix of data frame. Must
#'   have the same columns in the same order as the input data used to generate
#'   the \code{model}.
#' @param model Data associated with an existing embedding.
#' @param init_weighted If \code{TRUE}, then initialize the embedded coordinates
#'   of \code{X} using a weighted average of the coordinates of the nearest
#'   neighbors from the original embedding in \code{model}, where the weights
#'   used are the edge weights from the UMAP smoothed knn distances. Otherwise,
#'   use an unweighted average.
#' @param search_k Number of nodes to search during the neighbor retrieval. The
#'   larger k, the more the accurate results, but the longer the search takes.
#'   Default is the value used in building the \code{model} is used.
#' @param tmpdir Temporary directory to store nearest neighbor indexes during
#'   nearest neighbor search. Default is \code{\link{tempdir}}. The index is
#'   only written to disk if \code{n_threads > 1}; otherwise, this parameter is
#'   ignored.
#' @param n_epochs Number of epochs to use during the optimization of the
#'   embedded coordinates. A value between \code{30 - 100} is a reasonable trade
#'   off between speed and thoroughness. By default, this value is set to one
#'   third the number of epochs used to build the \code{model}.
#' @param n_threads Number of threads to use, (except during stochastic gradient
#'   descent). Default is half the number of concurrent threads supported by the
#'   system.
#' @param n_sgd_threads Number of threads to use during stochastic gradient
#'   descent. If set to > 1, then results will not be reproducible, even if
#'   `set.seed` is called with a fixed seed before running.
#' @param grain_size Minimum batch size for multithreading. If the number of
#'   items to process in a thread falls below this number, then no threads will
#'   be used. Used in conjunction with \code{n_threads} and
#'   \code{n_sgd_threads}.
#' @param verbose If \code{TRUE}, log details to the console.
#' @return A matrix of coordinates for \code{X} transformed into the space
#'   of the \code{model}.
#' @examples
#'
#' iris_train <- iris[1:100, ]
#' iris_test <- iris[101:150, ]
#'
#' # You must set ret_model = TRUE to return extra data needed
#' iris_train_umap <- umap(iris_train, ret_model = TRUE)
#' iris_test_umap <- umap_transform(iris_test, iris_train_umap)
#' @export
umap_transform <- function(X, model,
                           init_weighted = TRUE,
                           search_k = NULL,
                           tmpdir = tempdir(),
                           n_epochs = NULL,
                           n_threads = NULL,
                           n_sgd_threads = 0,
                           grain_size = 1,
                           verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (!all_nn_indices_are_loaded(model)) {
    stop("cannot use model: NN index is unloaded." ,
         " Try reloading with `load_uwot`")
  }
  
  if (is.null(n_epochs)) {
    n_epochs <- model$n_epochs
  }
  if (is.null(search_k)) {
    search_k <- model$search_k
  }

  nn_index <- model$nn_index
  n_neighbors <- model$n_neighbors
  local_connectivity <- model$local_connectivity
  train_embedding <- model$embedding
  method <- model$method
  scale_info <- model$scale_info
  metric <- model$metric
  pca_models <- model$pca_models

  a <- model$a
  b <- model$b
  gamma <- model$gamma
  alpha <- model$alpha
  negative_sample_rate <- model$negative_sample_rate
  approx_pow <- model$approx_pow
  norig_col <- model$norig_col
  pcg_rand <- model$pcg_rand
  if (is.null(pcg_rand)) {
    tsmessage("Using PCG for random number generation")
    pcg_rand <- TRUE
  }

  if (ncol(X) != norig_col) {
    stop("Incorrect dimensions: X must have ", norig_col, " columns")
  }

  if (methods::is(X, "data.frame")) {
    indexes <- which(vapply(X, is.numeric, logical(1)))
    if (length(indexes) == 0) {
      stop("No numeric columns found")
    }
    X <- as.matrix(X[, indexes])
  }
  n_vertices <- nrow(X)
  tsmessage(
    "Read ", n_vertices, " rows and found ", ncol(X),
    " numeric columns"
  )

  if (!is.null(scale_info)) {
    X <- apply_scaling(X, scale_info = scale_info, verbose = verbose)
  }

  adjusted_local_connectivity <- max(0, local_connectivity - 1.0)

  nblocks <- length(metric)
  graph <- NULL
  embedding <- NULL
  for (i in 1:nblocks) {
    tsmessage("Processing block ", i, " of ", nblocks)
    if (nblocks == 1) {
      ann <- nn_index
      Xsub <- X
    }
    else {
      ann <- nn_index[[i]]
      subset <- metric[[i]]
      if (is.list(subset)) {
        subset <- lsplit_unnamed(subset)$unnamed[[1]]
      }
      Xsub <- X[, subset, drop = FALSE]
    }

    if (!is.null(pca_models) && !is.null(pca_models[[as.character(i)]])) {
      Xsub <- apply_pca(
        X = Xsub, pca_res = pca_models[[as.character(i)]],
        verbose = verbose
      )
    }
    nn <- annoy_search(Xsub,
      k = n_neighbors, ann = ann, search_k = search_k,
      tmpdir = tmpdir,
      n_threads = n_threads, grain_size = grain_size,
      verbose = verbose
    )
    graph_block <- smooth_knn(nn,
      local_connectivity = adjusted_local_connectivity,
      n_threads = n_threads,
      grain_size = grain_size,
      verbose = verbose
    )

    embedding_block <- init_new_embedding(train_embedding, nn, graph_block,
      weighted = init_weighted,
      n_threads = n_threads,
      grain_size = grain_size, verbose = verbose
    )
    if (is.null(embedding)) {
      embedding <- embedding_block
    }
    else {
      embedding <- embedding + embedding_block
    }

    graph_block <- nn_to_sparse(nn$idx, as.vector(graph_block),
      self_nbr = FALSE,
      max_nbr_id = nrow(train_embedding)
    )
    if (is.null(graph)) {
      graph <- graph_block
    }
    else {
      graph <- set_intersect(graph, graph_block, weight = 0.5, reset = TRUE)
    }
  }

  if (nblocks > 1) {
    embedding <- embedding / nblocks
  }

  if (is.null(n_epochs)) {
    if (ncol(graph) <= 10000) {
      n_epochs <- 100
    }
    else {
      n_epochs <- 30
    }
  }
  else {
    n_epochs <- max(2, round(n_epochs / 3))
  }

  if (n_epochs > 0) {
    graph@x[graph@x < max(graph@x) / n_epochs] <- 0
    graph <- Matrix::drop0(graph)
    epochs_per_sample <- make_epochs_per_sample(graph@x, n_epochs)

    positive_head <- graph@i
    positive_tail <- Matrix::which(graph != 0, arr.ind = TRUE)[, 2] - 1

    tsmessage(
      "Commencing optimization for ", n_epochs, " epochs, with ",
      length(positive_head), " positive edges",
      pluralize("thread", n_sgd_threads, " using")
    )

    embedding <- t(embedding)
    train_embedding <- t(train_embedding)
    if (tolower(method) == "umap") {
      embedding <- optimize_layout_umap(
        head_embedding = embedding,
        tail_embedding = train_embedding,
        positive_head = positive_head,
        positive_tail = positive_tail,
        n_epochs = n_epochs,
        n_vertices = n_vertices,
        epochs_per_sample = epochs_per_sample,
        a = a, b = b, gamma = gamma,
        initial_alpha = alpha, negative_sample_rate,
        approx_pow = approx_pow,
        pcg_rand = pcg_rand,
        n_threads = n_sgd_threads,
        grain_size = grain_size,
        move_other = FALSE,
        verbose = verbose
      )
    }
    else {
      embedding <- optimize_layout_tumap(
        head_embedding = embedding,
        tail_embedding = train_embedding,
        positive_head = positive_head,
        positive_tail = positive_tail,
        n_epochs = n_epochs,
        n_vertices, epochs_per_sample,
        initial_alpha = alpha,
        negative_sample_rate = negative_sample_rate,
        pcg_rand = pcg_rand,
        n_threads = n_sgd_threads,
        grain_size = grain_size,
        move_other = FALSE,
        verbose = verbose
      )
    }
  }
  tsmessage("Finished")
  t(embedding)
}

init_new_embedding <- function(train_embedding, nn, graph, weighted = TRUE,
                               n_threads = NULL,
                               grain_size = 1, verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (weighted) {
    tsmessage(
      "Initializing by weighted average of neighbor coordinates",
      pluralize("thread", n_threads, " using")
    )
    embedding <- init_transform_parallel(train_embedding, nn$idx, graph,
      n_threads = n_threads,
      grain_size = grain_size
    )
  }
  else {
    tsmessage(
      "Initializing by average of neighbor coordinates",
      pluralize("thread", n_threads, " using")
    )
    embedding <- init_transform_av_parallel(train_embedding, nn$idx,
                                            n_threads = n_threads,
                                            grain_size = grain_size
    )
  }

  embedding
}


# Pure R implementation of (weighted) average. Superceded by C++ implementations
init_transform <- function(train_embedding, nn_index, weights = NULL) {
  nr <- nrow(nn_index)
  nc <- ncol(train_embedding)

  embedding <- matrix(nrow = nr, ncol = nc)
  if (is.null(weights)) {
    for (i in 1:nr) {
      nbr_embedding <- train_embedding[nn_index[i, ], ]
      embedding[i, ] <- apply(nbr_embedding, 2, mean)
    }
  }
  else {
    for (i in 1:nr) {
      nbr_embedding <- train_embedding[nn_index[i, ], ]
      nbr_weights <- weights[nn_index[i, ], i]
      embedding[i, ] <- apply(
        nbr_embedding, 2,
        function(x) {
          stats::weighted.mean(x, nbr_weights)
        }
      )
    }
  }

  embedding
}

apply_scaling <- function(X, scale_info, verbose = FALSE) {
  if (!is.null(scale_info[["scaled:range:min"]])) {
    tsmessage("Applying training data range scaling")
    X <- X - scale_info[["scaled:range:min"]]
    X <- X / scale_info[["scaled:range:max"]]
  }
  else if (!is.null(scale_info[["scaled:maxabs"]])) {
    tsmessage("Applying training data max-abs scaling")
    X <- scale(X, center = scale_info[["scaled:center"]], scale = FALSE)
    X <- X / scale_info[["scaled:maxabs"]]
  }
  else if (!is.null(scale_info[["scaled:colrange:min"]])) {
    tsmessage("Applying training data column range scaling")
    X <- sweep(X, 2, scale_info[["scaled:colrange:min"]])
    X <- sweep(X, 2, scale_info[["scaled:colrange:max"]], `/`)
  }
  else {
    tsmessage("Applying training data column filtering/scaling")
    X <- X[, scale_info[["scaled:nzvcols"]]]
    X <- scale(X,
      center = scale_info[["scaled:center"]],
      scale = scale_info[["scaled:scale"]]
    )
  }

  X
}
# Apply a previously calculated set of PCA rotations
apply_pca <- function(X, pca_res, verbose = FALSE) {
  tsmessage("Applying PCA reducing to ", ncol(X), " dimensions")
  if (!is.null(pca_res$center)) {
    X <- sweep(X, 2, pca_res$center)
  }
  X %*% pca_res$rotation
}

all_nn_indices_are_loaded <- function(model) {
  if (is.null(model$nn_index)) {
    stop("Invalid model: has no 'nn_index'")
  }
  if (is.list(model$nn_index)) {
    for (i in 1:length(model$nn_index)) {
      if (model$nn_index$getNTrees() == 0) {
        return(FALSE)
      }
    }
  }
  else if (model$nn_index$getNTrees() == 0) {
    return(FALSE)
  }
  TRUE
}
