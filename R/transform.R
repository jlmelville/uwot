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
#' @param nn_method Optional pre-calculated nearest neighbor data. 
#' 
#' The format is 
#'   a list consisting of two elements:
#'   \itemize{
#'     \item \code{"idx"}. A \code{n_vertices x n_neighbors} matrix where 
#'     \code{n_vertices} is the number of items to be transformed. The contents 
#'     of the matrix should be the integer indexes of the data used to generate
#'     the \code{model}, which are the \code{n_neighbors}-nearest neighbors of
#'     the data to be transformed.
#'     \item \code{"dist"}. A \code{n_vertices x n_neighbors} matrix
#'     containing the distances of the nearest neighbors.
#'   }
#'   Multiple nearest neighbor data (e.g. from two different pre-calculated
#'   metrics) can be passed by passing a list containing the nearest neighbor
#'   data lists as items.
#'   The \code{X} parameter is ignored when using pre-calculated nearest 
#'   neighbor data.
#' @param init_weighted If \code{TRUE}, then initialize the embedded coordinates
#'   of \code{X} using a weighted average of the coordinates of the nearest
#'   neighbors from the original embedding in \code{model}, where the weights
#'   used are the edge weights from the UMAP smoothed knn distances. Otherwise,
#'   use an un-weighted average.
#'   This parameter will be deprecated and removed at version 1.0 of this
#'   package. Use the \code{init} parameter as a replacement, replacing
#'   \code{init_weighted = TRUE} with \code{init = "weighted"} and
#'   \code{init_weighted = FALSE} with \code{init = "average"}.
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
#'   descent. If set to > 1, then be aware that if \code{batch = FALSE}, results
#'   will \emph{not} be reproducible, even if \code{set.seed} is called with a
#'   fixed seed before running. Set to \code{"auto"} to use the same value as
#'   \code{n_threads}.
#' @param grain_size Minimum batch size for multithreading. If the number of
#'   items to process in a thread falls below this number, then no threads will
#'   be used. Used in conjunction with \code{n_threads} and
#'   \code{n_sgd_threads}.
#' @param verbose If \code{TRUE}, log details to the console.
#' @param init how to initialize the transformed coordinates. One of:
#'   \itemize{
#'     \item \code{"weighted"} (The default). Use a weighted average of the
#'     coordinates of the nearest neighbors from the original embedding in
#'     \code{model}, where the weights used are the edge weights from the UMAP
#'     smoothed knn distances. Equivalent to \code{init_weighted = TRUE}.
#'     \item \code{"average"}. Use the mean average of the coordinates of 
#'     the nearest neighbors from the original embedding in \code{model}.
#'     Equivalent to \code{init_weighted = FALSE}.
#'     \item A matrix of user-specified input coordinates, which must have
#'     dimensions the same as \code{(nrow(X), ncol(model$embedding))}.
#'   }
#' This parameter should be used in preference to \code{init_weighted}.
#' @param batch If \code{TRUE}, then embedding coordinates are updated at the
#'   end of each epoch rather than during the epoch. In batch mode, results are
#'   reproducible with a fixed random seed even with \code{n_sgd_threads > 1},
#'   at the cost of a slightly higher memory use. You may also have to modify
#'   \code{learning_rate} and increase \code{n_epochs}, so whether this provides
#'   a speed increase over the single-threaded optimization is likely to be
#'   dataset and hardware-dependent.
#' @param learning_rate Initial learning rate used in optimization of the
#'   coordinates. This overrides the value associated with the \code{model}.
#'   This should be left unspecified under most circumstances.
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
umap_transform <- function(X = NULL, model = NULL,
                           nn_method = NULL,
                           init_weighted = TRUE,
                           search_k = NULL,
                           tmpdir = tempdir(),
                           n_epochs = NULL,
                           n_threads = NULL,
                           n_sgd_threads = 0,
                           grain_size = 1,
                           verbose = FALSE,
                           init = "weighted",
                           batch = FALSE,
                           learning_rate = NULL,
                           opt_args = NULL
                           ) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (is.null(nn_method)) {
    if (is.null(X)) {
      stop('argument "X" is missing, with no default')
    }
    if (is.null(model)) {
      stop('argument "model" is missing, with no default')
    }
    if (!all_nn_indices_are_loaded(model)) {
      stop("cannot use model: NN index is unloaded." ,
           " Try reloading with `load_uwot`")
    }
  }
  else {
    if (!is.null(X)) {
      tsmessage('argument "nn_method" is provided, ignoring argument "X"')
    }
  }
  
  if (is.null(n_epochs)) {
    n_epochs <- model$n_epochs
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
  nblocks <- length(metric)
  pca_models <- model$pca_models

  a <- model$a
  b <- model$b
  gamma <- model$gamma
  if (is.null(learning_rate)) {
    alpha <- model$alpha
  }
  else {
    alpha <- learning_rate
  }
  if (! is.numeric(alpha) || length(alpha) > 1 || alpha < 0) {
    stop("learning rate should be a positive number, not ", alpha)
  }
  negative_sample_rate <- model$negative_sample_rate
  approx_pow <- model$approx_pow
  norig_col <- model$norig_col
  pcg_rand <- model$pcg_rand
  if (is.null(pcg_rand)) {
    tsmessage("Using PCG for random number generation")
    pcg_rand <- TRUE
  }

  # the number of model vertices
  n_vertices <- NULL
  Xnames <- NULL
  if (!is.null(X)){
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
    if (!is.null(row.names(X))) {
      Xnames <- row.names(X)
    }
    checkna(X)
  } else if (nn_is_precomputed(nn_method)) {
    if (nblocks == 1) {
      if (length(nn_method) == 1) {
        graph <- nn_method[[1]]
      }
      else {
        graph <- nn_method
      }
      n_vertices <- nrow(graph$idx)
      check_graph(graph, n_vertices, n_neighbors)
      if (is.null(Xnames)) {
        Xnames <- nn_graph_row_names(graph)
      }
    }
    else {
      stopifnot(length(nn_method) == nblocks)
      for (i in 1:nblocks) {
        graph <- nn_method[[i]]
        if (is.null(n_vertices)) {
          n_vertices <- nrow(graph$idx)
        }
        check_graph(graph, n_vertices, n_neighbors)
        if (is.null(Xnames)) {
          Xnames <- nn_graph_row_names(graph)
        }
      }
    }
  }
  
  if (!is.null(init)) {
    if (is.logical(init)) {
      init_weighted <- init
    }
    else if (is.character(init)) {
      init <- tolower(init)
      if (init == "average") {
        init_weighted <- FALSE
      }
      else if (init == "weighted") {
        init_weighted <- TRUE
      }
      else {
        stop("Unknown option for init: '", init, "'")
      }
    }
    else if (is.matrix(init)) {
      indim <- c(nrow(init), ncol(init))
      xdim <- c(n_vertices, ncol(train_embedding))
      if (!all(indim == xdim)) {
        stop("Initial embedding matrix has wrong dimensions, expected (",
             xdim[1], ", ", xdim[2], "), but was (",
             indim[1], ", ", indim[2], ")")
      }
      if (is.null(Xnames) && !is.null(row.names(init))) {
        Xnames <- row.names(init)
      }
      init_weighted <- NULL
    }
    else {
      stop("Invalid input format for 'init'")
    }
  }
  if (verbose) {
    x_is_matrix <- methods::is(X, "matrix")
    tsmessage("Read ", n_vertices, " rows", appendLF = !x_is_matrix)
    if (x_is_matrix) {
      tsmessage(" and found ", ncol(X), " numeric columns", time_stamp = FALSE)
    }
  }

  if (!is.null(scale_info)) {
    X <- apply_scaling(X, scale_info = scale_info, verbose = verbose)
  }

  adjusted_local_connectivity <- max(0, local_connectivity - 1.0)

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
    if (!is.null(X)) {
      nn <- annoy_search(Xsub,
                         k = n_neighbors, ann = ann, search_k = search_k,
                         prep_data = TRUE,
                         tmpdir = tmpdir,
                         n_threads = n_threads, grain_size = grain_size,
                         verbose = verbose
      )
    } else if (nn_is_precomputed(nn_method)) {
      if (nblocks == 1 && !is.null(nn_method$idx)) {
        # When there's only one block, the NN graph can be passed directly
        nn <- nn_method
      }
      else {
        # otherwise we expect a list of NN graphs
        nn <- nn_method[[i]]
      }
    }
    graph_block <- smooth_knn(nn,
      local_connectivity = adjusted_local_connectivity,
      n_threads = n_threads,
      grain_size = grain_size,
      verbose = verbose
    )

    if (is.logical(init_weighted)) {
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

  if (is.logical(init_weighted)) {
    if (nblocks > 1) {
      embedding <- embedding / nblocks
    }
  }
  else {
    tsmessage("Initializing from user-supplied matrix")
    embedding <- init
  }

  if (n_epochs > 0) {
    graph@x[graph@x < max(graph@x) / n_epochs] <- 0
    graph <- Matrix::drop0(graph)

    # Edges are (i->j) where i (head) is from the new data and j (tail) is
    # in the model data
    # Unlike embedding of initial data, the edge list is therefore NOT symmetric
    # i.e. the presence of (i->j) does NOT mean (j->i) is also present because
    # i and j now come from different data
    if (batch) {
      # This is the same arrangement as Python UMAP
      graph <- Matrix::t(graph)
      # ordered indices of the new data nodes. Coordinates are updated 
      # during optimization
      positive_head <- Matrix::which(graph != 0, arr.ind = TRUE)[, 2] - 1
      # unordered indices of the model nodes (some may not have any incoming
      # edges), these coordinates will NOT update during the optimization
      positive_tail <- graph@i
    }
    else {
      # unordered indices of the new data nodes. Coordinates are updated 
      # during optimization
      positive_head <- graph@i
      # ordered indices of the model nodes (some may not have any incoming edges)
      # these coordinates will NOT update during the optimization
      positive_tail <- Matrix::which(graph != 0, arr.ind = TRUE)[, 2] - 1
    }
    
    n_head_vertices <- nrow(embedding)
    n_tail_vertices <- nrow(train_embedding)
    
    # if batch = TRUE points into the head (length == n_tail_vertices)
    # if batch = FALSE, points into the tail (length == n_head_vertices)
    positive_ptr <- graph@p

    epochs_per_sample <- make_epochs_per_sample(graph@x, n_epochs)
    
    tsmessage(
      "Commencing optimization for ", n_epochs, " epochs, with ",
      length(positive_head), " positive edges",
      pluralize("thread", n_sgd_threads, " using")
    )

    method <- tolower(method)
    if (method == "umap") {
      method_args <- list(a = a, b = b, gamma = gamma, approx_pow = approx_pow)
    }
    else {
      method_args <- list()
    }
    
    if (is.null(opt_args)) {
      opt_args <- list("sgd", alpha = alpha / 4, momentum = 0, alpha_update = "linear_decay")
    }
    embedding <- t(embedding)
    row.names(train_embedding) <- NULL
    train_embedding <- t(train_embedding)
    embedding <- optimize_layout_r(
      head_embedding = embedding,
      tail_embedding = train_embedding,
      positive_head = positive_head,
      positive_tail = positive_tail,
      positive_ptr = positive_ptr,
      n_epochs = n_epochs,
      n_head_vertices = n_head_vertices,
      n_tail_vertices = n_tail_vertices,
      epochs_per_sample = epochs_per_sample,
      method = tolower(method),
      method_args = method_args,
      initial_alpha = alpha / 4.0,
      opt_args = opt_args,
      negative_sample_rate = negative_sample_rate,
      pcg_rand = pcg_rand,
      batch = batch,
      n_threads = n_sgd_threads,
      grain_size = grain_size,
      move_other = FALSE,
      verbose = verbose
    )
    embedding <- t(embedding)
  }
  tsmessage("Finished")
  if (!is.null(Xnames)) {
    row.names(embedding) <- Xnames
  }
  embedding
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
