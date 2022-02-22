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
#'   dataset and hardware-dependent. If \code{NULL}, the transform will use the
#'   value provided in the \code{model}, if available. Default: \code{FALSE}. 
#' @param learning_rate Initial learning rate used in optimization of the
#'   coordinates. This overrides the value associated with the \code{model}.
#'   This should be left unspecified under most circumstances.
#' @param opt_args A list of optimizer parameters, used when 
#'   \code{batch = TRUE}. The default optimization method used is Adam (Kingma
#'   and Ba, 2014).
#'   \itemize{
#'     \item \code{method} The optimization method to use. Either \code{"adam"} 
#'     or \code{"sgd"} (stochastic gradient descent). Default: \code{"adam"}.
#'     \item \code{beta1} (Adam only). The weighting parameter for the
#'     exponential moving average of the first moment estimator. Effectively the
#'     momentum parameter. Should be a floating point value between 0 and 1.
#'     Higher values can smooth oscillatory updates in poorly-conditioned
#'     situations and may allow for a larger \code{learning_rate} to be
#'     specified, but too high can cause divergence. Default: \code{0.5}.
#'     \item \code{beta2} (Adam only). The weighting parameter for the
#'     exponential moving average of the uncentered second moment estimator.
#'     Should be a floating point value between 0 and 1. Controls the degree of
#'     adaptivity in the step-size. Higher values put more weight on previous
#'     time steps. Default: \code{0.9}.
#'     \item \code{eps} (Adam only). Intended to be a small value to prevent
#'     division by zero, but in practice can also affect convergence due to its
#'     interaction with \code{beta2}. Higher values reduce the effect of the
#'     step-size adaptivity and bring the behavior closer to stochastic gradient
#'     descent with momentum. Typical values are between 1e-8 and 1e-3. Default:
#'     \code{1e-7}.
#'     \item \code{alpha} The initial learning rate. Default: the value of the 
#'     \code{learning_rate} parameter.
#'   }
#'   If \code{NULL}, the transform will use the value provided in the 
#'   \code{model}, if available.
#' @param epoch_callback A function which will be invoked at the end of every
#'   epoch. Its signature should be:
#'   \code{(epoch, n_epochs, coords, fixed_coords)}, where:
#'   \itemize{
#'     \item \code{epoch} The current epoch number (between \code{1} and 
#'     \code{n_epochs}).
#'     \item \code{n_epochs} Number of epochs to use during the optimization of
#'     the embedded coordinates.
#'     \item \code{coords} The embedded coordinates as of the end of the current
#'     epoch, as a matrix with dimensions (N, \code{n_components}).
#'     \item \code{fixed_coords} The originally embedded coordinates from the
#'     \code{model}. These are fixed and do not change. A matrix with dimensions 
#'     (Nmodel, \code{n_components}) where \code{Nmodel} is the number of
#'     observations in the original data.
#'   }
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
                           batch = NULL,
                           learning_rate = NULL,
                           opt_args = NULL,
                           epoch_callback = NULL
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
  n_train_vertices <- nrow(train_embedding)
  ndim <- ncol(train_embedding)
  row.names(train_embedding) <- NULL
  # uwot model format should be changed so train embedding is stored transposed
  train_embedding <- t(train_embedding)

  method <- model$method
  scale_info <- model$scale_info
  metric <- model$metric
  nblocks <- length(metric)
  pca_models <- model$pca_models
  
  if (method == "leopold") {
    dens_scale <- model$dens_scale
    ai <- model$ai
    rad_coeff <- model$rad_coeff
  }
  
  if (is.null(batch)) {
    if (!is.null(model$batch)) {
      batch <- model$batch
    }
    else {
      batch <- FALSE
    }
  }
  
  if (is.null(opt_args)) {
    if (!is.null(model$opt_args)) {
      opt_args <- model$opt_args
    }
    else {
      opt_args <- list()
    }
  }

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
      indim <- dim(init)
      xdim <- c(n_vertices, ndim)
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
  localr <- NULL
  need_sigma <- method == "leopold" && nblocks == 1
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
    
    nnt <- nn_graph_t(nn)
    sknn_res <- smooth_knn(nnt$dist,
      local_connectivity = adjusted_local_connectivity,
      n_threads = n_threads,
      grain_size = grain_size,
      verbose = verbose,
      ret_sigma = need_sigma
    )
    graph_block <- sknn_res$matrix
    graph_blockv <- as.vector(graph_block)
    graph_block <- nn_to_sparse(nnt$idx, graph_blockv,
      self_nbr = FALSE,
      max_nbr_id = n_train_vertices,
      by_row = FALSE
    )
    
    if (is.null(localr) && need_sigma) {
      # because of the adjusted local connectivity rho is too small compared
      # to that used to generate the "training" data but sigma is larger, so
      # let's just stick with sigma + rho even though it tends to be an 
      # underestimate 
      localr <- sknn_res$sigma + sknn_res$rho
    }
    
    
    if (is.logical(init_weighted)) {
      embedding_block <-
        init_new_embedding(
          train_embedding = train_embedding,
          nn_idx = as.vector(nnt$idx),
          n_test_vertices = ncol(nnt$idx),
          graph = graph_blockv,
          weighted = init_weighted,
          n_threads = n_threads,
          grain_size = grain_size,
          verbose = verbose
        )
      if (is.null(embedding)) {
        embedding <- embedding_block
      }
      else {
        embedding <- embedding + embedding_block
      }
    }
    
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
    embedding <- t(init)
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
    
    n_head_vertices <- ncol(embedding)
    n_tail_vertices <- n_train_vertices
    
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
    if (method == "leopold") {
      # Use the linear model 2 log ai = -m log(localr) + c
      aj <- exp(0.5 * ((-log(localr) * rad_coeff[2]) + rad_coeff[1]))
      # Prevent too-small aj
      aj[aj < 0.1] <- 0.1
      method <- "leopold2"
    }

    method_args <- switch(method, 
      umap = list(a = a, b = b, gamma = gamma, approx_pow = approx_pow),
      leopold2 = list(ai = ai, aj = aj, b = b, ndim = ndim),
      list()
    )

    full_opt_args <- get_opt_args(opt_args, alpha)
    
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
      opt_args = full_opt_args,
      negative_sample_rate = negative_sample_rate,
      pcg_rand = pcg_rand,
      batch = batch,
      n_threads = n_sgd_threads,
      grain_size = grain_size,
      move_other = FALSE,
      verbose = verbose,
      epoch_callback = epoch_callback
    )
  }
  embedding <- t(embedding)
  tsmessage("Finished")
  if (!is.null(Xnames)) {
    row.names(embedding) <- Xnames
  }
  embedding
}

init_new_embedding <-
  function(train_embedding,
           nn_idx,
           n_test_vertices,
           graph,
           weighted = TRUE,
           n_threads = NULL,
           grain_size = 1,
           verbose = FALSE) {
    if (is.null(n_threads)) {
      n_threads <- default_num_threads()
    }
    avtype <- ifelse(weighted, "weighted ", "")
    tsmessage(
      "Initializing by ",
      avtype,
      "average of neighbor coordinates",
      pluralize("thread", n_threads, " using")
    )
    nn_weights <- NULL
    if (weighted) {
      nn_weights <- graph
    }
    init_transform_parallel(
      train_embedding = train_embedding,
      nn_index = nn_idx,
      n_test_vertices = n_test_vertices,
      nn_weights = nn_weights,
      n_threads = n_threads,
      grain_size = grain_size
    )
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
