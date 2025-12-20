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
#' @param nn_method Optional pre-calculated nearest neighbor data. There are
#'   two supported formats. The first is a list consisting of two elements:
#'   \itemize{
#'     \item \code{"idx"}. A \code{n_vertices x n_neighbors} matrix where
#'     \code{n_vertices} is the number of observations in \code{X}. The contents
#'     of the matrix should be the integer indexes of the data used to generate
#'     the \code{model}, which are the \code{n_neighbors}-nearest neighbors of
#'     the data to be transformed.
#'     \item \code{"dist"}. A \code{n_vertices x n_neighbors} matrix
#'     containing the distances of the nearest neighbors.
#'   }
#'   The second supported format is a sparse distance matrix of type
#'   \code{dgCMatrix}, with dimensions \code{n_model_vertices x n_vertices}.
#'   where \code{n_model_vertices} is the number of observations in the original
#'   data that generated the model. Distances should be arranged by column, i.e.
#'   a non-zero entry in row \code{j} of the \code{i}th column indicates that
#'   the \code{j}th observation in the original data used to generate the
#'   \code{model} is a nearest neighbor of the \code{i}th observation in the new
#'   data, with the distance given by the value of that element. In this format,
#'   a different number of neighbors is allowed for each observation, i.e.
#'   each column can contain a different number of non-zero values.
#'   Multiple nearest neighbor data (e.g. from two different pre-calculated
#'   metrics) can be passed by passing a list containing the nearest neighbor
#'   data lists as items.
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
#' @param ret_extra A vector indicating what extra data to return. May contain
#'   any combination of the following strings:
#'   \itemize{
#'     \item \code{"fgraph"} the high dimensional fuzzy graph (i.e. the fuzzy
#'       simplicial set of the merged local views of the input data). The graph
#'       is returned as a sparse matrix of class \link[Matrix]{dgCMatrix-class}
#'       with dimensions \code{NX} x \code{Nmodel}, where \code{NX} is the number
#'       of items in the data to transform in \code{X}, and \code{NModel} is
#'       the number of items in the data used to build the UMAP \code{model}.
#'       A non-zero entry (i, j) gives the membership strength of the edge
#'       connecting the vertex representing the ith item in \code{X} to the
#'       jth item in the data used to build the \code{model}. Note that the
#'       graph is further sparsified by removing edges with sufficiently low
#'       membership strength that they would not be sampled by the probabilistic
#'       edge sampling employed for optimization and therefore the number of
#'       non-zero elements in the matrix is dependent on \code{n_epochs}. If you
#'       are only interested in the fuzzy input graph (e.g. for clustering),
#'       setting \code{n_epochs = 0} will avoid any further sparsifying.
#'    \item \code{"nn"} the nearest neighbor graph for \code{X} with respect to
#'      the observations in the \code{model}. The graph will be returned as a
#'      list of two items: \code{idx} a matrix of indices, with as many rows
#'      as there are items in \code{X} and as many columns as there are nearest
#'      neighbors to be computed (this value is determined by the \code{model}).
#'      The indices are those of the rows of the data used to build the
#'      \code{model}, so they're not necessarily of much use unless you have
#'      access to that data. The second item, \code{dist} is a matrix of the
#'      equivalent distances, with the same dimensions as \code{idx}.
#'   }
#' @param seed Integer seed to use to initialize the random number generator
#'   state. Combined with \code{n_sgd_threads = 1} or \code{batch = TRUE}, this
#'   should give consistent output across multiple runs on a given installation.
#'   Setting this value is equivalent to calling \code{\link[base]{set.seed}},
#'   but it may be more convenient in some situations than having to call a
#'   separate function. The default is to not set a seed, in which case this
#'   function uses the behavior specified by the supplied \code{model}: If the
#'   model specifies a seed, then the model seed will be used to seed then
#'   random number generator, and results will still be consistent (if
#'   \code{n_sgd_threads = 1}). If you want to force the seed to not be set,
#'   even if it is set in \code{model}, set \code{seed = FALSE}.
#' @return A matrix of coordinates for \code{X} transformed into the space
#'   of the \code{model}, or if \code{ret_extra} is specified, a list
#'   containing:
#'   \itemize{
#'     \item \code{embedding} the matrix of optimized coordinates.
#'     \item if \code{ret_extra} contains \code{"fgraph"}, an item of the same
#'     name containing the high-dimensional fuzzy graph as a sparse matrix, of
#'     type \link[Matrix]{dgCMatrix-class}.
#'     \item if \code{ret_extra} contains \code{"sigma"}, returns a vector of
#'     the smooth knn distance normalization terms for each observation as
#'     \code{"sigma"} and a vector \code{"rho"} containing the largest
#'     distance to the locally connected neighbors of each observation.
#'     \item if \code{ret_extra} contains \code{"localr"}, an item of the same
#'     name containing a vector of the estimated local radii, the sum of
#'     \code{"sigma"} and \code{"rho"}.
#'     \item if \code{ret_extra} contains \code{"nn"}, an item of the same name
#'     containing the nearest neighbors of each item in \code{X} (with respect
#'     to the items that created the \code{model}).
#'   }
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
                           epoch_callback = NULL,
                           ret_extra = NULL,
                           seed = NULL) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (is.character(n_sgd_threads) && n_sgd_threads == "auto") {
    n_sgd_threads <- n_threads
  }
  if (!is.numeric(n_sgd_threads)) {
    stop("Unknown n_sgd_threads value: ", n_sgd_threads, " should be a positive
         integer or 'auto'")
  }
  nn_method <- normalize_nn_method(nn_method)
  if (is.null(nn_method)) {
    if (is.null(X)) {
      stop('argument "X" is missing, with no default')
    }
    if (is.null(model)) {
      stop('argument "model" is missing, with no default')
    }
    if (!all_nn_indices_are_loaded(model)) {
      stop(
        "cannot use model: NN index is unloaded.",
        " Try reloading with `load_uwot`"
      )
    }
  } else {
    if (!is.null(X)) {
      tsmessage('argument "nn_method" is provided, ignoring argument "X"')
      X <- NULL
    }
  }

  if (is.character(model$nn_method) &&
      model$nn_method == "hnsw" && !is_installed("RcppHNSW")) {
    stop(
      "This model requires the RcppHNSW package to be installed."
    )
  }

  if (is.character(model$nn_method) &&
      model$nn_method == "nndescent" && !is_installed("rnndescent")) {
    stop(
      "This model requires the rnndescent package to be installed."
    )
  }

  if (is.null(n_epochs)) {
    n_epochs <- model$n_epochs
    if (is.null(n_epochs)) {
      if (ncol(graph) <= 10000) {
        n_epochs <- 100
      } else {
        n_epochs <- 30
      }
    } else {
      n_epochs <- max(2, round(n_epochs / 3))
    }
  }

  # Handle setting the random number seed internally:
  # 1. If the user specifies seed = FALSE, definitely don't set the seed, even
  # if the model has a seed.
  # 2. If the user specifies seed = integer, then use that seed, even if the
  # model has a seed.
  # 3. If the user does not specify a seed, then use the model seed, if it
  # exists. Otherwise don't set a seed. Also use this code path if the user
  # sets seed = TRUE
  if (is.logical(seed) && !seed) {
    # do nothing
  }
  # handle the seed = TRUE case in this clause too
  else if (is.logical(seed) || is.null(seed)) {
    if (!is.null(model$seed)) {
      tsmessage("Setting model random seed ", model$seed)
      set.seed(model$seed)
    }
    # otherwise no model seed, so do nothing
  } else {
    tsmessage("Setting random seed ", seed)
    set.seed(seed)
  }

  if (is.null(search_k)) {
    search_k <- model$search_k
  }

  nn_index <- model$nn_index
  n_neighbors <- model$n_neighbors
  local_connectivity <- model$local_connectivity

  train_embedding <- model$embedding
  if (!is.matrix(train_embedding)) {
    # this should only happen if the user set
    # `n_epochs = 0, init = NULL, ret_model = TRUE`
    stop(
      "Invalid embedding coordinates: should be a matrix, but got ",
      paste0(class(train_embedding), collapse = " ")
    )
  }
  if (any(is.na(train_embedding))) {
    stop("Model embedding coordinates contains NA values")
  }
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
    aj <- model$ai
    rad_coeff <- model$rad_coeff
  }

  if (is.null(batch)) {
    if (!is.null(model$batch)) {
      batch <- model$batch
    } else {
      batch <- FALSE
    }
  }

  if (is.null(opt_args)) {
    if (!is.null(model$opt_args)) {
      opt_args <- model$opt_args
    } else {
      opt_args <- list()
    }
  }

  a <- model$a
  b <- model$b
  gamma <- model$gamma
  if (is.null(learning_rate)) {
    alpha <- model$alpha
  } else {
    alpha <- learning_rate
  }
  if (!is.numeric(alpha) || length(alpha) > 1 || alpha < 0) {
    stop("learning rate should be a positive number, not ", alpha)
  }
  negative_sample_rate <- model$negative_sample_rate
  approx_pow <- model$approx_pow
  norig_col <- model$norig_col

  rng_type <- model$rng_type
  if (is.null(rng_type)) {
    pcg_rand <- model$pcg_rand
    if (is.null(pcg_rand)) {
      rng_type <- "pcg"
    }
    else {
      if (pcg_rand) {
        rng_type <- "pcg"
      } else {
        rng_type <- "tausworthe"
      }
    }
  }
  num_precomputed_nns <- model$num_precomputed_nns
  binary_edge_weights <- model$binary_edge_weights
  if (is.null(binary_edge_weights)) {
    binary_edge_weights <- FALSE
  }

  # the number of model vertices
  n_vertices <- NULL
  Xnames <- NULL
  if (!is.null(X)) {
    if (!(methods::is(X, "data.frame") ||
          methods::is(X, "matrix") || is_sparse_matrix(X))) {
      stop("Unknown input data format")
    }
    if (!is.null(norig_col) && ncol(X) != norig_col) {
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
    if (n_vertices < 1) {
      stop("Not enough rows in X")
    }
    if (!is.null(row.names(X))) {
      Xnames <- row.names(X)
    }
    checkna(X)
  } else if (nn_is_precomputed(nn_method)) {
    # https://github.com/jlmelville/uwot/issues/97
    # In the case where the training model didn't use pre-computed neighbors
    # we treat it like it had one block
    if (num_precomputed_nns == 0) {
      num_precomputed_nns <- 1
    }
    # store single nn graph as a one-item list
    if (num_precomputed_nns == 1 && nn_is_single(nn_method)) {
      nn_method <- list(nn_method)
    }
    if (length(nn_method) != num_precomputed_nns) {
      stop(
        "Wrong # pre-computed neighbor data blocks, expected: ",
        num_precomputed_nns, " but got: ", length(nn_method)
      )
    }
    if (length(n_neighbors) != num_precomputed_nns) {
      stop(
        "Wrong # n_neighbor values (one per neighbor block), expected: ",
        num_precomputed_nns, " but got: ", length(n_neighbors)
      )
    }
    for (i in 1:num_precomputed_nns) {
      graph <- nn_method[[i]]

      if (is.list(graph)) {
        check_graph(graph,
          expected_rows = n_vertices,
          expected_cols = n_neighbors[[i]], bipartite = TRUE
        )
        if (is.null(n_vertices)) {
          n_vertices <- nrow(graph$idx)
        }
        if (is.null(Xnames)) {
          Xnames <- nn_graph_row_names(graph)
        }
      } else if (is_sparse_matrix(graph)) {
        # nn graph should have dims n_train_obs x n_test_obs
        graph <- Matrix::drop0(graph)
        if (is.null(n_vertices)) {
          n_vertices <- ncol(graph)
        }
        if (is.null(Xnames)) {
          Xnames <- colnames(graph)
        }
      } else {
        stop("Error: unknown neighbor graph format")
      }
    }
    nblocks <- num_precomputed_nns
  }

  if (!is.null(init)) {
    if (is.logical(init)) {
      init_weighted <- init
    } else if (is.character(init)) {
      init <- tolower(init)
      if (init == "average") {
        init_weighted <- FALSE
      } else if (init == "weighted") {
        init_weighted <- TRUE
      } else {
        stop("Unknown option for init: '", init, "'")
      }
    } else if (is.matrix(init)) {
      indim <- dim(init)
      xdim <- c(n_vertices, ndim)
      if (!all(indim == xdim)) {
        stop(
          "Initial embedding matrix has wrong dimensions, expected (",
          xdim[1], ", ", xdim[2], "), but was (",
          indim[1], ", ", indim[2], ")"
        )
      }
      if (any(is.na(init))) {
        stop("Initial embedding matrix coordinates contains NA values")
      }
      if (is.null(Xnames) && !is.null(row.names(init))) {
        Xnames <- row.names(init)
      }
      init_weighted <- NULL
    } else {
      stop("Invalid input format for 'init'")
    }
  }
  if (is.null(n_vertices)) {
    stop("Failed to read input correctly: invalid input format")
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
  sigma <- NULL
  rho <- NULL
  export_nns <- NULL
  ret_nn <- FALSE
  if ("nn" %in% ret_extra) {
    ret_nn <- TRUE
    export_nns <- list()
  }
  need_sigma <- (method == "leopold" && nblocks == 1) || "sigma" %in% ret_extra
  for (i in 1:nblocks) {
    tsmessage("Processing block ", i, " of ", nblocks)
    if (!is.null(X)) {
      if (nblocks == 1) {
        Xsub <- X
        ann <- nn_index
      } else {
        subset <- metric[[i]]
        if (is.list(subset)) {
          subset <- lsplit_unnamed(subset)$unnamed[[1]]
        }
        Xsub <- X[, subset, drop = FALSE]
        ann <- nn_index[[i]]
      }

      if (!is.null(pca_models) && !is.null(pca_models[[as.character(i)]])) {
        Xsub <- apply_pca(
          X = Xsub, pca_res = pca_models[[as.character(i)]],
          verbose = verbose
        )
      }
      if (is.null(ann$type) || startsWith(ann$type, "annoy")) {
        nn <- annoy_search(
          Xsub,
          k = n_neighbors,
          ann = ann,
          search_k = search_k,
          prep_data = TRUE,
          tmpdir = tmpdir,
          n_threads = n_threads,
          grain_size = grain_size,
          verbose = verbose
        )
      }
      else if (startsWith(ann$type, "hnsw")) {
        if (is.list(model$nn_args)) {
          nn_args <- model$nn_args
          nn_args_names <- names(nn_args)

          if ("ef" %in% nn_args_names) {
            ef <- nn_args$ef
          }
          else {
            ef <- 10
          }
        }

        nn <- hnsw_search(
          X,
          k = n_neighbors,
          ann = ann,
          ef = ef,
          n_threads = n_threads,
          verbose = verbose
        )

        # We use the L2 HNSW index for Euclidean so we need to process the
        # distances
        if (names(model$metric)[[1]] == "euclidean") {
          nn$dist <- sqrt(nn$dist)
        }
      }
      else if (startsWith(ann$type, "nndescent")) {
        nn <-
          nndescent_search(
            X,
            k = n_neighbors,
            ann = ann,
            nn_args = model$nn_args,
            n_threads = n_threads,
            verbose = verbose
          )
      }
      else {
        stop("Unknown nn method: ", ann$type)
      }
      if (ret_nn) {
        export_nns[[i]] <- nn
        names(export_nns)[[i]] <- ann$metric
      }
    } else if (is.list(nn_method)) {
      # otherwise we expect a list of NN graphs
      nn <- nn_method[[i]]
      if (ret_nn) {
        export_nns[[i]] <- nn
        names(export_nns)[[i]] <- "precomputed"
      }
    } else {
      stop(
        "Can't transform new data if X is NULL ",
        "and no sparse distance matrix available"
      )
    }

    osparse <- NULL
    if (is_sparse_matrix(nn)) {
      nn <- Matrix::drop0(nn)
      osparse <- order_sparse(nn)
      nn_idxv <- osparse$i + 1
      nn_distv <- osparse$x
      nn_ptr <- osparse$p
      n_nbrs <- diff(nn_ptr)
      if (any(n_nbrs < 1)) {
        stop("All observations need at least one neighbor")
      }
      target <- log2(n_nbrs)
      skip_first <- TRUE
    } else {
      nnt <- nn_graph_t(nn)
      if (length(n_neighbors) == nblocks) {
        # if model came from multiple different external neighbor data
        n_nbrs <- n_neighbors[[i]]
      } else {
        # multiple internal blocks
        n_nbrs <- n_neighbors
      }
      if (is.na(n_nbrs) || n_nbrs != nrow(nnt$idx)) {
        # original neighbor data was sparse, but we are using dense knn format
        # or n_neighbors doesn't match
        n_nbrs <- nrow(nnt$idx)
        tsmessage(
          "Possible mismatch with original vs new neighbor data ",
          "format, using ", n_nbrs, " nearest neighbors"
        )
      }
      target <- log2(n_nbrs)
      nn_ptr <- n_nbrs
      nn_distv <- as.vector(nnt$dist)
      nn_idxv <- as.vector(nnt$idx)
      skip_first <- TRUE
    }

    sknn_res <- smooth_knn(
      nn_dist = nn_distv,
      nn_ptr = nn_ptr,
      skip_first = skip_first,
      target = target,
      local_connectivity = adjusted_local_connectivity,
      n_threads = n_threads,
      grain_size = grain_size,
      verbose = verbose,
      ret_sigma = TRUE
    )
    if (is.null(localr) && need_sigma) {
      # because of the adjusted local connectivity rho is too small compared
      # to that used to generate the "training" data but sigma is larger, so
      # let's just stick with sigma + rho even though it tends to be an
      # underestimate
      sigma <- sknn_res$sigma
      rho <- sknn_res$rho
      localr <- sknn_res$sigma + sknn_res$rho
    }

    graph_blockv <- sknn_res$matrix
    if (is_sparse_matrix(nn)) {
      graph_block <- Matrix::sparseMatrix(
        j = osparse$i, p = osparse$p, x = graph_blockv,
        dims = rev(osparse$dims), index1 = FALSE
      )
    } else {
      graph_block <- nn_to_sparse(nn_idxv, n_vertices, graph_blockv,
        self_nbr = FALSE,
        max_nbr_id = n_train_vertices,
        by_row = FALSE
      )
    }

    if (is.logical(init_weighted)) {
      embedding_block <-
        init_new_embedding(
          train_embedding = train_embedding,
          nn_idx = nn_idxv,
          n_test_vertices = n_vertices,
          graph = graph_blockv,
          weighted = init_weighted,
          n_threads = n_threads,
          grain_size = grain_size,
          verbose = verbose
        )
      if (is.null(embedding)) {
        embedding <- embedding_block
      } else {
        embedding <- embedding + embedding_block
      }
    }

    if (is.null(graph)) {
      graph <- graph_block
    } else {
      graph <- set_intersect(graph, graph_block,
        weight = 0.5,
        reset_connectivity = FALSE
      )
    }
  }

  if (is.logical(init_weighted)) {
    if (nblocks > 1) {
      embedding <- embedding / nblocks
    }
  } else {
    tsmessage("Initializing from user-supplied matrix")
    embedding <- t(init)
  }

  if (binary_edge_weights) {
    tsmessage("Using binary edge weights")
    graph@x <- rep(1, length(graph@x))
  }

  if (batch) {
    # This is the same arrangement as Python UMAP
    graph <- Matrix::t(graph)
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
      # ordered indices of the new data nodes. Coordinates are updated
      # during optimization
      positive_head <- Matrix::which(graph != 0, arr.ind = TRUE)[, 2] - 1
      # unordered indices of the model nodes (some may not have any incoming
      # edges), these coordinates will NOT update during the optimization
      positive_tail <- graph@i
    } else {
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
      ai <- exp(0.5 * ((-log(localr) * rad_coeff[2]) + rad_coeff[1]))
      # Prevent too-small/large aj
      min_ai <- min(sqrt(a * 10^(-2 * dens_scale)), 0.1)
      ai[ai < min_ai] <- min_ai
      max_ai <- sqrt(a * 10^(2 * dens_scale))
      ai[ai > max_ai] <- max_ai
      method <- "leopold2"
    }

    method_args <- switch(method,
      umap = list(a = a, b = b, gamma = gamma, approx_pow = approx_pow),
      tumap = list(gamma = gamma),
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
      rng_type = rng_type,
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
  if (length(ret_extra) > 0) {
    res <- list(embedding = embedding)
    for (name in ret_extra) {
      if (name == "fgraph") {
        if (batch) {
          # #129: we transposed graph in the batch=TRUE case (#118) but need to
          # transpose back for export
          graph <- Matrix::t(graph)
        }
        res$fgraph <- graph
      }
      if (name == "sigma") {
        res$sigma <- sigma
        res$rho <- rho
      }
      if (name == "localr" && !is.null(localr)) {
        res$localr <- localr
      }
      if (ret_nn && !is.null(export_nns)) {
        res$nn <- export_nns
      }
    }
  } else {
    res <- embedding
  }

  res
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
  } else {
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
  } else if (!is.null(scale_info[["scaled:maxabs"]])) {
    tsmessage("Applying training data max-abs scaling")
    X <- scale(X, center = scale_info[["scaled:center"]], scale = FALSE)
    X <- X / scale_info[["scaled:maxabs"]]
  } else if (!is.null(scale_info[["scaled:colrange:min"]])) {
    tsmessage("Applying training data column range scaling")
    X <- sweep(X, 2, scale_info[["scaled:colrange:min"]])
    X <- sweep(X, 2, scale_info[["scaled:colrange:max"]], `/`)
  } else {
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
