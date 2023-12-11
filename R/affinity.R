# set_op_mix_ratio = between 0 and 1 mixes in fuzzy set intersection
# set to 0 for intersection only
#' @import Matrix
fuzzy_set_union <- function(X, set_op_mix_ratio = 1) {
  XX <- X * Matrix::t(X)
  if (set_op_mix_ratio == 0) {
    Matrix::drop0(XX)
  } else if (set_op_mix_ratio == 1) {
    Matrix::drop0(X + Matrix::t(X) - XX)
  } else {
    Matrix::drop0(
      set_op_mix_ratio * (X + Matrix::t(X) - XX) + (1 - set_op_mix_ratio) * XX
    )
  }
}

# Calculate the (asymmetric) affinity matrix based on the nearest neighborhoods
# default target for calibration is the sum of affinities = log2(n_nbrs)
# nn distances should be stored column-wise
smooth_knn <- function(nn_dist,
                       nn_ptr = NULL,
                       skip_first = TRUE,
                       target = NULL,
                       local_connectivity = 1.0,
                       n_threads = NULL,
                       grain_size = 1,
                       ret_sigma = FALSE,
                       verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  tsmessage(
    "Commencing smooth kNN distance calibration",
    pluralize("thread", n_threads, " using"),
    appendLF = FALSE
  )
  if (length(target) == 1) {
    tsmessage(" with target n_neighbors = ", formatC(2^target), time_stamp = FALSE)
  } else {
    tsmessage(time_stamp = FALSE)
  }
  affinity_matrix_res <- smooth_knn_distances_parallel(
    nn_dist = nn_dist,
    nn_ptr = nn_ptr,
    skip_first = skip_first,
    target = target,
    n_iter = 64,
    local_connectivity = local_connectivity,
    tol = 1e-5,
    min_k_dist_scale = 1e-3,
    n_threads = n_threads,
    grain_size = grain_size,
    ret_sigma = ret_sigma
  )
  if (verbose && affinity_matrix_res$n_failures > 0) {
    tsmessage(affinity_matrix_res$n_failures, " smooth knn distance failures")
  }
  affinity_matrix_res
}

smooth_knn_matrix <- function(nn,
                              target = NULL,
                              local_connectivity = 1.0,
                              bandwidth = 1.0,
                              ret_sigma = FALSE,
                              n_threads = NULL,
                              grain_size = 1,
                              verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }

  osparse <- NULL
  if (is_sparse_matrix(nn)) {
    nn <- Matrix::drop0(nn)
    osparse <- order_sparse(nn)
    nn_dist <- osparse$x
    nn_ptr <- osparse$p
    n_nbrs <- diff(nn_ptr)
    if (any(n_nbrs < 1)) {
      stop("All observations need at least one neighbor")
    }
    if (is.null(target)) {
      # add 1 to n_nbrs to account for implicit self neighbor
      target <- log2(n_nbrs + 1) * bandwidth
    }
    skip_first <- FALSE
  } else {
    nnt <- nn_graph_t(nn)
    n_nbrs <- nrow(nnt$dist)
    if (is.null(target)) {
      target <- log2(n_nbrs) * bandwidth
    }
    nn_ptr <- n_nbrs
    nn_dist <- as.vector(nnt$dist)
    skip_first <- TRUE
  }
  affinity_matrix_res <- smooth_knn(
    nn_dist = nn_dist,
    nn_ptr = nn_ptr,
    skip_first = skip_first,
    target = target,
    local_connectivity = local_connectivity,
    ret_sigma = ret_sigma,
    n_threads = n_threads,
    grain_size = grain_size,
    verbose = verbose
  )

  v <- affinity_matrix_res$matrix
  if (is_sparse_matrix(nn)) {
    # use j instead of i to transpose it
    v <- Matrix::sparseMatrix(
      j = osparse$i, p = osparse$p, x = v,
      dims = osparse$dims, index1 = FALSE
    )
    Matrix::diag(v) <- 0.0
    v <- Matrix::drop0(v)
  } else {
    v <- nng_to_sparse(nnt$idx, v, self_nbr = TRUE, by_row = FALSE)
  }
  affinity_matrix_res$matrix <- v
  affinity_matrix_res
}

# Given nearest neighbor data and a measure of distance compute
# the fuzzy simplicial set (here represented as a fuzzy graph in the form of a
# sparse matrix) associated to the data. This is done by locally approximating
# geodesic distance at each point, creating a fuzzy simplicial set for each such
# point, and then combining all the local fuzzy simplicial sets into a global
# one via a fuzzy union
fuzzy_simplicial_set <- function(nn,
                                 target = NULL,
                                 set_op_mix_ratio = 1.0,
                                 local_connectivity = 1.0, bandwidth = 1.0,
                                 ret_sigma = FALSE,
                                 n_threads = NULL,
                                 grain_size = 1,
                                 verbose = FALSE) {
  affinity_matrix_res <- smooth_knn_matrix(
    nn = nn,
    target = target,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
    ret_sigma = ret_sigma,
    n_threads = n_threads,
    grain_size = grain_size,
    verbose = verbose
  )
  res <- fuzzy_set_union(affinity_matrix_res$matrix, set_op_mix_ratio = set_op_mix_ratio)
  if (ret_sigma) {
    res <- list(matrix = res)
    res$sigma <- affinity_matrix_res$sigma
    res$rho <- affinity_matrix_res$rho
  }
  res
}

symmetrize <- function(P) {
  0.5 * (P + Matrix::t(P))
}

perplexity_similarities <- function(nn, perplexity = NULL, ret_sigma = FALSE,
                                    n_threads = NULL,
                                    grain_size = 1,
                                    kernel = "gauss",
                                    verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  if (is.null(perplexity) && kernel != "knn") {
    stop("Must provide perplexity")
  }

  sigma <- NULL
  if (kernel == "gauss") {
    tsmessage(
      "Commencing calibration for perplexity = ", formatC(perplexity),
      pluralize("thread", n_threads, " using")
    )

    nnt <- nn_graph_t(nn)
    n_vertices <- ncol(nnt$dist)
    affinity_matrix_res <- calc_row_probabilities_parallel(
      nn_dist = as.vector(nnt$dist),
      n_vertices = n_vertices,
      perplexity = perplexity,
      ret_sigma = ret_sigma,
      n_threads = n_threads,
      grain_size = grain_size
    )
    if (verbose && affinity_matrix_res$n_failures > 0) {
      tsmessage(affinity_matrix_res$n_failures, " perplexity failures")
    }

    dint <- NULL
    if (ret_sigma && !is.null(affinity_matrix_res$sigma)) {
      # An analytical version of the "soft" correlation dimension estimate of
      # intrinsic dimensionality from multi-scale SNE by Lee et al (2015).
      # http://jlmelville.github.io/sneer/dimensionality.html
      d <- nnt$dist
      p <- affinity_matrix_res$matrix
      logp <- log(p + .Machine$double.eps)
      s <- affinity_matrix_res$sigma
      h <- -colSums(p * logp)
      lph <- sweep(logp, 2, h, `+`)
      dhdb <- colSums(d * d * p * lph)
      dint <- -2 * dhdb / (s * s)
    }

    affinity_matrix <- nng_to_sparse(nnt$idx, as.vector(affinity_matrix_res$matrix),
      self_nbr = TRUE, by_row = FALSE
    )
    if (!is.null(affinity_matrix_res$sigma)) {
      sigma <- affinity_matrix_res$sigma
    }
  } else {
    # knn kernel
    tsmessage("Using knn graph for input weights with k = ", ncol(nn$idx))
    # Make each row sum to 1, ignoring the self-index
    # i.e. diagonal will be zero
    affinity_matrix <- nng_to_sparse(nn$idx, val = 1 / (ncol(nn$idx) - 1))
    Matrix::diag(affinity_matrix) <- 0
    affinity_matrix <- Matrix::drop0(affinity_matrix)
  }
  res <- list(matrix = symmetrize(affinity_matrix))
  if (ret_sigma && !is.null(sigma)) {
    res$sigma <- sigma
    if (!is.null(dint)) {
      res$dint <- dint
    }
  }
  res
}

# Convert the matrix of NN indices to a sparse asymmetric matrix where each
# edge has a weight of val (scalar or vector)
# return a sparse matrix with dimensions of nrow(nn_idx) x max_nbr_id
nn_to_sparse <- function(nn_idxv, n_obs, val = 1, self_nbr = FALSE,
                         max_nbr_id = NULL, by_row = TRUE) {
  n_nbrs <- length(nn_idxv) / n_obs


  if (is.null(max_nbr_id)) {
    max_nbr_id <- ifelse(self_nbr, n_obs, max(nn_idxv))
  }

  if (length(val) == 1) {
    xs <- rep(val, n_obs * n_nbrs)
  } else {
    xs <- val
  }
  if (by_row) {
    is <- rep(1:n_obs, times = n_nbrs)
  } else {
    is <- rep(1:n_obs, each = n_nbrs)
  }

  dims <- c(n_obs, max_nbr_id)
  res <- Matrix::sparseMatrix(i = is, j = nn_idxv, x = xs, dims = dims)

  if (self_nbr) {
    Matrix::diag(res) <- 0
    res <- Matrix::drop0(res)
  }
  res
}

nng_to_sparse <- function(nn_idx, val = 1, self_nbr = FALSE,
                          max_nbr_id = NULL, by_row = TRUE) {
  if (by_row) {
    n_obs <- nrow(nn_idx)
  } else {
    n_obs <- ncol(nn_idx)
  }

  nn_to_sparse(as.vector(nn_idx), n_obs,
    val = val, self_nbr = self_nbr,
    max_nbr_id = max_nbr_id, by_row = by_row
  )
}

# transpose the index and distance matrix
nn_graph_t <- function(nn_graph) {
  list(idx = t(nn_graph$idx), dist = t(nn_graph$dist))
}

order_sparse <- function(spm) {
  x <- spm@x
  i <- spm@i
  p <- spm@p

  x_sort <- rep(0, length(x))
  i_sort <- rep(0, length(i))

  n_vertices <- length(p) - 1

  for (v in 1:n_vertices) {
    p_begin <- p[v]
    p_end <- p[v + 1]
    if (p_end - p_begin == 0) {
      next
    }
    pb1 <- p_begin + 1
    x_order <- order(x[pb1:p_end])
    x_sort[pb1:p_end] <- x[x_order + p_begin]
    i_sort[pb1:p_end] <- i[x_order + p_begin]
  }

  list(i = i_sort, p = p, x = x_sort, order = x_order, dims = spm@Dim)
}
