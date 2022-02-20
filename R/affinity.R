# set_op_mix_ratio = between 0 and 1 mixes in fuzzy set intersection
# set to 0 for intersection only
#' @import Matrix
fuzzy_set_union <- function(X, set_op_mix_ratio = 1) {
  XX <- X * Matrix::t(X)
  if (set_op_mix_ratio == 0) {
    Matrix::drop0(XX)
  }
  else if (set_op_mix_ratio == 1) {
    Matrix::drop0(X + Matrix::t(X) - XX)
  }
  else {
    Matrix::drop0(
      set_op_mix_ratio * (X + Matrix::t(X) - XX) + (1 - set_op_mix_ratio) * XX
    )
  }
}

# Calculate the (asymmetric) affinity matrix based on the nearest neighborhoods
# default target for calibration is the sum of affinities = log2(n_nbrs)
# nn distances should be stored column-wise
smooth_knn <- function(nn_distc,
                       target = NULL,
                       local_connectivity = 1.0,
                       bandwidth = 1.0,
                       n_threads = NULL,
                       grain_size = 1,
                       ret_sigma = FALSE,
                       verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  tsmessage(
    "Commencing smooth kNN distance calibration",
    pluralize("thread", n_threads, " using")
  )
  if (is.null(target)) {
    n_nbrs <- nrow(nn_distc)
    target <- log2(n_nbrs)
  }
  affinity_matrix_res <- smooth_knn_distances_parallel(
    nn_dist = nn_distc,
    target = target,
    n_iter = 64,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
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

# Given nearest neighbor data and a measure of distance compute
# the fuzzy simplicial set (here represented as a fuzzy graph in the form of a
# sparse matrix) associated to the data. This is done by locally approximating
# geodesic distance at each point, creating a fuzzy simplicial set for each such
# point, and then combining all the local fuzzy simplicial sets into a global
# one via a fuzzy union
fuzzy_simplicial_set <- function(nn,
                                 set_op_mix_ratio = 1.0,
                                 local_connectivity = 1.0, bandwidth = 1.0,
                                 ret_sigma = FALSE,
                                 n_threads = NULL,
                                 grain_size = 1,
                                 verbose = FALSE) {
  if (is.null(n_threads)) {
    n_threads <- default_num_threads()
  }
  
  nnt <- nn_graph_t(nn)
  affinity_matrix_res <- smooth_knn(nnt$dist,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
    ret_sigma = ret_sigma,
    n_threads = n_threads,
    grain_size = grain_size,
    verbose = verbose
  )
  affinity_matrix <- affinity_matrix_res$matrix

  affinity_matrix <- nn_to_sparse(nnt$idx, as.vector(affinity_matrix),
    self_nbr = TRUE, by_row = FALSE
  )
  res <- list(matrix = fuzzy_set_union(affinity_matrix, set_op_mix_ratio = set_op_mix_ratio))
  if (ret_sigma) {
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
    affinity_matrix_res <- calc_row_probabilities_parallel(
      nn_dist = nnt$dist,
      perplexity = perplexity,
      ret_sigma = ret_sigma,
      n_threads = n_threads,
      grain_size = grain_size
    )
    if (verbose && affinity_matrix_res$n_failures > 0) {
      tsmessage(affinity_matrix_res$n_failures, " perplexity failures")
    }
    affinity_matrix <- nn_to_sparse(nnt$idx, as.vector(affinity_matrix_res$matrix),
      self_nbr = TRUE, by_row = FALSE
    )
    if (!is.null(affinity_matrix_res$sigma)) {
      sigma <- affinity_matrix_res$sigma
    }
  }
  else {
    # knn kernel
    tsmessage("Using knn graph for input weights with k = ", ncol(nn$idx))
    # Make each row sum to 1, ignoring the self-index
    # i.e. diagonal will be zero
    affinity_matrix <- nn_to_sparse(nn$idx, val = 1 / (ncol(nn$idx) - 1))
    Matrix::diag(affinity_matrix) <- 0
    affinity_matrix <- Matrix::drop0(affinity_matrix)
  }
  res <- list(matrix = symmetrize(affinity_matrix))
  if (ret_sigma && !is.null(sigma)) {
    res$sigma <- sigma
  }
  res
}

# Convert the matrix of NN indices to a sparse asymmetric matrix where each
# edge has a weight of val (scalar or vector)
# return a sparse matrix with dimensions of nrow(nn_idx) x max_nbr_id
nn_to_sparse <- function(nn_idx, val = 1, self_nbr = FALSE,
                         max_nbr_id = NULL, by_row = TRUE) {
  
  if (by_row) {
    n_obs <- nrow(nn_idx)
    n_nbrs <- ncol(nn_idx)
  }
  else {
    n_obs <- ncol(nn_idx)
    n_nbrs <- nrow(nn_idx)
  }
  
  if (is.null(max_nbr_id)) {
    max_nbr_id <- ifelse(self_nbr, n_obs, max(nn_idx))
  }
  
  if (length(val) == 1) {
    xs <- rep(val, n_obs * n_nbrs)
  }
  else {
    xs <- val
  }
  if (by_row) {
    is <- rep(1:n_obs, times = n_nbrs)
  }
  else {
    is <- rep(1:n_obs, each = n_nbrs)
  }
  js <- as.vector(nn_idx)

  dims <- c(n_obs, max_nbr_id)
  res <- Matrix::sparseMatrix(i = is, j = js, x = xs, dims = dims)
  if (self_nbr) {
    Matrix::diag(res) <- 0
    res <- Matrix::drop0(res)
  }
  res
}
