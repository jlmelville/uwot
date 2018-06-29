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
    Matrix::drop0(set_op_mix_ratio * (X + Matrix::t(X) - XX) + (1 - set_op_mix_ratio) * XX)
  }
}

# Given a set of data X, a neighborhood size, and a measure of distance compute
# the fuzzy simplicial set (here represented as a fuzzy graph in the form of a
# sparse matrix) associated to the data. This is done by locally approximating
# geodesic distance at each point, creating a fuzzy simplicial set for each such
# point, and then combining all the local fuzzy simplicial sets into a global
# one via a fuzzy union
fuzzy_simplicial_set <- function(X, n_neighbors, set_op_mix_ratio = 1.0,
                                 local_connectivity = 1.0, bandwidth = 1.0,
                                 nn_method = "fnn",
                                 n_trees = 50,
                                 search_k = 2 * n_neighbors * n_trees,
                                 n_threads = RcppParallel::defaultNumThreads() / 2,
                                 grain_size = 1000,
                                 verbose = FALSE) {
  nn <- find_nn(X, n_neighbors, method = nn_method, n_trees = n_trees,
                n_threads = n_threads, grain_size = grain_size,
                search_k = search_k, verbose = verbose)

  gc()

  tsmessage("Commencing smooth kNN distance calibration for k = ",
            formatC(n_neighbors))
  affinity_matrix <- smooth_knn_distances_cpp(nn_dist = nn$dist,
                                              nn_idx = nn$idx,
                                              n_iter = 64,
                                              local_connectivity = local_connectivity,
                                              bandwidth = bandwidth,
                                              tol = 1e-5,
                                              min_k_dist_scale = 1e-3,
                                              verbose = verbose)

  fuzzy_set_union(affinity_matrix, set_op_mix_ratio = set_op_mix_ratio)
}

symmetrize <- function(P) {
  0.5 * (P + Matrix::t(P))
}

perplexity_similarities <- function(X, n_neighbors, perplexity,
                                 nn_method = "fnn",
                                 n_trees = 50,
                                 search_k = 2 * n_neighbors * n_trees,
                                 n_threads = RcppParallel::defaultNumThreads() / 2,
                                 grain_size = 1000,
                                 verbose = FALSE) {
  nn <- find_nn(X, n_neighbors, method = nn_method, n_trees = n_trees,
                search_k = search_k, n_threads = n_threads,
                grain_size = grain_size,
                verbose = verbose)

  gc()

  if (n_threads > 0) {
    tsmessage("Commencing perplexity calibration for perplexity = ", formatC(perplexity),
              " k = ", formatC(n_neighbors), " using ", pluralize("thread", n_threads))
    affinity_matrix <- calc_row_probabilities_parallel(nn_dist = nn$dist,
                                                  nn_idx = nn$idx,
                                                  perplexity = perplexity,
                                                  grain_size = grain_size,
                                                  verbose = verbose)
  }
  else {
    tsmessage("Commencing perplexity calibration for perplexity = ", formatC(perplexity),
              " k = ", formatC(n_neighbors))
    affinity_matrix <- calc_row_probabilities_cpp(nn_dist = nn$dist,
                                              nn_idx = nn$idx,
                                              perplexity = perplexity,
                                              verbose = verbose)
  }
  symmetrize(affinity_matrix)
}


# Obsolete pure R function ------------------------------------------------

# The following function is now obsolete because a C++ function is available.
# Remains for testing purposes.
#
# The UMAP equivalent of perplexity calibration.
# Some differences:
# 1. The target value is the log2 of k, not the Shannon entropy associated
# with the desired perplexity.
# 2. Input weights are exponential, rather than Gaussian, with respect to the
# distances. The distances are also centered with respect to the smoothed
# distance to the nearest (non-zero distance) neighbor. A non-integral
# 'local_connectivity' value can result in this shortest distance between an
# interpolated value between two distances.
# 3. The weights are not normalized. Their raw sum is compared to the target
# value.
# 4. Distances beyond the k-nearest neighbors are not used in the calibration.
# The equivalent weights are set to 0.
# 5. Weights associated with distances shorter than the smoothed nearest
# neighbor distance are clamped to 1.
# This code has been converted from the original Python and may not be very
# idiomatic (or vectorizable).
# tol is SMOOTH_K_TOLERANCE in the Python code.
#' @importFrom methods new
#' @import Matrix
smooth_knn_distances <-
  function(nn_dist,
           nn_idx,
           n_iter = 64,
           local_connectivity = 1.0,
           bandwidth = 1.0,
           tol = 1e-5,
           min_k_dist_scale = 1e-3,
           ret_extra = FALSE,
           verbose = FALSE) {

    k <- ncol(nn_dist)
    n <- nrow(nn_dist)

    # In the python code the target is multiplied by the bandwidth, but fuzzy_simplicial_set
    # doesn't pass the user-supplied version on purpose, so it's always 1
    target <- log2(k)

    if (ret_extra) {
      rhos <- rep(0, n)
      sigmas <- rep(0, n)
    }

    mean_distances <- mean(nn_dist)

    progress <- Progress$new(n, display = verbose)
    for (i in 1:n) {
      lo <- 0.0
      hi <- Inf
      mid <- 1.0

      ith_distances <- nn_dist[i, ]
      non_zero_dists <- ith_distances[ith_distances > 0.0]
      if (length(non_zero_dists) >= local_connectivity) {
        index <- floor(local_connectivity)
        interpolation <- local_connectivity - index
        if (index > 0) {
          if (interpolation <= tol) {
            rho <- non_zero_dists[index]
          }
          else {
            rho <- non_zero_dists[index] + interpolation *
              (non_zero_dists[index + 1] - non_zero_dists[index])
          }
        }
        else {
          rho <- interpolation * non_zero_dists[1]
        }
      } else if (length(non_zero_dists) > 0) {
        rho <- max(non_zero_dists)
      }
      else {
        rho <- 0.0
      }

      for (iter in 1:n_iter) {
        psum <- 0.0
        for (j in 2:ncol(nn_dist)) {
          dist <- max(0, (nn_dist[i, j] - rho))
          psum <- psum + exp(-(dist / mid))
        }
        val <- psum

        if (abs(val - target) < tol) {
          break
        }

        if (val > target) {
          hi <- mid
          mid <- (lo + hi) / 2.0
        }
        else {
          lo <- mid
          if (is.infinite(hi)) {
            mid <- mid * 2
          }
          else {
            mid <- (lo + hi) / 2.0
          }
        }
      }
      sigma <- mid

      if (rho > 0.0) {
        sigma <- max(sigma, min_k_dist_scale * mean(ith_distances))
      }
      else {
        sigma <- max(sigma, min_k_dist_scale * mean_distances)
      }

      prow <- exp(-(nn_dist[i, ] - rho) / (sigma * bandwidth))
      prow[nn_dist[i, ] - rho <= 0] <- 1
      nn_dist[i, ] <- prow

      if (ret_extra) {
        rhos[i] <- rho
        sigmas[i] <- sigma
      }

      progress$increment()
    }
    P <- Matrix::sparseMatrix(i = rep(1:n, times = k), j = as.vector(nn_idx),
                              x = as.vector(nn_dist))
    Matrix::diag(P) <- 0
    P <- Matrix::drop0(P)

    if (ret_extra) {
      if (verbose) {
        summarize(sigmas, "sigma summary")
      }
      list(sigma = sigmas, rho = rhos, P = P)
    }
    else {
      P
    }
  }


calc_row_probabilities <- function(nn_dist,
                                   nn_idx,
                                   perplexity,
                                   n_iter = 200,
                                   tol = 1e-5,
                                   ret_extra = FALSE,
                                   verbose = FALSE)
{
  k <- ncol(nn_dist)
  n <- nrow(nn_dist)

  target <- log(perplexity)

  if (ret_extra) {
    sigmas <- rep(0, n)
  }

  progress <- Progress$new(n, display = verbose)
  for (i in 1:n) {

    lo <- 0.0
    hi <- Inf
    mid <- 1.0

    Di <- nn_dist[i, -1] ^ 2

    for (iter in 1:n_iter) {
      sres <- shannon(Di, mid)
      val <- sres$H

      if (abs(val - target) < tol) {
        break
      }

      if (val < target) {
        hi <- mid
        mid <- (lo + hi) / 2.0
      }
      else {
        lo <- mid
        if (is.infinite(hi)) {
          mid <- mid * 2
        }
        else {
          mid <- (lo + hi) / 2.0
        }
      }
    }
    beta <- mid
    sres <- shannon(Di, beta)
    prow <- sres$W / sres$Z
    nn_dist[i, -1] <- prow

    if (ret_extra) {
      sigmas[i] <-  sqrt(1 / beta)
    }

    progress$increment()
  }
  P <- Matrix::sparseMatrix(i = rep(1:n, times = k), j = as.vector(nn_idx),
                            x = as.vector(nn_dist))
  Matrix::diag(P) <- 0
  P <- Matrix::drop0(P)

  if (ret_extra) {
    if (verbose) {
      summarize(sigmas, "sigma summary")
    }
    list(sigma = sigmas, P = P)
  }
  else {
    P
  }
}

shannon <- function(D2, beta) {
  W <- exp(-D2 * beta)
  Z <- sum(W)

  if (Z == 0) {
    H <- 0
  }
  else {
    H <- log(Z) + beta * sum(D2 * W) / Z
  }
  list(
    W = W,
    Z = Z,
    H = H
  )
}
