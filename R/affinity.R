# The UMAP equivalent of perplexity calibration in x2aff. k is continuous rather
# than integral and so is analogous to perplexity.
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
smooth_knn_distances <-
  function(nn,
           n_iter = 64,
           local_connectivity = 1.0,
           bandwidth = 1.0,
           tol = 1e-5,
           min_k_dist_scale = 1e-3,
           verbose = FALSE) {

    nn_dist <- nn$dist
    nn_idx <- nn$idx
    k <- ncol(nn_dist)

    if (verbose) {
      tsmessage("Commencing smooth kNN distance calibration for k = ", formatC(k))
    }

    n <- nrow(nn_dist)
    target <- log2(k) * bandwidth
    rho <- rep(0, n)
    sigma <- rep(0, n)

    mean_distances <- NULL

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
            rho[i] <- non_zero_dists[index]
          }
          else {
            rho[i] <- non_zero_dists[index] + interpolation *
              (non_zero_dists[index + 1] - non_zero_dists[index])
          }
        }
        else {
          rho[i] <- interpolation * non_zero_dists[1]
        }
      } else if (length(non_zero_dists) > 0) {
        rho[i] <- max(non_zero_dists)
      }
      else {
        rho[i] <- 0.0
      }

      for (iter in 1:n_iter) {
        psum <- 0.0
        for (j in 2:ncol(nn_dist)) {
          dist <- max(0, (nn_dist[i, j] - rho[i]))
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
      sigma[i] <- mid

      if (rho[i] > 0.0) {
        sigma[i] <- max(sigma[i], min_k_dist_scale * mean(ith_distances))
      }
      else {
        if (is.null(mean_distances)) {
          mean_distances <- mean(nn_dist)
        }
        sigma[i] <- max(sigma[i], min_k_dist_scale * mean_distances)
      }

      prow <- exp(-(nn_dist[i, ] - rho[i]) / (sigma[i] * bandwidth))
      prow[nn_dist[i, ] - rho[i] <= 0] <- 1
      nn_dist[i, ] <- prow
    }
    P <- Matrix::sparseMatrix(i = rep(1:n, times = k), j = as.vector(nn_idx),
                              x = as.vector(nn_dist))
    Matrix::diag(P) <- 0

    if (verbose) {
      summarize(sigma, "sigma summary")
    }
    list(sigma = sigma, rho = rho, P = Matrix::drop0(P))
  }

# set_op_mix_ratio = between 0 and 1 mixes in fuzzy set intersection
# set to 0 for intersection only
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
