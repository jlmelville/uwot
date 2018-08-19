library(uwot)
context("Smooth kNN distances")

nn_8 <- find_nn(iris10, k = 8)
res <- smooth_knn_distances(nn_8$dist, nn_8$idx, ret_extra = TRUE)
expect_equal(res$sigma, c(
  0.2567215, 0.22098923, 0.08285332, 0.09981823,
  0.28608322, 0.17873764, 0.15968704, 0.17134094,
  0.25434113, 0.19572449
))
expect_equal(res$rho, c(
  0.14142136, 0.17320508, 0.24494897, 0.24494897,
  0.14142136, 0.6164414, 0.26457513, 0.17320508,
  0.3, 0.17320508
))


nn_4 <- find_nn(iris10, k = 4)
res <- smooth_knn_distances(nn_4$dist, nn_4$idx, ret_extra = TRUE)
expect_equal(res$sigma, c(
  0.17993927, 0.20488739, 0.0493803, 0.09060478,
  0.24940491, 0.00390625, 0.15367126, 0.13551712,
  0.24542618, 0.20633698
))
expect_equal(res$rho, c(
  0.14142136, 0.17320508, 0.24494897, 0.24494897,
  0.14142136, 0.6164414, 0.26457513, 0.17320508, 0.3,
  0.17320508
))

# distance matrix
res <- smooth_knn_distances(nn_4$dist, nn_4$idx, ret_extra = TRUE)
expect_equal(res$sigma, c(
  0.17993927, 0.20488739, 0.0493803, 0.09060478,
  0.24940491, 0.00390625, 0.15367126, 0.13551712,
  0.24542618, 0.20633698
))
expect_equal(res$rho, c(
  0.14142136, 0.17320508, 0.24494897, 0.24494897,
  0.14142136, 0.6164414, 0.26457513, 0.17320508, 0.3,
  0.17320508
))

### Various fuzzy set matrices are defined in helper_fuzzy_sets.R
# unsymmetrized fuzzy set
res$P <- nn_to_sparse(nn_4$idx, as.vector(res$P), self_nbr = TRUE)
expect_equal(res$P, V_asymm, tol = 1e-4)

# Fuzzy Set Union
expect_equal(fuzzy_set_union(res$P), V_union, tol = 1e-4)

# mix intersection with union
expect_equal(fuzzy_set_union(res$P, set_op_mix_ratio = 0.5), V_mix,
  tol = 1e-4
)

# intersection
expect_equal(fuzzy_set_union(res$P, set_op_mix_ratio = 0), V_intersect,
  tol = 1e-4
)

# local connectivity
res <- smooth_knn_distances(nn_4$dist, nn_4$idx,
  local_connectivity = 1.5,
  ret_extra = TRUE
)
res$P <- nn_to_sparse(nn_4$idx, as.vector(res$P), self_nbr = TRUE)
expect_equal(res$P, V_asymm_local, tol = 1e-4)

### C++ tests
res_cpp_conn1 <- smooth_knn_distances_parallel(nn_4$dist, nn_4$idx,
  n_iter = 64, local_connectivity = 1.0,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  parallelize = FALSE, verbose = FALSE
)
expect_equal(nn_to_sparse(nn_4$idx, as.vector(res_cpp_conn1),
  self_nbr = TRUE
), V_asymm, tol = 1e-4)

res_cpp_conn1.5 <- smooth_knn_distances_parallel(nn_4$dist, nn_4$idx,
  n_iter = 64, local_connectivity = 1.5,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  parallelize = FALSE, verbose = FALSE
)
expect_equal(nn_to_sparse(nn_4$idx, as.vector(res_cpp_conn1.5),
  self_nbr = TRUE
), V_asymm_local, tol = 1e-4)


RcppParallel::setThreadOptions(numThreads = 1)
res_cpp_conn1 <- smooth_knn_distances_parallel(nn_4$dist, nn_4$idx,
  n_iter = 64, local_connectivity = 1.0,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  grain_size = 1, verbose = FALSE
)
expect_equal(nn_to_sparse(nn_4$idx, as.vector(res_cpp_conn1),
  self_nbr = TRUE
), V_asymm, tol = 1e-4)

res_cpp_conn1.5 <- smooth_knn_distances_parallel(nn_4$dist, nn_4$idx,
  n_iter = 64, local_connectivity = 1.5,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  grain_size = 1, verbose = FALSE
)
expect_equal(nn_to_sparse(nn_4$idx, as.vector(res_cpp_conn1.5),
  self_nbr = TRUE
), V_asymm_local, tol = 1e-4)


# Test cross-distances
V_asymm_local_cross <- V_asymm_local
diag(V_asymm_local_cross) <- 1
V_asymm_local_cross <- cbind(
  V_asymm_local_cross,
  matrix(0, nrow = 10, ncol = 2)
)

res_cpp_conn1.5_cross <- smooth_knn_distances_parallel(nn_4$dist, nn_4$idx,
  n_iter = 64, local_connectivity = 1.5,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  parallelize = FALSE, verbose = FALSE
)
expect_equal(nn_to_sparse(
  nn_4$idx, as.vector(res_cpp_conn1.5_cross),
  self_nbr = FALSE, max_nbr_id = 12
),
V_asymm_local_cross,
tol = 1e-4
)

res_cpp_conn1.5_cross <- smooth_knn_distances_parallel(nn_4$dist, nn_4$idx,
  n_iter = 64, local_connectivity = 1.5,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3, verbose = FALSE
)
expect_equal(nn_to_sparse(
  nn_4$idx, as.vector(res_cpp_conn1.5_cross),
  self_nbr = FALSE, max_nbr_id = 12
),
V_asymm_local_cross,
tol = 1e-4
)
