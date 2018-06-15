library(uwot)
context("fuzzy simplicial set")

### Various fuzzy set matrices are defined in helper_fuzzy_sets.R

# matrix
res <- fuzzy_simplicial_set(iris10, n_neighbors = 4, set_op_mix_ratio = 1, local_connectivity = 1,
                            bandwidth = 1, nn_method = "fnn", verbose = FALSE)
expect_equal(res, V_union, tol = 1e-4)

# dist
res <- fuzzy_simplicial_set(diris10, n_neighbors = 4, set_op_mix_ratio = 1, local_connectivity = 1,
                            bandwidth = 1, nn_method = "fnn", verbose = FALSE)
expect_equal(res, V_union, tol = 1e-4)

# sparse distance matrix
res <- fuzzy_simplicial_set(dmiris10z, n_neighbors = 4, set_op_mix_ratio = 1, local_connectivity = 1,
                            bandwidth = 1, nn_method = "fnn", verbose = FALSE)
expect_equal(res, V_union, tol = 1e-4)

# mix union + intersection
res <- fuzzy_simplicial_set(iris10, n_neighbors = 4, set_op_mix_ratio = 0.5, local_connectivity = 1,
                            bandwidth = 1, nn_method = "fnn", verbose = FALSE)
expect_equal(res, V_mix, tol = 1e-4)

# intersection
res <- fuzzy_simplicial_set(iris10, n_neighbors = 4, set_op_mix_ratio = 0, local_connectivity = 1,
                            bandwidth = 1, nn_method = "fnn", verbose = FALSE)
expect_equal(res, V_intersect, tol = 1e-4)

# Union + local_connectivity
res <- fuzzy_simplicial_set(iris10, n_neighbors = 4, set_op_mix_ratio = 1, local_connectivity = 1.5,
                            bandwidth = 1, nn_method = "fnn", verbose = FALSE)
expect_equal(res, V_union_local, tol = 1e-4)


# Union + bandwidth
res <- fuzzy_simplicial_set(iris10, n_neighbors = 4, set_op_mix_ratio = 1, local_connectivity = 1,
                            bandwidth = 0.5, nn_method = "fnn", verbose = FALSE)
expect_equal(res, V_union_bandwidth, tol = 1e-4)

# intersect + local + bandwidth
res <- fuzzy_simplicial_set(iris10, n_neighbors = 4, set_op_mix_ratio = 0, local_connectivity = 1.5,
                            bandwidth = 0.5, nn_method = "fnn", verbose = FALSE)
expect_equal(res, V_intersect_local_bandwidth, tol = 1e-4)
