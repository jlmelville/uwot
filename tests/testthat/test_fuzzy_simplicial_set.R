library(uwot)
context("fuzzy simplicial set")

### Various fuzzy set matrices are defined in helper_fuzzy_sets.R

# matrix
# same as
# umap.umap_.fuzzy_simplicial_set(iris10, 4, random_state=42, metric="euclidean")
# as of 0.5 (and even earlier)
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 0, nn = nn
)
expect_equal(res, V_union, tol = 1e-4)

# nnsp0 <- Matrix::sparseMatrix(
#   i = c(0, 4, 7, 9, 1, 2, 3, 9, 1, 2, 3, 6, 2, 3, 8, 9, 0, 4, 6, 7, 0, 4, 5, 7, 2, 3, 6, 7, 0, 4, 7, 9, 1, 2, 3, 8, 1, 2, 3, 9),
#   p = c(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40),
#   x = c(0, 0.141421356237309, 0.173205080756888, 0.469041575982343, 0, 0.3,
#         0.331662479035541, 0.173205080756888, 0.3, 0, 0.244948974278318,
#         0.264575131106459, 0.244948974278318, 0, 0.3, 0.316227766016839,
#         0.141421356237309, 0, 0.458257569495584, 0.223606797749979,
#         0.616441400296898, 0.616441400296898, 0, 0.7, 0.264575131106459,
#         0.33166247903554, 0, 0.424264068711929, 0.173205080756888,
#         0.223606797749979, 0, 0.33166247903554, 0.509901951359279,
#         0.435889894354067, 0.3, 0, 0.173205080756888, 0.316227766016838,
#         0.316227766016839, 0),
#   dims = c(10, 10),
#   index1 = FALSE
# )
# res <- fuzzy_simplicial_set(
#   set_op_mix_ratio = 1, local_connectivity = 1,
#   bandwidth = 1, verbose = FALSE, n_threads = 0, nn = Matrix::drop0(nnsp0)
# )$matrix
# expect_equal(res, V_union, tol = 1e-4)

# mix union + intersection
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0.5, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 0, nn = nn
)
expect_equal(res, V_mix, tol = 1e-4)


# intersection
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 0, nn = nn
)
expect_equal(res, V_intersect, tol = 1e-4)

# Union + local_connectivity
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1.5,
  bandwidth = 1, verbose = FALSE, n_threads = 0, nn = nn
)
expect_equal(res, V_union_local, tol = 1e-4)

# use unique iris nbrs to make comparison with Python UMAP easier
# Union + bandwidth
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1,
  bandwidth = 0.5, verbose = FALSE, n_threads = 0, nn = self_unn4
)
expect_equal(res, V_union_bandwidth, tol = 1e-4)

# intersect + local + bandwidth
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0, local_connectivity = 1.5,
  bandwidth = 0.5, verbose = FALSE, n_threads = 0, nn = self_unn4
)
expect_equal(res, V_intersect_local_bandwidth, tol = 1e-4)

# parallel code path
# matrix
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 1, nn = nn
)
expect_equal(res, V_union, tol = 1e-4)

# mix union + intersection
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0.5, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 1, nn = nn
)
expect_equal(res, V_mix, tol = 1e-4)

# intersection
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 1, nn = nn
)
expect_equal(res, V_intersect, tol = 1e-4)

# Union + local_connectivity
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1.5,
  bandwidth = 1, verbose = FALSE, n_threads = 1, nn = nn
)
expect_equal(res, V_union_local, tol = 1e-4)


# Union + bandwidth
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1,
  bandwidth = 0.5, verbose = FALSE, n_threads = 1, nn = self_unn4
)
expect_equal(res, V_union_bandwidth, tol = 1e-4)

# intersect + local + bandwidth
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0, local_connectivity = 1.5,
  bandwidth = 0.5, verbose = FALSE, n_threads = 1, nn = self_unn4
)
expect_equal(res, V_intersect_local_bandwidth, tol = 1e-4)

# parallel code path
# matrix
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 1, nn = nn
)
expect_equal(res, V_union, tol = 1e-4)

# mix union + intersection
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0.5, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 1, nn = nn
)
expect_equal(res, V_mix, tol = 1e-4)

# intersection
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0, local_connectivity = 1,
  bandwidth = 1, verbose = FALSE, n_threads = 1, nn = nn
)
expect_equal(res, V_intersect, tol = 1e-4)

# Union + local_connectivity
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1.5,
  bandwidth = 1, verbose = FALSE, n_threads = 1, nn = nn
)
expect_equal(res, V_union_local, tol = 1e-4)


# Union + bandwidth
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 1, local_connectivity = 1,
  bandwidth = 0.5, verbose = FALSE, n_threads = 1, nn = self_unn4
)
expect_equal(res, V_union_bandwidth, tol = 1e-4)

# intersect + local + bandwidth
res <- fuzzy_simplicial_set(
  set_op_mix_ratio = 0, local_connectivity = 1.5,
  bandwidth = 0.5, verbose = FALSE, n_threads = 1, nn = self_unn4
)
expect_equal(res, V_intersect_local_bandwidth, tol = 1e-4)
