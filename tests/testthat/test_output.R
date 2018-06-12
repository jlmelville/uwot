library(uwot)
context("UMAP output")

# No way to compare with the Python implementation due to differences in
# random number implementations as well as floating point comparison.
# These tests characterize the C++ output, which also matched a pure R
# implementation.

norm_vec_out <- c(
  1.38167422, -0.56341620, -0.52414563, -0.44066463, 0.35588609, 0.92339921,
  -0.10439442, 1.27274590, -0.29202394, -0.78932333, 0.50988506, -0.80345216,
  2.37696125,  0.61399186, -1.75754661,  0.03745652, 0.60281106, -0.39475789,
  0.16222497, -1.09888769)
set.seed(1337)
res <- umap(iris10, n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
            init = "normlaplacian", verbose = FALSE)
expect_equal(res, c2y(norm_vec_out), tol = 1e-6)

# Distance matrix input
set.seed(1337)
res <- umap(dist(iris10), n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
            init = "normlaplacian", verbose = FALSE)
expect_equal(res, c2y(norm_vec_out), tol = 1e-6)

# Normalized Laplacian + Noise initializaton
set.seed(1337)
res <- umap(iris10, n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
            init = "spectral", verbose = FALSE)
expect_equal(res,
             c2y(c(
                6.2398524, -3.5898787, -5.2273663, -5.2708939,  7.3886341,
                7.7185850, -1.7308497,  5.0974633, -4.7344445, -4.2000233,
               -0.9458299, -8.5785094,  3.8055032,  2.6296196,  0.1233900,
               -1.1830738, 10.0006417,  0.3684023,  1.5181026, -8.3615143
             )), tol = 1e-6)


