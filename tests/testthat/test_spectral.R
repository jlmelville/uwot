library(uwot)
context("Spectral")

test_that("normalized laplacian", {
  # These numbers come from running UMAP Python code:
  # spectral_layout(pairwise_distances(iris.data[0:10, :]))
  # NB:
  # 1. iris data in scikit-learn is currently from UCI repo, which has errors
  #   (although this doesn't affect the first ten entries)
  # 2. eigenvector calculation is not that converged and specifies a starting
  #   vector that we can't supply with either RSpectra or eigen.
  # 3. The eigenvectors are only identical up to a sign, so we take the absolute
  #   values.
  expected_norm_lap <-
    c2y(
      0.7477, -0.1292, -0.03001, 0.02127, -0.563, -0.01149, 0.1402,
      -0.2725, -0.01241, 0.1084, -0.106, -0.5723, 0.2024, -0.3082,
      0.1642, -5.549e-05, -0.04843, -0.1747, 0.1684, 0.6611
    )

  res <- normalized_laplacian_init(x2d(iris[1:10, ]))
  expect_equal(abs(res), abs(expected_norm_lap), tolerance = 1e-2)
})

test_that("laplacian eigenmap", {
  expected_lap_eig <-
    c2y(
      0.3964, -0.2585, -0.297, -0.3923, 0.3905, 0.3581, -0.1268,
      0.2687, -0.356, -0.1954, 0.2775, 0.3298, 0.1282, -0.09545, 0.1503,
      -0.4656, -0.1417, 0.4416, -0.3753, 0.4397
    )

  # Test with distance matrix (simple and symmetric)
  res <- laplacian_eigenmap(x2d(iris[1:10, ]))
  expect_equal(abs(res), abs(expected_lap_eig), tolerance = 1e-4)
})


test_that("1 dimensional output gives a matrix", {
  expect_ok_matrix(spectral_init(V_union, ndim = 1, verbose = FALSE))
  expect_ok_matrix(normalized_laplacian_init(V_union,
                                             ndim = 1,
                                             verbose = FALSE))
  expect_ok_matrix(laplacian_eigenmap(V_union, ndim = 1, verbose = FALSE))
})
