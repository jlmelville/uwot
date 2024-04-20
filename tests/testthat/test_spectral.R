library(uwot)
library(RSpectra)
context("Spectral")

test_that("1 dimensional output gives a matrix", {
  expect_ok_matrix(spectral_init(V_union, ndim = 1, verbose = FALSE), nc = 1)
  expect_ok_matrix(normalized_laplacian_init(V_union,
    ndim = 1,
    verbose = FALSE
  ), nc = 1)
  expect_ok_matrix(laplacian_eigenmap(V_union, ndim = 1, verbose = FALSE),
    nc = 1
  )
  # 23: ndim was always 2
  expect_ok_matrix(agspectral_init(V_union, n_neg_nbrs = 2, ndim = 1, verbose = FALSE),
    nc = 1
  )
})

test_that("connected components", {
  # Example from doc of scipy.sparse.csgraph.connected_components
  graph <- as(Matrix::drop0(matrix(
    c(
      0, 1, 1, 0, 0,
      0, 0, 1, 0, 0,
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 1,
      0, 0, 0, 0, 0
    ),
    nrow = 5, byrow = TRUE
  )), "generalMatrix")

  cc_res <- connected_components(graph)
  expect_equal(cc_res$n_components, 2)
  expect_equal(cc_res$labels, c(0, 0, 0, 1, 1))

  # Slightly more complicated example validated by running the Python version
  graph100 <- matrix(0, nrow = 10, ncol = 10)
  graph100[cbind(c(2, 6, 7, 8), c(5, 3, 7, 6))] <- 1
  graph100 <- Matrix::drop0(graph100)

  g100_nc <- 7
  g100_labels <- c(0, 1, 2, 3, 1, 2, 4, 2, 5, 6)
  cc_res <- connected_components(graph100)
  expect_equal(cc_res$n_components, g100_nc)
  expect_equal(cc_res$labels, g100_labels)

  # test recursive initialization of components
  sgraph <- graph + Matrix::t(graph)
  expect_ok_matrix(spectral_init(sgraph), nr = 5, nc = 2)
})
