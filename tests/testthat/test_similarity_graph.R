library(uwot)
library(RSpectra)

test_that("graph-only UMAP agrees with similarity_graph and precomputed neighbors", {
  # #96: more convenient way to just get the high dimensional similarity graph

  # Hard way first (by using umap function)
  # allow for no initialization if n_epochs = 0
  # Allowable but returns nothing of value
  expect_warning(
    res <- umap(iris10, n_neighbors = 4, init = NULL, n_epochs = 0),
    "will be returned"
  )
  expect_null(res)

  # More sensibly, return high-dimensional data
  res <- umap(
    iris10,
    n_neighbors = 4,
    init = NULL,
    n_epochs = 0,
    ret_extra = c("fgraph")
  )
  expect_type(res, "list")
  expect_null(res$embedding)
  expect_s4_class(res$fgraph, "sparseMatrix")

  # also a bad idea but maybe ok to extract something from the return value
  # manually
  expect_warning(
    res_with_model <- umap(
      iris10,
      n_neighbors = 4,
      init = NULL,
      n_epochs = 0,
      ret_model = TRUE
    ),
    "will not be valid for transforming"
  )
  expect_type(res_with_model, "list")
  expect_null(res_with_model$embedding)
  # but you cannot use umap_transform with the model
  expect_error(
    umap_transform(iris10, res_with_model),
    "(?i)invalid embedding coordinates"
  )

  # Simpler way (by using similarity_graph function)
  set.seed(42)
  sim_graph <- similarity_graph(iris10, n_neighbors = 4)
  expect_equal(res$fgraph, sim_graph)

  # can return extra data
  set.seed(42)
  sim_graph_extra <- similarity_graph(
    iris10,
    n_neighbors = 4,
    ret_extra = c("sigma", "nn")
  )
  expect_type(sim_graph_extra, "list")
  expect_equal(sim_graph, sim_graph_extra$similarity_graph)
  expect_equal(length(sim_graph_extra$sigma), 10)
  expect_equal(length(sim_graph_extra$rho), 10)
  expect_type(sim_graph_extra$nn$euclidean, "list")

  # can use pre-computed nn instead of data
  sim_graph_nn <- similarity_graph(
    nn_method = sim_graph_extra$nn,
    ret_extra = c("sigma")
  )
  expect_equal(sim_graph, sim_graph_nn$similarity_graph)
  expect_equal(sim_graph_extra$sigma, sim_graph_nn$sigma)
})

test_that("LargeVis returns the similarity graph and bandwidths", {
  # can use largevis for t-SNE-like graph
  sim_graph_largevis <- similarity_graph(
    iris10,
    method = "largevis",
    perplexity = 5,
    ret_extra = c("sigma")
  )
  expect_type(sim_graph_largevis, "list")
  expect_s4_class(sim_graph_largevis$similarity_graph, "sparseMatrix")
  expect_equal(dim(sim_graph_largevis$similarity_graph), c(10, 10))
  expect_equal(length(sim_graph_largevis$sigma), 10)
})

test_that("precomputed neighbors respect the set mix and local connectivity", {
  # specific use case of bbknnR
  fss <- fuzzy_simplicial_set(
    nn = nn,
    set_op_mix_ratio = 0.5,
    local_connectivity = 2
  )

  sim_graph_bbknnR <- similarity_graph(
    nn_method = nn,
    set_op_mix_ratio = 0.5,
    local_connectivity = 2
  )
  expect_equal(sim_graph_bbknnR, fss)
})

test_that("supervision and binary edge weights work with similarity_graph", {
  # supervised
  sim_graphy <- similarity_graph(iris10, n_neighbors = 4, y = ynum)
  expect_s4_class(sim_graphy, "sparseMatrix")

  # binary edge weights
  sim_graphb <- similarity_graph(
    iris10,
    n_neighbors = 4,
    binary_edge_weights = TRUE
  )
  expect_true(all(sim_graphb@x == 1))
})


test_that("optimize graph layout", {
  iris30 <- iris[c(1:10, 51:60, 101:110), ]
  iris30_sim_graph <- similarity_graph(iris30, n_neighbors = 10)
  set.seed(42)
  iris30_opt <- optimize_graph_layout(iris30_sim_graph, X = iris30)

  set.seed(42)
  iris30_umap <- umap(iris30, n_neighbors = 10)
  expect_equal(iris30_opt, iris30_umap)
})

test_that("graph operations expand triangular storage before reading edges", {
  g <- methods::as(
    Matrix::Matrix(exp(-abs(outer(1:16, 1:16, "-"))), sparse = TRUE),
    "generalMatrix"
  )
  diag(g) <- 0
  g <- Matrix::drop0(g)
  layout <- function(x) {
    set.seed(3)
    optimize_graph_layout(x, init = "random", n_epochs = 3, n_sgd_threads = 0)
  }
  expected_layout <- layout(g)
  for (uplo in c("U", "L")) {
    s <- Matrix::forceSymmetric(g, uplo = uplo)
    expect_equal(layout(s), expected_layout)
    for (combine in list(simplicial_set_union, simplicial_set_intersect)) {
      expected <- combine(g, g, n_threads = 0)
      expect_equal(combine(s, s, n_threads = 0), expected)
      expect_equal(combine(g, s, n_threads = 0), expected)
    }
  }
})

test_that("empty graph operands use the finite missing-edge weight", {
  empty <- Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    x = numeric(),
    dims = c(3L, 3L)
  )
  graph <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3),
    j = c(2, 1, 3, 2),
    x = c(1, 1, 0.5, 0.5)
  )
  for (combine in list(simplicial_set_intersect, simplicial_set_union)) {
    left <- combine(empty, graph, n_threads = 0)
    expect_equal(left, combine(graph, empty, n_threads = 0))
    expect_true(all(is.finite(left@x)))
    expect_equal(left > 0, graph > 0)
    expect_length(combine(empty, empty, n_threads = 0)@x, 0)
  }
  nn <- Matrix::sparseMatrix(i = c(2, 3, 1), j = 1:3, x = 1)
  x <- similarity_graph(
    NULL,
    nn_method = nn,
    n_neighbors = 2,
    set_op_mix_ratio = 0,
    y = c(1, 2, 4),
    n_threads = 0
  )
  expect_true(length(x@x) > 0 && all(is.finite(x@x)))
})

test_that("optimizer column indices follow CSC order across empty columns", {
  x <- Matrix::sparseMatrix(
    i = c(1, 3, 1),
    j = c(2, 2, 4),
    x = 1,
    dims = c(4, 5)
  )
  for (graph in list(x, Matrix::t(x))) {
    expected <- Matrix::which(graph != 0, arr.ind = TRUE)[, 2] - 1L
    expect_identical(edge_columns(graph), unname(expected))
  }
})

test_that("graph layout rejects invalid shapes and isolated vertices", {
  g <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = 1)
  expect_error(optimize_graph_layout(as.matrix(g)), "sparse")
  expect_error(optimize_graph_layout(g[, 1, drop = FALSE]), "square")
  expect_error(optimize_graph_layout(Matrix::triu(g)), "symmetric")
  isolated <- Matrix::sparseMatrix(
    i = c(1, 2),
    j = c(2, 1),
    x = 1,
    dims = c(3, 3)
  )
  expect_error(optimize_graph_layout(isolated), "at least one neighbor")
})

test_that("graph combinations require conformable sparse operands and a valid weight", {
  g <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = 1)
  for (combine in list(simplicial_set_union, simplicial_set_intersect)) {
    expect_error(combine(as.matrix(g), g), "x must be a sparse")
    expect_error(combine(g, as.matrix(g)), "y must be a sparse")
    expect_error(combine(g, g[, 1, drop = FALSE]), "identical dimensions")
  }
  expect_error(simplicial_set_intersect(g, g, weight = 2), "weight")
})

test_that("scaling and neighbor defaults are independent for UMAP and LargeVis", {
  x <- as.matrix(iris[1:32, 1:4])
  x[, 1] <- 100 * x[, 1]
  graph <- function(method, ...) {
    similarity_graph(
      x,
      method = method,
      perplexity = 5,
      nn_method = "fnn",
      n_threads = 0,
      ret_extra = "nn",
      ...
    )
  }

  expect_equal(
    graph("umap", scale = "Z"),
    graph("umap", n_neighbors = 15, scale = "Z")
  )
  expect_equal(
    graph("largevis", n_neighbors = 15),
    graph("largevis", scale = "maxabs")
  )
  raw <- graph("largevis", n_neighbors = 15, scale = NULL)
  expect_equal(graph("largevis", scale = NULL), raw)
  expect_false(isTRUE(all.equal(raw$nn, graph("largevis")$nn)))
})
