library(uwot)
context("rnndescent integration")

test_that("rnndescent accepts integer-valued doubles from uwot", {
  testthat::skip_if_not_installed("rnndescent")

  model <- umap(
    iris10,
    n_neighbors = 4,
    nn_method = "nndescent",
    nn_args = list(
      n_trees = 2,
      leaf_size = 5,
      max_tree_depth = 20,
      n_iters = 2,
      max_candidates = 5,
      n_search_trees = 1,
      epsilon = 0,
      max_search_fraction = 1
    ),
    n_threads = 1,
    n_build_threads = 1,
    ret_model = TRUE,
    ret_extra = "nn",
    n_epochs = 0,
    verbose = FALSE
  )

  expect_is_nn(model$nn$euclidean)
  expect_type(model$nn$euclidean$idx, "integer")
  expect_true(all(model$nn$euclidean$idx %in% seq_len(nrow(iris10))))
  expect_true(all(is.finite(model$nn$euclidean$dist)))

  query <- umap_transform(
    iris10[1:3, , drop = FALSE],
    model,
    n_threads = 1,
    n_epochs = 0,
    ret_extra = "nn",
    verbose = FALSE
  )

  expect_is_nn(query$nn$euclidean, nr = 3)
  expect_type(query$nn$euclidean$idx, "integer")
  expect_true(all(query$nn$euclidean$idx %in% seq_len(nrow(iris10))))
  expect_true(all(is.finite(query$nn$euclidean$dist)))
})

test_that("rnndescent supports sparse umap2 build and transform", {
  testthat::skip_if_not_installed("rnndescent")

  sparse_iris10 <- Matrix::Matrix(iris10, sparse = TRUE)
  model <- umap2(
    sparse_iris10,
    n_neighbors = 4,
    nn_method = "nndescent",
    n_threads = 1,
    n_build_threads = 1,
    ret_model = TRUE,
    ret_extra = "nn",
    n_epochs = 0,
    verbose = FALSE
  )

  expect_is_nn(model$nn$euclidean)
  expect_type(model$nn$euclidean$idx, "integer")
  expect_true(all(model$nn$euclidean$idx %in% seq_len(nrow(iris10))))
  expect_true(all(is.finite(model$nn$euclidean$dist)))

  query <- umap_transform(
    sparse_iris10[1:3, , drop = FALSE],
    model,
    n_threads = 1,
    n_epochs = 0,
    ret_extra = "nn",
    verbose = FALSE
  )

  expect_is_nn(query$nn$euclidean, nr = 3)
  expect_type(query$nn$euclidean$idx, "integer")
  expect_true(all(query$nn$euclidean$idx %in% seq_len(nrow(iris10))))
  expect_true(all(is.finite(query$nn$euclidean$dist)))
})
