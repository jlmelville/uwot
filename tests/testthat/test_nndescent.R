library(uwot)

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

test_that("sparse default uses NN-Descent and retains only its search index", {
  skip_if_not_installed("rnndescent")
  x <- Matrix::Matrix(iris10, sparse = TRUE)
  fit <- function(...) {
    umap2(
      x,
      n_neighbors = 4,
      init = "random",
      n_epochs = 0,
      ret_model = TRUE,
      n_threads = 0,
      seed = 42,
      verbose = FALSE,
      ...
    )
  }
  auto <- fit()
  explicit <- fit(nn_method = "nndescent", ret_nn = TRUE)
  expect_identical(auto$nn_method, "nndescent")
  expect_null(auto$nn_index$ann$graph)
  expect_null(explicit$nn_index$ann$graph)
  expect_equal(dim(explicit$nn[[1]]$idx), c(nrow(x), 4))
  set.seed(42)
  q <- query_neighbors(explicit, x)
  set.seed(42)
  expect_equal(query_neighbors(auto, x), q)
})

test_that("each NN-Descent feature block remains searchable", {
  skip_if_not_installed("rnndescent")
  m <- fit_small_model(
    nn_method = "nndescent",
    metric = list(euclidean = 1:2, euclidean = 3:4)
  )
  query <- function(m) {
    umap_transform(
      iris10,
      m,
      n_epochs = 0,
      n_threads = 0,
      seed = 42,
      ret_extra = "nn",
      verbose = FALSE
    )$nn
  }
  before <- query(m)
  expect_length(before, 2)
  for (i in 1:2) {
    cols <- m$metric[[i]]
    d <- as.matrix(dist(iris10[, cols]))
    idx <- before[[i]]$idx
    expected <- matrix(d[cbind(as.vector(row(idx)), as.vector(idx))], nrow(idx))
    expect_equal(unname(before[[i]]$dist), expected, tolerance = 1e-5)
  }
  archive <- tempfile()
  on.exit(unlink(archive))
  save_uwot(m, archive, unload = TRUE)
  m <- load_uwot(archive)
  on.exit(unload_uwot(m), add = TRUE, after = FALSE)
  expect_equal(query(m), before)
})
