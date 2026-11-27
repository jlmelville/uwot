library(uwot)

test_that("standardization centers and rescales each feature", {
  iris10_sd <- apply(iris10, 2, sd)
  iris10_mean <- apply(iris10, 2, mean)

  iris10_none <- scale_input(iris10, scale_type = FALSE)
  expect_equal(apply(iris10_none, 2, sd), iris10_sd)
  expect_equal(apply(iris10_none, 2, mean), iris10_mean)

  iris10_scale <- scale_input(iris10, scale_type = TRUE)
  expect_equal(unname(apply(iris10_scale, 2, sd)), rep(1, 4))
  expect_equal(
    unname(apply(iris10_scale, 2, mean)),
    rep(0, 4)
  )

  # "scale" and "z" and TRUE are synonyms
  expect_equal(scale_input(iris10, scale_type = "scale"), iris10_scale)
  expect_equal(scale_input(iris10, scale_type = "Z"), iris10_scale)
})

test_that("range and maximum-absolute scaling use the requested bounds", {
  iris10_maxabs <- scale_input(iris10, scale_type = "maxabs")
  expect_equal(
    unname(apply(iris10_maxabs, 2, mean)),
    rep(0, 4)
  )
  expect_equal(max(abs(iris10_maxabs)), 1)

  iris10_range <- scale_input(iris10, scale_type = "range")
  expect_equal(max(iris10_range), 1)
  expect_equal(min(iris10_range), 0)

  iris10_colrange <- scale_input(iris10, scale_type = "colrange")
  expect_equal(
    unname(apply(iris10_colrange, 2, max)),
    rep(1, 4)
  )
  expect_equal(
    unname(apply(iris10_colrange, 2, min)),
    rep(0, 4)
  )
})

test_that("scaling applied outside umap should not appear in the model", {
  iris10s <- scale(iris10)
  iris10s_umap0 <-
    umap(
      iris10s,
      n_neighbors = 4,
      n_epochs = 0,
      init = "rand",
      ret_model = TRUE,
      scale = FALSE
    )
  expect_null(iris10s_umap0$scale_info)

  iris10s_umap0s <-
    umap(
      iris10s,
      n_neighbors = 4,
      n_epochs = 0,
      init = "rand",
      ret_model = TRUE,
      scale = TRUE
    )
  expect_equal(
    names(iris10s_umap0s$scale_info),
    c("scaled:center", "scaled:scale", "scaled:nzvcols")
  )
})

test_that("metric blocks select the original features across scaling", {
  x <- cbind(constant = 1, trend = 0:3, zig = c(0, 10, 1, 11))
  neighbors <- function(cols) {
    umap(
      x,
      metric = list(euclidean = cols),
      scale = TRUE,
      nn_method = "fnn",
      n_neighbors = 2,
      n_components = 1,
      init = "random",
      n_epochs = 0,
      n_threads = 0,
      ret_nn = TRUE
    )$nn[[1]]
  }
  expect_equal(neighbors(2L), neighbors("trend"))
  expect_equal(neighbors(list(3L, pca = NULL)), neighbors("zig"))
  expect_equal(neighbors(2L)$idx[1, 2], 2L)
  expect_error(neighbors(1:2), "zero-variance column")
  expect_error(scale_input(x[, 1, drop = FALSE], TRUE), "zero variance")

  m <- fit_small_model(x, scale = TRUE, metric = list(euclidean = 2L))
  ref <- fit_small_model(x, scale = TRUE, metric = list(euclidean = "trend"))
  expect_equal(query_neighbors(m, x), query_neighbors(ref, x))
})

test_that("scaling preserves the dimensions of individual features and queries", {
  for (pc in list(NULL, 2L)) {
    m <- fit_small_model(scale = TRUE, pca = pc, pca_method = "svd")
    one <- query_neighbors(m, iris10[1, , drop = FALSE])
    two <- query_neighbors(m, iris10[1:2, , drop = FALSE])
    expect_equal(one$idx, two$idx[1, , drop = FALSE])
    expect_equal(one$dist, two$dist[1, , drop = FALSE])
  }
  x <- cbind(constant = 1, varying = seq_len(10))
  m <- fit_small_model(x, scale = TRUE, n_components = 1)
  expect_equal(query_neighbors(m, x[5, , drop = FALSE])$idx[1, 1], 5L)
})

test_that("column range scaling ignores training constants in every query", {
  x <- cbind(constant = 7, iris10)
  m <- fit_small_model(x, scale = "colrange", ret_nn = TRUE)
  ref <- fit_small_model(iris10, scale = "colrange", ret_nn = TRUE)
  expect_equal(m$nn, ref$nn)
  query <- x[1:2, , drop = FALSE]
  query[, "constant"] <- c(-100, 100)
  expected <- query_neighbors(ref, iris10[1:2, , drop = FALSE])
  expect_equal(query_neighbors(m, query), expected)
  file <- tempfile()
  on.exit(unlink(file))
  save_uwot(m, file, unload = TRUE)
  restored <- load_uwot(file)
  on.exit(unload_uwot(restored), add = TRUE)
  expect_equal(query_neighbors(restored, query), expected)

  constants <- scale_input(
    x[, "constant", drop = FALSE],
    "colrange",
    ret_model = TRUE
  )
  expect_true(all(constants == 0))
  expect_equal(
    as.numeric(apply_scaling(matrix(99), attr_to_scale_info(constants))),
    0
  )
})

test_that("global scaling reports zero training scale factors", {
  flat <- matrix(3, 10, 2)
  columns <- cbind(rep(1, 10), rep(2, 10))
  expect_error(
    fit_small_model(flat, scale = "range"),
    "training scale factor is zero"
  )
  expect_error(
    fit_small_model(columns, scale = "maxabs"),
    "training scale factor is zero"
  )
  expect_equal(range(scale_input(columns, "range")), c(0, 1))
})
