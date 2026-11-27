library(uwot)

iris10prcomp <- prcomp(iris10, retx = TRUE, center = TRUE, scale. = FALSE)

test_that("PCA initialization", {
  iris10_pca_scores <- pca_init(iris10, ndim = 2)
  suppressWarnings(iris10_irlba_scores <- irlba_scores(iris10, ncol = 2))

  expect_equal(
    abs(iris10prcomp$x[, 1:2]),
    abs(iris10_pca_scores),
    ignore_attr = "dimnames"
  )
  expect_equal(
    abs(iris10prcomp$x[, 1:2]),
    abs(iris10_irlba_scores),
    ignore_attr = "dimnames"
  )

  suppressWarnings(iris10_svdr_scores <- irlba_svdr_scores(iris10, ncol = 2))
  expect_equal(
    abs(iris10prcomp$x[, 1:2]),
    abs(iris10_svdr_scores),
    ignore_attr = "dimnames"
  )
})

test_that("1 component initialization works", {
  expect_ok_matrix(pca_init(iris10, ndim = 1), nc = 1)
})

test_that("the PCA shortcut preserves one embedding column through optimization", {
  for (epochs in c(0, 3)) {
    embedding <- umap(
      iris10,
      n_neighbors = 4,
      pca = 2,
      n_components = 1,
      init = "pca",
      n_epochs = epochs,
      n_threads = 0
    )
    expect_ok_matrix(embedding, nc = 1)
  }
})

test_that("PCA returns model data", {
  iris10_pca_scores <- pca_init(iris10, ndim = 2, ret_extra = TRUE)
  expect_equal(
    abs(iris10prcomp$x[, 1:2]),
    abs(iris10_pca_scores$scores),
    ignore_attr = "dimnames"
  )
  expect_equal(
    abs(iris10prcomp$rotation[, 1:2]),
    abs(iris10_pca_scores$rotation),
    ignore_attr = "dimnames"
  )
  expect_equal(
    abs(iris10prcomp$center),
    abs(iris10_pca_scores$center)
  )

  suppressWarnings(
    iris10_irlba_scores <- irlba_scores(iris10, ncol = 2, ret_extra = TRUE)
  )
  expect_equal(
    abs(iris10prcomp$x[, 1:2]),
    abs(iris10_irlba_scores$scores),
    ignore_attr = "dimnames"
  )
  expect_equal(
    abs(iris10prcomp$rotation[, 1:2]),
    abs(iris10_irlba_scores$rotation),
    ignore_attr = "dimnames"
  )
  expect_equal(
    abs(iris10prcomp$center),
    abs(iris10_irlba_scores$center)
  )
})

test_that("logical pca ok", {
  set.seed(1337)
  random_logical <- matrix(rnorm(1000), nrow = 100) > 0.5
  random_int <- random_logical * 1

  expect_equal(
    abs(irlba_scores(random_logical, ncol = 2)),
    abs(irlba_scores(random_int, ncol = 2))
  )
})

test_that("bigstatsr models project training rows in their fitted PCA space", {
  skip_if_not_installed("bigstatsr")
  set.seed(32)
  x <- matrix(rnorm(48 * 20), nrow = 48)
  model <- suppressMessages(umap(
    x,
    nn_method = "annoy",
    n_neighbors = 4,
    pca = 2,
    pca_method = "bigstatsr",
    init = "random",
    n_epochs = 0,
    n_threads = 0,
    ret_model = TRUE,
    verbose = TRUE
  ))
  on.exit(unload_uwot(model, cleanup = FALSE))
  expect_identical(dim(model$pca_models[[1]]$rotation), c(ncol(x), 2L))
  q <- umap_transform(
    x[1:3, ],
    model,
    n_epochs = 0,
    n_threads = 0,
    ret_extra = "nn"
  )
  expect_equal(unname(q$nn[[1]]$idx[, 1]), 1:3)
  expect_ok_matrix(q$embedding, nr = 3)
})
