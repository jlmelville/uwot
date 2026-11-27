library(uwot)

test_that("single-block metrics accept column names and indices", {
  set.seed(1337)
  res <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    verbose = FALSE,
    n_threads = 0
  )
  expect_ok_matrix(res)

  set.seed(1337)
  resmli <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = list("euclidean" = 1:4),
    verbose = FALSE,
    n_threads = 0
  )
  expect_equal(resmli, res)

  set.seed(1337)
  resmls <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = list(
      "euclidean" = c(
        "Sepal.Length",
        "Sepal.Width",
        "Petal.Length",
        "Petal.Width"
      )
    ),
    verbose = FALSE,
    n_threads = 0
  )
  expect_equal(resmls, res)
})

test_that("two Euclidean blocks retain separate neighbors for transformation", {
  set.seed(1337)
  jiris10 <- jitter(iris10)
  metric2 <- list(
    "euclidean" = c(1, 2),
    "euclidean" = c("Petal.Length", "Petal.Width")
  )
  reseuc2 <- umap(
    jiris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = metric2,
    verbose = FALSE,
    n_threads = 0,
    ret_nn = TRUE,
    ret_model = TRUE
  )
  expect_ok_matrix(reseuc2$embedding)
  expect_equal(reseuc2$metric, metric2)
  expect_type(reseuc2$nn, "list")
  expect_equal(names(reseuc2$nn), c("euclidean", "euclidean"))
  expect_is_nn(reseuc2$nn[[1]], 10, 4)
  expect_is_nn(reseuc2$nn[[2]], 10, 4)
  expect_ok_matrix(umap_transform(jiris10, reseuc2))
})

test_that("categorical and numeric data can be combined", {
  i10factor <- factor(c(rep("foo", 3), rep("bar", 3), rep("baz", 4)))

  res_y2 <- umap(
    iris10[, -1],
    y = cbind(i10factor, iris$Sepal.Length[1:10]),
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    verbose = FALSE,
    n_threads = 0
  )
  expect_ok_matrix(res_y2)

  nafactor <- as.factor(c(
    levels(iris$Species)[
      c(rep(1, 3), rep(2, 3), rep(3, 3))
    ],
    NA
  ))
  iris10c <- cbind(data.frame(iris10), nafactor)
  rescat <- umap(
    iris10c,
    metric = list("euclidean" = 1:4, "categorical" = "nafactor"),
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    verbose = FALSE,
    n_threads = 0
  )
  expect_ok_matrix(rescat)
})


test_that("mixed metrics can be fitted and transformed", {
  irismixed <- data.frame(iris10, ynum, ynum2, ycat, ycat2)
  resmixed <- umap(
    irismixed,
    n_neighbors = 4,
    n_epochs = 2,
    learning_rate = 0.5,
    min_dist = 0.001,
    init = "rand",
    verbose = FALSE,
    n_threads = 1,
    metric = list(
      "euclidean" = 1:4,
      "euclidean" = 5,
      "cosine" = 6,
      "categorical" = c("ycat", "ycat2")
    )
  )
  expect_ok_matrix(resmixed)

  irismixed <- data.frame(iris10, ynum, ynum2, ycat, ycat2)
  resmixed <- umap(
    irismixed,
    n_neighbors = 4,
    n_epochs = 2,
    learning_rate = 0.5,
    min_dist = 0.001,
    init = "rand",
    verbose = FALSE,
    n_threads = 1,
    metric = list(
      "euclidean" = 1:4,
      "euclidean" = 5,
      "cosine" = 6,
      "categorical" = c("ycat", "ycat2")
    ),
    ret_model = TRUE,
    ret_nn = TRUE
  )
  expect_ok_matrix(resmixed$embedding)
  expect_ok_matrix(umap_transform(irismixed, resmixed, n_threads = 1))
  expect_equal(names(resmixed$nn), c("euclidean", "euclidean", "cosine"))
})

test_that("an integer column can supply a categorical metric", {
  # #20: allow matrix column for categorical
  int_column <- c(1, 2, 3, 4, 4, 4, 2, 1, 2, 1)
  irisic <- cbind(iris10, int_column)
  resic <- umap(
    irisic,
    n_neighbors = 4,
    n_epochs = 2,
    learning_rate = 0.5,
    init = "spca",
    verbose = FALSE,
    n_threads = 0,
    metric = list("euclidean" = 1:4, "categorical" = 5)
  )
  expect_ok_matrix(resic)
})

test_that("HNSW validates the names of feature and target metric blocks", {
  skip_if_not_installed("RcppHNSW")
  blocks <- list(euclidean = 1:2, cosine = 3:4)
  model <- fit_small_model(nn_method = "hnsw", metric = blocks, ret_nn = TRUE)
  for (i in seq_along(blocks)) {
    ref <- fit_small_model(
      iris10[, blocks[[i]]],
      nn_method = "hnsw",
      metric = names(blocks)[i],
      ret_nn = TRUE
    )
    expect_equal(model$nn[[i]], ref$nn[[1]])
  }
  single <- fit_small_model(
    nn_method = "hnsw",
    metric = blocks[1],
    ret_nn = TRUE
  )
  expect_equal(single$nn[[1]], model$nn[[1]])
  supervised <- fit_small_model(
    nn_method = "hnsw",
    y = iris10,
    target_metric = blocks,
    ret_extra = "fgraph"
  )
  expect_true(all(is.finite(supervised$fgraph@x)))
  expect_error(
    fit_small_model(nn_method = "hnsw", metric = list(manhattan = 1L)),
    "bad metric"
  )
  expect_error(
    fit_small_model(
      nn_method = "hnsw",
      y = iris10,
      target_metric = list(manhattan = 1L)
    ),
    "bad target_metric"
  )
})

test_that("PCA skips narrow blocks and preserves full-rank block rotations", {
  blocks <- list(euclidean = 1L, euclidean = 1:2, euclidean = 2:4)
  explicit <- list(
    euclidean = list(1L, pca = NULL),
    euclidean = list(1:2, pca = 2L),
    euclidean = list(2:4, pca = 2L)
  )
  fit <- function(metric, pca = NULL) {
    fit_small_model(
      metric = metric,
      pca = pca,
      pca_method = "svd",
      ret_nn = TRUE
    )
  }
  m <- fit(blocks, pca = 2L)
  ref <- fit(explicit)
  expect_identical(names(m$pca_models), c("2", "3"))
  expect_equal(m$nn, ref$nn)
  project <- function(m) {
    umap_transform(iris10, m, n_epochs = 0, n_threads = 0, ret_extra = "nn")$nn
  }
  expect_equal(project(m), project(ref))
})
