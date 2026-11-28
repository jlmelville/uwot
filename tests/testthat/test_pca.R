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

test_that("irlba PCA accepts logical and integer matrices without warnings", {
  set.seed(1337)
  x <- matrix(rnorm(1000), nrow = 100) > 0.5
  expected <- irlba_scores(x * 1, ncol = 2)
  for (input in list(x, x * 1L)) {
    expect_no_warning(scores <- irlba_scores(input, ncol = 2))
    expect_equal(abs(scores), abs(expected))
  }
})

test_that("distance PCA preserves labels and planar distances", {
  x <- iris10[, 1:2]
  rownames(x) <- paste0("point", seq_len(nrow(x)))
  coords <- umap(
    dist(x),
    n_neighbors = 4,
    init = "pca",
    n_epochs = 0,
    n_threads = 0
  )
  expect_identical(rownames(coords), rownames(x))
  expect_equal(as.vector(dist(coords)), as.vector(dist(x)))
})

for (method in c("svdr", "bigstatsr")) {
  for (center in c(TRUE, FALSE)) {
    test_that(
      paste(method, "projects training rows with pca_center =", center),
      {
        if (method == "bigstatsr") {
          skip_if_not_installed("bigstatsr")
        }
        set.seed(32)
        x <- matrix(rnorm(48 * 20), nrow = 48)
        model <- suppressMessages(umap(
          x,
          nn_method = "annoy",
          n_neighbors = 4,
          pca = 2,
          pca_method = method,
          pca_center = center,
          init = "random",
          n_epochs = 0,
          n_threads = 0,
          ret_model = TRUE,
          verbose = TRUE
        ))
        on.exit(unload_uwot(model, cleanup = FALSE))
        expect_identical(dim(model$pca_models[[1]]$rotation), c(ncol(x), 2L))
        if (center) {
          expect_equal(
            unname(model$pca_models[[1]]$center),
            unname(colMeans(x))
          )
        } else {
          expect_true(all(model$pca_models[[1]]$center == 0))
        }
        q <- umap_transform(
          x[1:3, ],
          model,
          n_epochs = 0,
          n_threads = 0,
          ret_extra = "nn"
        )
        expect_equal(unname(q$nn[[1]]$idx[, 1]), 1:3)
        expect_ok_matrix(q$embedding, nr = 3)
      }
    )
  }
}

test_that("PCA skips Jaccard and still reduces a neighboring Euclidean block", {
  skip_if_not_installed("rnndescent")
  set.seed(11)
  x <- matrix(rbinom(40 * 16, 1, 0.3), nrow = 40)
  fit <- function(metric) {
    fit_small_model(
      x,
      nn_method = "nndescent",
      metric = metric,
      pca = 2,
      ret_nn = TRUE,
      seed = 42
    )
  }
  check_jaccard <- function(nn, binary) {
    idx <- nn$idx[, 2]
    other <- binary[idx, , drop = FALSE]
    expected <- rowSums(binary != other) / pmax(1, rowSums(binary | other))
    expect_gt(sum(expected), 0)
    expect_equal(nn$dist[, 2], expected, tolerance = 1e-6)
  }
  single <- fit("jaccard")
  expect_length(single$pca_models, 0)
  check_jaccard(single$nn$jaccard, x)
  mixed <- fit(list(euclidean = 1:8, jaccard = 9:16))
  expect_equal(names(mixed$pca_models), "1")
  expect_equal(dim(mixed$pca_models[[1]]$rotation), c(8, 2))
  check_jaccard(mixed$nn$jaccard, x[, 9:16])
})
