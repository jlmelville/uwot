library(uwot)
library(RSpectra)

test_that("invalid data, dimensions and neighbor counts are rejected", {
  expect_error(
    umap(X = list(X = "bad", Y = "wrong", nn = "what")),
    "data format"
  )
  expect_error(umap(iris10, n_neighbors = 1, n_threads = 0), "n_neighbors")
  expect_error(umap(iris10, n_neighbors = 15, n_threads = 0), "n_neighbors")
  expect_error(
    umap(iris10, set_op_mix_ratio = 10, n_threads = 0),
    "set_op_mix_ratio"
  )
  expect_error(
    umap(iris10, set_op_mix_ratio = -10, n_threads = 0),
    "set_op_mix_ratio"
  )
  expect_error(
    umap(iris10, local_connectivity = 0.5, n_threads = 0),
    "local_connectivity"
  )
  expect_error(umap(diris10, ret_model = TRUE, n_threads = 0), "models")
  expect_error(umap(dmiris10z, ret_model = TRUE, n_threads = 0), "models")
  expect_error(umap(dmiris10z[, 1:9], n_threads = 0), "distance")
  expect_error(umap(iris[, "Species", drop = FALSE], n_threads = 0), "numeric")
})
test_that("invalid neighbor graphs and unsupported FNN options are rejected", {
  expect_error(umap(iris10, n_threads = 0, nn_method = list()), "precalculated")
  expect_error(
    umap(
      iris10,
      n_threads = 0,
      nn_method = list(
        idx = matrix(1:4, nrow = 2),
        dist = matrix(1:4, nrow = 2)
      )
    ),
    "rows"
  )
  expect_error(
    umap(
      iris10,
      n_threads = 0,
      nn_method = list(idx = matrix(1:40, nrow = 10))
    ),
    "dist"
  )
  expect_error(
    umap(
      iris10,
      n_threads = 0,
      nn_method = list(
        idx = matrix(1:40, nrow = 10),
        dist = matrix(1:4, nrow = 2)
      )
    ),
    "dimensions"
  )
  expect_error(
    umap(
      iris10,
      n_threads = 0,
      n_neighbors = 4,
      nn_method = "fnn",
      metric = "cosine"
    ),
    "FNN"
  )
  expect_error(
    umap(
      iris10,
      n_threads = 0,
      n_neighbors = 4,
      nn_method = "fnn",
      ret_model = TRUE
    ),
    "FNN"
  )
})
test_that("invalid embedding and execution options are rejected", {
  expect_error(lvish(iris10, n_threads = 0, perplexity = 50), "perplexity")
  expect_error(tumap(iris10, n_components = 0), "n_components")
  expect_error(umap(iris10, pca = 1), "'pca' must be >=")
  expect_error(umap(iris10, pca_method = "bad-pca-package"))

  expect_error(
    umap(iris10, n_threads = 0, n_neighbors = 4, y = c(1:9, NA)),
    "numeric y"
  )
  expect_error(
    umap(
      X = NULL,
      n_threads = 0,
      n_neighbors = 4,
      nn_method = nn,
      init = "spca"
    ),
    "spca"
  )
  # add an extra column to nn
  nn5 <- nn
  nn5$idx <- cbind(nn5$idx, rep(100, nrow(nn5$idx)))
  nn5$dist <- cbind(nn5$dist, rep(100.0, nrow(nn5$dist)))
  expect_error(
    umap(X = NULL, n_threads = 0, nn_method = list(nn, nn5)),
    "Invalid neighbor"
  )
  expect_error(umap(iris10, n_threads = 0, pca = 0), "positive integer")

  expect_error(umap(iris10, n_threads = -1), "n_threads")
  expect_error(umap(iris10, n_sgd_threads = -1), "n_sgd_threads")
})

test_that("callback setup remains safe on method validation errors", {
  callback_called <- FALSE
  callback <- function(...) {
    callback_called <<- TRUE
  }

  expect_error(
    uwot:::optimize_layout_r(
      head_embedding = matrix(c(0, 1, 0, 1), ncol = 2),
      tail_embedding = NULL,
      positive_head = integer(),
      positive_tail = integer(),
      positive_ptr = c(0L, 0L, 0L),
      n_epochs = 1L,
      n_head_vertices = 2L,
      n_tail_vertices = 2L,
      epochs_per_sample = numeric(),
      method = "umap",
      method_args = list(),
      initial_alpha = 1,
      opt_args = list(method = "sgd"),
      epoch_callback = callback,
      negative_sample_rate = 5,
      rng_type = "tausworthe",
      batch = FALSE,
      n_threads = 0L,
      grain_size = 1L,
      move_other = TRUE,
      verbose = FALSE
    ),
    "Missing embedding method argument: a"
  )

  expect_false(callback_called)
})

test_that("transform dimensions and initialization types are checked", {
  model <- umap(iris10, n_neighbors = 2, ret_model = TRUE, n_epochs = 2)
  expect_error(umap_transform(iris10[, 1:2], model), "Incorrect dimensions")

  # #42: check init is a matrix or a string; complain otherwise
  expect_error(
    umap(iris10, n_neighbors = 4, init = as.matrix(iris[, 1:3])),
    "(10, 2)"
  )
  expect_error(umap(iris10, n_neighbors = 4, init = iris), "matrix or string")
})

# Don't use data with NA in it
test_that("Detect data with NA in", {
  diris10na <- diris10
  diris10na[1] <- NA
  expect_error(umap(diris10na), "missing", ignore.case = TRUE)

  dmiris10zna <- dmiris10z
  dmiris10zna[2, 1] <- NA
  expect_error(
    umap(dmiris10zna, n_neighbors = 4),
    "missing",
    ignore.case = TRUE
  )

  iris10na <- iris10
  iris10na[1, 1] <- NA
  expect_error(umap(iris10na, n_neighbors = 4), "missing", ignore.case = TRUE)
})

test_that("sparse distances must supply enough neighbors", {
  set.seed(42)
  nnsp10 <- Matrix::drop0(matrix(runif(100), nrow = 10)^2, 0.5)
  expect_error(
    umap(iris10, n_neighbors = 4, nn_method = nnsp10[, -10]),
    "same number"
  )
  expect_error(
    umap(iris10, n_neighbors = 4, nn_method = nnsp10[-10, -10]),
    "unexpected number of rows"
  )

  # give obs 5 0 neighbors
  nnsp10_nbr0 <- nnsp10
  nnsp10_nbr0[, 5] <- 0
  nnsp10_nbr0 <- Matrix::drop0(nnsp10_nbr0)
  expect_error(
    umap(X = NULL, n_neighbors = 4, nn_method = nnsp10_nbr0),
    "at least one neighbor"
  )
})

test_that("embedding and transform initialization are validated", {
  # 76: umap_transform does not validate input sufficiently
  model <- umap(iris[1:10, ], n_neighbors = 4, n_epochs = 0, ret_model = TRUE)
  expect_error(
    trans <- umap_transform(iris[0, ], model = model),
    "Not enough rows"
  )

  # bad min_dist/spread
  expect_error(umap(iris, spread = 1, min_dist = 20), "a, b")

  # n_components too high
  expect_warning(
    umap(
      iris10,
      n_components = 50,
      ret_model = TRUE,
      init = "rand",
      n_neighbors = 4,
      n_epochs = 0
    ),
    "n_components >"
  )

  suppressWarnings(expect_error(
    umap(iris[1:100, ], n_components = 10),
    "Initial data contains NA"
  ))

  for (value in c(NA_real_, Inf, -Inf)) {
    transform_init <- model$embedding[1:5, , drop = FALSE]
    transform_init[1, 1] <- value
    expect_error(
      umap_transform(iris[51:55, ], model, init = transform_init, n_epochs = 0),
      "Initial embedding.*finite"
    )
    invalid_model <- model
    invalid_model$embedding[1, 1] <- value
    expect_error(
      umap_transform(iris[51:55, ], invalid_model, n_epochs = 0),
      "Model embedding.*finite"
    )
  }
})

test_that("large initialization standard deviations produce a warning", {
  # 110: warn if standard deviation of initial input could create small gradients
  expect_warning(
    umap(iris10, init_sdev = 100.0, n_neighbors = 4),
    "embedding standard deviation"
  )
})

test_that("training and projection validate sampler and graph dimensions", {
  for (rate in c(0, -1, Inf, NaN, NA_real_)) {
    expect_error(
      fit_small_model(negative_sample_rate = rate),
      "negative_sample_rate"
    )
  }
  m <- fit_small_model()
  m$negative_sample_rate <- 0
  expect_error(query_neighbors(m), "negative_sample_rate")
  graph <- list(idx = matrix(integer(), 3, 0), dist = matrix(numeric(), 3, 0))
  expect_error(
    similarity_graph(NULL, nn_method = graph, verbose = FALSE),
    "too few"
  )
  expect_error(fit_small_model(y = factor(1:2)), "one row per")
  m <- fit_small_model()
  s <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = c(1, 1), dims = c(3, 2))
  expect_error(
    umap_transform(model = m, nn_method = s, verbose = FALSE),
    "unexpected number of rows"
  )
})

test_that("an epoch callback can run another optimizer in the same R process", {
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script))
  nested_fit <- function() {
    x <- as.matrix(iris[1:10, 1:4])
    fit <- function(...) {
      uwot::umap(
        x,
        init = "random",
        n_neighbors = 4,
        n_threads = 0,
        n_epochs = 3,
        verbose = FALSE,
        ...
      )
    }
    calls <- 0L
    callback <- function(epoch, n_epochs, coords) {
      calls <<- calls + 1L
      if (epoch == 1L) stopifnot(all(is.finite(fit())))
    }
    stopifnot(all(is.finite(fit(epoch_callback = callback))), calls == 3L)
  }
  # Use the same installed build that the parent process is testing.
  lib <- dirname(find.package("uwot"))
  skip_if(
    !file.exists(file.path(lib, "uwot", "Meta", "package.rds")),
    "requires an installed package"
  )
  writeLines(
    c(
      ".libPaths(c(commandArgs(TRUE)[1], .libPaths()))",
      deparse(body(nested_fit))
    ),
    script
  )
  log <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script), shQuote(lib)),
    stdout = TRUE,
    stderr = TRUE,
    timeout = 60
  )
  expect_null(attr(log, "status"), info = paste(log, collapse = "\n"))
})
