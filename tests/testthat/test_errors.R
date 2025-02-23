library(uwot)
library(RSpectra)
context("Input validation")

expect_error(umap(X = list(X = "bad", Y = "wrong", nn = "what")), "data format")
expect_error(umap(iris10, n_neighbors = 1, n_threads = 0), "n_neighbors")
expect_error(umap(iris10, n_neighbors = 15, n_threads = 0), "n_neighbors")
expect_error(umap(iris10, set_op_mix_ratio = 10, n_threads = 0), "set_op_mix_ratio")
expect_error(umap(iris10, set_op_mix_ratio = -10, n_threads = 0), "set_op_mix_ratio")
expect_error(umap(iris10, local_connectivity = 0.5, n_threads = 0), "local_connectivity")
expect_error(umap(diris10, ret_model = TRUE, n_threads = 0), "models")
expect_error(umap(dmiris10z, ret_model = TRUE, n_threads = 0), "models")
expect_error(umap(dmiris10z[, 1:9], n_threads = 0), "distance")
expect_error(umap(dmiris10z[, 1:9], n_threads = 0), "distance")
expect_error(umap(iris[, "Species", drop = FALSE], n_threads = 0), "numeric")
expect_error(umap(iris10, n_threads = 0, nn_method = list()), "precalculated")
expect_error(umap(iris10, n_threads = 0, nn_method = list(idx = matrix(1:4, nrow = 2), dist = matrix(1:4, nrow = 2))), "rows")
expect_error(umap(iris10, n_threads = 0, nn_method = list(idx = matrix(1:40, nrow = 10))), "dist")
expect_error(umap(iris10, n_threads = 0, nn_method = list(
  idx = matrix(1:40, nrow = 10),
  dist = matrix(1:4, nrow = 2)
)), "dimensions")
expect_error(umap(iris10, n_threads = 0, n_neighbors = 4, nn_method = "fnn", metric = "cosine"), "FNN")
expect_error(umap(iris10, n_threads = 0, n_neighbors = 4, nn_method = "fnn", ret_model = TRUE), "FNN")
expect_error(lvish(iris10, n_threads = 0, perplexity = 50), "perplexity")
expect_error(tumap(iris10, n_components = 0), "n_components")
expect_error(umap(iris10, pca = 1), "'pca' must be >=")
expect_error(umap(iris10, pca_method = "bad-pca-package"))

expect_error(umap(iris10, n_threads = 0, n_neighbors = 4, y = c(1:9, NA)), "numeric y")
expect_error(umap(
  X = NULL, n_threads = 0, n_neighbors = 4, nn_method = nn,
  init = "spca"
), "spca")
# add an extra column to nn
nn5 <- nn
nn5$idx <- cbind(nn5$idx, rep(100, nrow(nn5$idx)))
nn5$dist <- cbind(nn5$dist, rep(100.0, nrow(nn5$dist)))
expect_error(umap(X = NULL, n_threads = 0, nn_method = list(nn, nn5)), "Invalid neighbor")
expect_error(umap(iris10, n_threads = 0, pca = 0), "positive integer")

expect_error(umap(iris10, n_threads = -1), "n_threads")
expect_error(umap(iris10, n_sgd_threads = -1), "n_sgd_threads")

model <- umap(iris10, n_neighbors = 2, ret_model = TRUE, n_epochs = 2)
expect_error(umap_transform(iris10[, 1:2], model), "Incorrect dimensions")

# #42: check init is a matrix or a string; complain otherwise
expect_error(umap(iris10, n_neighbors = 4, init = as.matrix(iris[, 1:3])), "(10, 2)")
expect_error(umap(iris10, n_neighbors = 4, init = iris), "matrix or string")

# Don't use data with NA in it
test_that("Detect data with NA in", {
  diris10na <- diris10
  diris10na[1] <- NA
  expect_error(umap(diris10na), "missing", ignore.case = TRUE)

  dmiris10zna <- dmiris10z
  dmiris10zna[2, 1] <- NA
  expect_error(umap(dmiris10zna, n_neighbors = 4), "missing", ignore.case = TRUE)

  iris10na <- iris10
  iris10na[1, 1] <- NA
  expect_error(umap(iris10na, n_neighbors = 4), "missing", ignore.case = TRUE)
})

set.seed(42)
nnsp10 <- Matrix::drop0(matrix(runif(100), nrow = 10)^2, 0.5)
expect_error(umap(iris10, n_neighbors = 4, nn_method = nnsp10[, -10]), "same number")
expect_error(umap(iris10, n_neighbors = 4, nn_method = nnsp10[-10, -10]), "unexpected number of rows")

# give obs 5 0 neighbors
nnsp10_nbr0 <- nnsp10
nnsp10_nbr0[, 5] <- 0
nnsp10_nbr0 <- Matrix::drop0(nnsp10_nbr0)
expect_error(umap(X = NULL, n_neighbors = 4, nn_method = nnsp10_nbr0), "at least one neighbor")

# 76: umap_transform does not validate input sufficiently
model <- umap(iris[1:10, ], n_neighbors = 4, n_epochs = 0, ret_model = TRUE)
expect_error(trans <- umap_transform(iris[0, ], model = model), "Not enough rows")

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

# user-supplied intialization should not contain NA
transform_init <- model$embedding[1:5, ]
transform_init[1, 1] <- NA
expect_error(trans <- umap_transform(iris[51:55, ], model = model, init = transform_init), "contains NA")

# model embedding coords should also not contain NA
old11 <- model$embedding[1, 1]
model$embedding[1, 1] <- NA
expect_error(trans <- umap_transform(iris[51:55, ], model = model), "contains NA")
model$embedding[1, 1] <- old11

# 110: warn if standard deviation of initial input could create small gradients
expect_warning(
  umap(iris10, init_sdev = 100.0, n_neighbors = 4),
  "embedding standard deviation"
)
