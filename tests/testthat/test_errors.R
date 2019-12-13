library(uwot)
context("Input validation")


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
expect_error(umap(iris10, n_threads = 0, nn_method = list(idx = matrix(1:4, nrow = 2))), "idx")
expect_error(umap(iris10, n_threads = 0, nn_method = list(idx = matrix(1:40, nrow = 10))), "dist")
expect_error(umap(iris10, n_threads = 0, nn_method = list(
  idx = matrix(1:40, nrow = 10),
  dist = matrix(1:4, nrow = 2)
)), "dist")
expect_error(umap(iris10, n_threads = 0, n_neighbors = 4, nn_method = "fnn", metric = "cosine"), "FNN")
expect_error(umap(iris10, n_threads = 0, n_neighbors = 4, nn_method = "fnn", ret_model = TRUE), "FNN")
expect_error(lvish(iris10, n_threads = 0, perplexity = 50), "perplexity")
expect_error(tumap(iris10, n_components = 0), "n_components")
expect_error(umap(iris10, pca = 1), "'pca' must be >=")

expect_error(umap(iris10, n_threads = 0, n_neighbors = 4, y = c(1:9, NA)), "numeric y")
expect_error(umap(
  X = NULL, n_threads = 0, n_neighbors = 4, nn_method = nn,
  init = "spca"
), "spca")
# add an extra column to nn
nn5 <- nn
nn5$idx <- cbind(nn5$idx, rep(100, nrow(nn5$idx)))
nn5$dist <- cbind(nn5$dist, rep(100.0, nrow(nn5$dist)))
expect_error(umap(X = NULL, n_threads = 0, nn_method = list(nn, nn5)), "Precalculated")
expect_error(umap(iris10, n_threads = 0, pca = 0), "positive integer")

expect_error(umap(iris10, n_threads = -1), "n_threads")
expect_error(umap(iris10, n_sgd_threads = -1), "n_sgd_threads")

model <- umap(iris10, n_neighbors = 2, ret_model = TRUE, n_epochs = 2)
expect_error(umap_transform(iris10[, 1:2], model), "Incorrect dimensions")

# #42: check init is a matrix or a string; complain otherwise
expect_error(umap(iris10, n_neighbors = 4, init = as.matrix(iris[, 1:3])), "(10, 2)")
expect_error(umap(iris10, n_neighbors = 4, init = iris), "matrix or string")
