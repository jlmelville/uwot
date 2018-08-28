library(uwot)
context("API output")

set.seed(1337)
# No way to compare with the Python implementation due to differences in
# random number implementations as well as floating point comparison
# and various architecture differences. So we'll just check that the output
# is ok
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
  init = "normlaplacian", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)

# Distance matrix input
res <- umap(dist(iris10),
  n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
  init = "laplacian", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)

# t-UMAP and cosine metric
res <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, alpha = 0.5, metric = "cosine",
  init = "spectral", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)


# UMAP and cosine metric n_threads = 1 issue #5
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, alpha = 0.5, metric = "cosine",
  init = "spectral", verbose = FALSE, n_threads = 1
)
expect_ok_matrix(res)

# metric = Manhattan
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, alpha = 0.5, metric = "manhattan",
  init = "rand", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, alpha = 0.5, metric = "manhattan",
  init = "spca", verbose = FALSE, n_threads = 1
)
expect_ok_matrix(res)

# lvish and force use of annoy
res <- lvish(iris10,
  perplexity = 4, n_epochs = 2, alpha = 0.5, nn_method = "annoy",
  init = "lvrand", verbose = FALSE, n_threads = 1
)
expect_ok_matrix(res)

# lvish with knn
res <- lvish(iris10,
  kernel = "knn", perplexity = 4, n_epochs = 2, alpha = 0.5,
  init = "lvrand", verbose = FALSE, n_threads = 1
)
expect_ok_matrix(res)

# return a model
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1,
  ret_model = TRUE
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)

res_test <- umap_transform(iris10, res, n_threads = 1, verbose = FALSE)
expect_ok_matrix(res_test)


res <- umap(iris10, n_components = 1, n_neighbors = 4, n_epochs = 2,
            n_threads = 1, verbose = FALSE)
expect_ok_matrix(res)

