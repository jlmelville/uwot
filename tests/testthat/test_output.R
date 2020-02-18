library(uwot)
context("API output")

set.seed(1337)
# No way to compare with the Python implementation due to differences in
# random number implementations as well as floating point comparison
# and various architecture differences. So we'll just check that the output
# is ok
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "normlaplacian", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)

# Results are repeatable with n_threads = 0 (or 1) and same seed
set.seed(1337)
res2 <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "normlaplacian", verbose = FALSE, n_threads = 0
)
expect_equal(res2, res)



# Distance matrix input
res <- umap(dist(iris10),
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "laplacian", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)

# t-UMAP and cosine metric
res <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, metric = "cosine",
  init = "spectral", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)


# UMAP and cosine metric n_threads = 1 issue #5
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, metric = "cosine",
  init = "spectral", verbose = FALSE, n_threads = 1
)
expect_ok_matrix(res)

# metric = Manhattan
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, metric = "manhattan",
  init = "rand", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, metric = "manhattan",
  init = "spca", verbose = FALSE, n_threads = 1
)
expect_ok_matrix(res)

# init with matrix
iris10_pca <- prcomp(iris10,
  retx = TRUE, center = TRUE,
  scale. = FALSE
)$x[, 1:2]

res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = iris10_pca, verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)
# Ensure that internal C++ code doesn't modify user-supplied initialization
expect_equal(iris10_pca, prcomp(iris10,
  retx = TRUE, center = TRUE,
  scale. = FALSE
)$x[, 1:2])

# return nn
# reset seed here so we can compare output with next test result
set.seed(1337)
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "spca", verbose = FALSE, n_threads = 0,
  ret_nn = TRUE
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)
expect_is(res$nn, "list")
expect_is(res$nn$euclidean, "list")
expect_ok_matrix(res$nn$euclidean$idx, nc = 4)
expect_ok_matrix(res$nn$euclidean$dist, nc = 4)

# Use pre-calculated nn: should be the same as previous result
set.seed(1337)
res_nn <- umap(iris10,
  nn_method = res$nn[[1]], n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "spca", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res_nn)
expect_equal(res_nn, res$embedding)

# X = NULL is ok if passing nn data and rand init
set.seed(1337)
res_nnxn <- umap(
  X = NULL,
  nn_method = nn, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 0
)

# Passing nn list directly is also ok
set.seed(1337)
res_nnl <- umap(iris10,
  nn_method = res$nn, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 0,
  ret_nn = TRUE
)
expect_ok_matrix(res_nnl$embedding)
expect_equal(res_nnl$nn[[1]], res$nn[[1]])
expect_equal(names(res_nnl$nn), "precomputed")
expect_equal(res_nnxn, res_nnl$embedding)

# Use multiple nn data
res_nn2 <- umap(iris10,
  nn_method = list(nn, nn), n_epochs = 2, learning_rate = 0.5,
  min_dist = 0.001,
  init = "spca", verbose = FALSE, n_threads = 0, ret_nn = TRUE
)
expect_ok_matrix(res_nn2$embedding)
expect_equal(names(res_nn2$nn), c("precomputed", "precomputed"))


# lvish and force use of annoy
res <- lvish(iris10,
  perplexity = 4, n_epochs = 2, learning_rate = 0.5, nn_method = "annoy",
  init = "lvrand", verbose = FALSE, n_threads = 1
)
expect_ok_matrix(res)

# lvish with knn
res <- lvish(iris10,
  kernel = "knn", perplexity = 4, n_epochs = 2, learning_rate = 0.5,
  init = "lvrand", verbose = FALSE, n_threads = 1
)
expect_ok_matrix(res)

# return a model
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1,
  ret_model = TRUE
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)

res_test <- umap_transform(iris10, res, n_threads = 1, verbose = FALSE)
expect_ok_matrix(res_test)

# return nn and a model
res <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "rand", verbose = FALSE, n_threads = 1,
  ret_model = TRUE, ret_nn = TRUE
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)

expect_is(res$nn, "list")
expect_is_nn(res$nn[[1]], k = 4)
expect_equal(names(res$nn), "euclidean")
res_test <- umap_transform(iris10, res, n_threads = 0, verbose = FALSE)
expect_ok_matrix(res_test)


# https://github.com/jlmelville/uwot/issues/6
res <- umap(iris10,
  n_components = 1, n_neighbors = 4, n_epochs = 2,
  n_threads = 1, verbose = FALSE
)
expect_ok_matrix(res, nc = 1)


# Supervised
set.seed(1337)
res_y <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  min_dist = 0.001, init = "spca", verbose = FALSE, n_threads = 0,
  y = 1 / (1:10)^2, target_n_neighbors = 2
)
expect_ok_matrix(res_y)

# Repeat using equivalent NN info for y
y_nn <- list(
  idx = matrix(c(
    1, 2,
    2, 3,
    3, 4,
    4, 5,
    5, 6,
    6, 7,
    7, 8,
    8, 9,
    9, 10,
    10, 9
  ), ncol = 2, byrow = TRUE),
  dist = matrix(c(
    0, 0.750000000,
    0, 0.138888896,
    0, 0.048611112,
    0, 0.022500001,
    0, 0.012222221,
    0, 0.007369615,
    0, 0.004783163,
    0, 0.003279321,
    0, 0.002345679,
    0, 0.002345679
  ), ncol = 2, byrow = TRUE)
)

set.seed(1337)
res_ynn <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  min_dist = 0.001, init = "spca", verbose = FALSE, n_threads = 0,
  y = y_nn
)
expect_ok_matrix(res_ynn)
# Should be the same result
expect_equal(res_ynn, res_y)

bin10 <- structure(c(
  0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
  0L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L,
  0L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L
), .Dim = c(10L, 4L))
res <- umap(bin10,
  n_neighbors = 4, metric = "hamming", verbose = FALSE,
  n_threads = 1
)
expect_ok_matrix(res)

# Multiple metrics
set.seed(1337)
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, n_threads = 0,
  metric = list(euclidean = c(1, 2), euclidean = c(3, 4)),
  ret_model = TRUE
)
res_trans <- umap_transform(iris10,
  model = res, verbose = FALSE, n_threads = 0,
  n_epochs = 2
)
expect_ok_matrix(res_trans)

# PCA dimensionality reduction
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, n_threads = 0, pca = 2,
  ret_model = TRUE
)
expect_ok_matrix(res$embedding)
expect_is(res$pca_models, "list")
expect_equal(length(res$pca_models), 1)
expect_ok_matrix(res$pca_models[["1"]]$rotation, nr = 4, nc = 2)
expect_equal(res$pca_models[["1"]]$center, c(4.86, 3.31, 1.45, 0.22),
  check.attributes = FALSE
)

# no centering
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, n_threads = 0, pca = 2,
  pca_center = FALSE, ret_model = TRUE
)
expect_ok_matrix(res$embedding)
expect_is(res$pca_models, "list")
expect_equal(length(res$pca_models), 1)
expect_ok_matrix(res$pca_models[["1"]]$rotation, nr = 4, nc = 2)
expect_null(res$pca_models[["1"]]$center)

res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  metric = list("euclidean" = 1:2, "euclidean" = 3:4),
  init = "spca", verbose = FALSE, n_threads = 0, pca = 2
)
expect_ok_matrix(res)

# Mixed metrics, PCA and transform
set.seed(1337)
ib10 <- cbind(iris10, bin10, bin10)
res <- umap(ib10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, n_threads = 0,
  metric = list(
    euclidean = c(1, 2), hamming = 5:12,
    euclidean = c(3, 4)
  ),
  pca = 2,
  ret_model = TRUE
)
expect_ok_matrix(res$embedding)
expect_is(res$pca_models, "list")
expect_equal(length(res$pca_models), 2)
expect_equal(names(res$pca_models), c("1", "3"))
expect_ok_matrix(res$pca_models[["1"]]$rotation, nr = 2, nc = 2)
expect_equal(res$pca_models[["1"]]$center, c(4.86, 3.31),
  check.attributes = FALSE
)
expect_ok_matrix(res$pca_models[["3"]]$rotation, nr = 2, nc = 2)
expect_equal(res$pca_models[["3"]]$center, c(1.45, 0.22),
  check.attributes = FALSE
)

res_trans <- umap_transform(ib10,
  model = res, verbose = FALSE, n_threads = 0,
  n_epochs = 2
)
expect_ok_matrix(res_trans)

# Override pca command in third block
set.seed(1337)
res <- umap(ib10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, n_threads = 0,
  metric = list(
    euclidean = c(1, 2),
    hamming = 5:8,
    euclidean = list(c(3, 4), pca = NULL)
  ),
  pca = 2,
  ret_model = TRUE
)
expect_ok_matrix(res$embedding)
expect_is(res$pca_models, "list")
expect_equal(length(res$pca_models), 1)
expect_equal(names(res$pca_models), "1")
expect_ok_matrix(res$pca_models[["1"]]$rotation, nr = 2, nc = 2)
expect_equal(res$pca_models[["1"]]$center, c(4.86, 3.31),
  check.attributes = FALSE
)

res_trans <- umap_transform(ib10,
  model = res, verbose = FALSE, n_threads = 0,
  n_epochs = 2
)
expect_ok_matrix(res_trans)

# Turn off PCA centering for binary data
set.seed(1337)
res <- umap(bin10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, n_threads = 0,
  metric = "manhattan", pca = 2, pca_center = FALSE,
  ret_model = TRUE
)
expect_ok_matrix(res$embedding)
expect_is(res$pca_models, "list")
expect_equal(length(res$pca_models), 1)
expect_equal(names(res$pca_models), "1")
expect_ok_matrix(res$pca_models[["1"]]$rotation, nr = 4, nc = 2)
expect_null(res$pca_models[["1"]]$center)

res_trans <- umap_transform(bin10,
  model = res, verbose = FALSE, n_threads = 0,
  n_epochs = 2
)
expect_ok_matrix(res_trans)


# shrunk spectral initialization
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "pca", verbose = FALSE, n_threads = 0, init_sdev = 2
)
expect_ok_matrix(res)

res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "laplacian", verbose = FALSE, n_threads = 0,
  init_sdev = 0.1
)
expect_ok_matrix(res)

res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spectral", verbose = FALSE, n_threads = 0,
  init_sdev = 5
)
expect_ok_matrix(res)


# umap transform when test datset size > train dataset size
set.seed(1337)
res <- umap(iris10[1:4, ],
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, ret_model = TRUE
)
expect_is(res, "list")
expect_ok_matrix(res$embedding, nr = 4)

res_test <- umap_transform(iris10[5:10, ], res, verbose = FALSE, n_epochs = 10)
expect_ok_matrix(res_test, nr = 6)
# 31: ensure we store the ndim for single-metric models
expect_equal(res$metric$euclidean$ndim, 4)

# taus88 prng
res <- umap(iris10,
  pcg_rand = FALSE,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spectral", verbose = FALSE, n_threads = 0,
  init_sdev = 5
)
expect_ok_matrix(res)

# https://github.com/jlmelville/uwot/issues/39
res <- umap(iris10, n_neighbors = 4, n_threads = 0.5)
expect_ok_matrix(res)
res <- umap(iris10, n_neighbors = 4, n_threads = 1.5)
expect_ok_matrix(res)
res <- umap(iris10, n_neighbors = 4, n_sgd_threads = 0.5)
expect_ok_matrix(res)
res <- umap(iris10, n_neighbors = 4, n_sgd_threads = 1.5)
expect_ok_matrix(res)

# https://github.com/jlmelville/uwot/issues/47
# return fuzzy graph
res <- umap(iris10,
            n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
            init = "spca", verbose = FALSE, n_threads = 0,
            ret_extra = c("fgraph")
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)
expect_is(res$fgraph, "Matrix")

res <- tumap(iris10,
            n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
            init = "spca", verbose = FALSE, n_threads = 0,
            ret_extra = c("fgraph")
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)
expect_is(res$fgraph, "Matrix")

# param is ret_P and returned value is P in lvish
res <- lvish(iris10,
             perplexity = 4, n_epochs = 2, learning_rate = 0.5,
             init = "spca", verbose = FALSE, n_threads = 0,
             ret_extra = c("P")
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)
expect_is(res$P, "Matrix")
