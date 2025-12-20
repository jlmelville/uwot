library(uwot)
library(RSpectra)
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
res <- umap(stats::dist(iris10),
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

nn_index <- list(index = res$nn[[1]]$idx, distance = res$nn[[1]]$dist)
set.seed(1337)
res_nn_index <- umap(iris10,
  nn_method = nn_index, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "spca", verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res_nn_index)
expect_equal(res_nn_index, res$embedding)

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

# Passing nn list directly and return a model
set.seed(1337)
res_nnl <- umap(iris10,
  nn_method = res$nn, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 0,
  ret_nn = TRUE, ret_model = TRUE
)
expect_ok_matrix(res_nnl$embedding)
expect_equal(res_nnl$nn[[1]], res$nn[[1]])
expect_equal(names(res_nnl$nn), "precomputed")
expect_equal(res_nnxn, res_nnl$embedding)
expect_equal(res_nnl$num_precomputed_nns, 1)

# Passing nn list directly and return a model and set X to NULL
set.seed(1337)
res_nnl <- umap(
  X = NULL,
  nn_method = res$nn, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 0,
  ret_nn = TRUE, ret_model = TRUE
)
expect_ok_matrix(res_nnl$embedding)
expect_equal(res_nnl$nn[[1]], res$nn[[1]])
expect_equal(names(res_nnl$nn), "precomputed")
expect_equal(res_nnxn, res_nnl$embedding)
expect_equal(res_nnl$num_precomputed_nns, 1)

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
  init = "lvrand", verbose = FALSE, n_threads = 1, ret_extra = c("sigma")
)
expect_ok_matrix(res$embedding)
expect_equal(res$sigma, sqrt(c(0.3039, 0.2063, 0.09489, 0.08811, 0.3091, 0.6789, 0.1743, 0.1686, 0.3445, 0.1671)), tol = 1e-4)

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
# #95: export min_dist and spread in returned model
expect_equal(res$min_dist, 0.001)
expect_equal(res$spread, 1)

resab <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, a = 1, b = 0.9,
  init = "rand", verbose = FALSE, n_threads = 1,
  ret_model = TRUE
)
expect_equal(resab$a, 1)
expect_equal(resab$b, 0.9)
# min_dist and spread in returned model are NULL if a and b are set
expect_null(resab$min_dist)
expect_null(resab$spread)


res_test <- umap_transform(iris10, res, n_threads = 1, verbose = FALSE)
expect_ok_matrix(res_test)


# test we can use 0 epochs
res_test0 <- umap_transform(iris10, res, n_epochs = 0, n_threads = 1, verbose = FALSE)
expect_ok_matrix(res_test)
expect_equal(dim(res_test0), c(10, 2))

# return nn and a model
set.seed(42)
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

# test sparse nn matrix exactly the same as knn graph with explicit 0s for
# self neighbors
sparse_nbr_matrix0 <- Matrix::sparseMatrix(
  i = c(
    0, 4, 7, 9, 1, 2, 3, 9, 1, 2, 3, 6, 2, 3, 8, 9, 0, 4, 6, 7, 0, 4, 5, 7, 2, 3, 6, 7, 0, 4, 7, 9, 1, 2, 3, 8, 1, 2, 3, 9
  ),
  p = c(
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40
  ),
  x = c(
    0, 0.141421228647, 0.173204988241, 0.469041585922,
    0, 0.300000220537, 0.331662625074, 0.173205047846, 0.300000220537,
    0, 0.244949027896, 0.264575153589, 0.244949027896,
    0, 0.2999997437, 0.316227942705, 0.141421228647,
    0, 0.458257555962, 0.223606646061, 0.616441547871, 0.616441547871,
    0, 0.700000047684, 0.264575153589, 0.331662654877,
    0, 0.424264162779, 0.173204988241, 0.223606646061,
    0, 0.331662625074, 0.509901940823, 0.435889661312, 0.2999997437,
    0, 0.173205047846, 0.31622800231, 0.316227942705,
    0
  ),
  index1 = FALSE
)
set.seed(42)
res_spnn0 <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "rand", verbose = FALSE, n_threads = 1,
  nn_method = sparse_nbr_matrix0, ret_nn = TRUE
)
expect_is(res, "list")
expect_ok_matrix(res_spnn0$embedding)
# should get same results as with internal nn calculation
expect_equal(res_spnn0$embedding, res$embedding)

sparse_nbr_matrix0_with_names <- sparse_nbr_matrix0
row.names(sparse_nbr_matrix0_with_names) <- row.names(iris10)
colnames(sparse_nbr_matrix0_with_names) <- row.names(iris10)

expect_equal(res_spnn0$nn$euclidean, sparse_nbr_matrix0_with_names)

# sparse neighbor matrix without explicit zeros
sparse_nbr_matrix <- Matrix::sparseMatrix(
  i = c(
    4, 7, 9, 2, 3, 9, 1, 3, 6, 2, 8, 9, 0, 6, 7, 0, 4, 7, 2, 3, 7, 0, 4, 9, 1, 2, 3, 1, 2, 3
  ),
  p = c(
    0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30
  ),
  x = c(
    0.141421228647, 0.173204988241, 0.469041585922,
    0.300000220537, 0.331662625074, 0.173205047846,
    0.300000220537, 0.244949027896, 0.264575153589,
    0.244949027896, 0.2999997437, 0.316227942705,
    0.141421228647, 0.458257555962, 0.223606646061,
    0.616441547871, 0.616441547871, 0.700000047684,
    0.264575153589, 0.331662654877, 0.424264162779,
    0.173204988241, 0.223606646061, 0.331662625074,
    0.509901940823, 0.435889661312, 0.2999997437,
    0.173205047846, 0.31622800231, 0.316227942705
  ),
  index1 = FALSE
)

set.seed(42)
res_spnn <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "rand", verbose = FALSE, n_threads = 1,
  nn_method = sparse_nbr_matrix
)
expect_ok_matrix(res_spnn)
# should get same results as with internal nn calculation
expect_equal(res_spnn, res$embedding)

# null X is ok with sparse nearest neighbors
set.seed(42)
res_spnn_nullX <- tumap(
  X = NULL,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "rand", verbose = FALSE, n_threads = 1,
  nn_method = sparse_nbr_matrix0_with_names
)
expect_ok_matrix(res_spnn_nullX)
# output picks up row names from input distance matrix
expect_equal(res_spnn_nullX, res$embedding)

# https://github.com/jlmelville/uwot/issues/6
res <- umap(iris10,
  n_components = 1, n_neighbors = 4, n_epochs = 2,
  n_threads = 1, verbose = FALSE
)
expect_ok_matrix(res, nc = 1)

# enforce irlba for spectral initialization even if RSpectra is present
res <- umap(iris10,
  n_components = 1, n_neighbors = 4, n_epochs = 2,
  n_threads = 1, verbose = FALSE, init = "irlba_spectral"
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

# test that init_sdev actually applies to the input matrix
# store old sd
res_sd <- apply(res, 2, sd)
res2 <- umap(iris10,
  n_neighbors = 4, n_epochs = 0, learning_rate = 0.5,
  init = res, verbose = FALSE, n_threads = 0,
  init_sdev = 5
)
expect_ok_matrix(res2)
expect_equal(apply(res2, 2, sd), rep(5, ncol(res2)))
# make sure input is unchanged
expect_equal(apply(res, 2, sd), res_sd)

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


# 22 Pearson correlation
set.seed(42)
res_cor <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  metric = "correlation", init = "spectral", verbose = FALSE,
  n_threads = 0, ret_model = TRUE
)
expect_ok_matrix(res_cor$embedding)

# Ensure cosine results are different from correlation
set.seed(42)
res_cos <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  metric = "cosine", init = "spectral", verbose = FALSE,
  n_threads = 0, ret_model = TRUE
)
expect_gt(sum((res_cor$embedding - res_cos$embedding)^2), 1e-3)


# Ensure correlation transform results differ from cosine
set.seed(42)
res_trans_cor <- umap_transform(x2m(iris[11:20, ]), res_cor, n_threads = 0, verbose = FALSE)
expect_ok_matrix(res_trans_cor)

# Switch metric and results should differ
res_cor$nn_index$metric <- "cosine"
set.seed(42)
res_trans_cor2 <- umap_transform(x2m(iris[11:20, ]), res_cor, n_threads = 0, verbose = FALSE)
expect_ok_matrix(res_trans_cor2)
expect_gt(sum((res_trans_cor - res_trans_cor2)^2), 1e-3)


# 81: Preserve row names
set.seed(42)
xnames <-
  data.frame(matrix(rnorm(10 * 4), nrow = 10), row.names = letters[1:10])
xumap <-
  umap(
    xnames,
    n_neighbors = 4,
    verbose = FALSE,
    n_threads = 0,
    ret_model = TRUE,
    ret_nn = TRUE
  )
expect_equal(row.names(xumap$embedding), row.names(xnames))
expect_equal(row.names(xumap$nn$euclidean$idx), row.names(xnames))
expect_equal(row.names(xumap$nn$euclidean$dist), row.names(xnames))

first_coords <- c()
test_callback <- function(epochs, n_epochs, coords) {
  first_coords <<- c(first_coords, coords[1, 1])
}
set.seed(42)
ibatch <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, batch = TRUE,
  n_threads = 0, n_sgd_threads = 0, ret_model = TRUE,
  epoch_callback = test_callback
)
expect_equal(length(first_coords), 2)

set.seed(42)
ibatch2 <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, batch = TRUE,
  n_threads = 0, n_sgd_threads = 2, ret_model = TRUE
)
expect_equal(ibatch$embedding, ibatch2$embedding)

itest <- x2m(iris[11:20, ])
first_coords <- c()
fixed_first_coords <- c()
test_transform_callback <- function(epochs, n_epochs, coords, fixed_coords) {
  first_coords <<- c(first_coords, coords[1, 1])
  fixed_first_coords <<- c(fixed_first_coords, fixed_coords[1, 1])
}
set.seed(42)
ibatchtest <- umap_transform(itest, ibatch, epoch_callback = test_transform_callback, n_epochs = 5)
expect_equal(length(first_coords), 5)
expect_equal(length(fixed_first_coords), 5)
# coords don't actually change on the first epoch
expect_equal(length(unique(first_coords)), 4)
# if coords are fixed they should be the same at each epoch
expect_equal(length(unique(fixed_first_coords)), 1)

set.seed(42)
ibatchtest2 <- umap_transform(itest, ibatch, n_sgd_threads = 2, n_epochs = 5)
expect_equal(ibatchtest, ibatchtest2)


oargs_umap <- tumap(iris10,
  n_neighbors = 4, n_epochs = 0, learning_rate = 0.5,
  init = "spca", verbose = FALSE, batch = TRUE,
  n_threads = 0, n_sgd_threads = 0, ret_model = TRUE,
  opt_args = list(alpha = 0.4, beta1 = 0.1, beta2 = 0.2, eps = 1e-3)
)
expect_equal(length(oargs_umap$opt_args), 5)
expect_equal(oargs_umap$opt_args$method, "adam")
expect_equal(oargs_umap$opt_args$alpha, 0.4)
expect_equal(oargs_umap$opt_args$beta1, 0.1)
expect_equal(oargs_umap$opt_args$beta2, 0.2)
expect_equal(oargs_umap$opt_args$eps, 1e-3)


oargs_umap <- tumap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "spca", verbose = FALSE, batch = TRUE,
  n_threads = 0, n_sgd_threads = 0, ret_model = TRUE,
  opt_args = list(method = "sgd", alpha = 0.4)
)
expect_equal(length(oargs_umap$opt_args), 2)
expect_equal(oargs_umap$opt_args$method, "sgd")
expect_equal(oargs_umap$opt_args$alpha, 0.4)

# Return sigma and rho
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "spca", verbose = FALSE, n_threads = 0,
  ret_extra = c("sigma")
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)

expected_sigma <- c(
  0.1799, 0.2049, 0.04938, 0.0906, 0.2494, 0.003906, 0.1537, 0.1355, 0.2454, 0.2063
)
sigma <- res$sigma
expect_equal(sigma, expected_sigma, tolerance = 1e-4)

expected_rho <- c(
  0.1414, 0.1732, 0.2449, 0.2449, 0.1414, 0.6164, 0.2646, 0.1732, 0.3, 0.1732
)
rho <- res$rho
expect_equal(rho, expected_rho, 1e-4)

res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "normlaplacian", verbose = FALSE, n_threads = 0, dens_scale = 1
)
expect_ok_matrix(res)
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "normlaplacian", verbose = FALSE, n_threads = 0, dens_scale = 1,
  ret_extra = c("sigma", "localr")
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)
sigma <- res$sigma
expect_equal(sigma, expected_sigma, tolerance = 1e-4)
rho <- res$rho
expect_equal(rho, expected_rho, tolerance = 1e-4)

expected_localr <- c(
  0.3214, 0.3781, 0.2943, 0.3356, 0.3908, 0.6203, 0.4182, 0.3087, 0.5454, 0.3795
)
localr <- res$localr
expect_equal(localr, expected_localr, tolerance = 1e-4)

res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "normlaplacian", verbose = FALSE, n_threads = 0, dens_scale = 1,
  ret_model = TRUE
)
expect_is(res, "list")
expect_ok_matrix(res$embedding)
expected_ai <- c(
  8.072, 2.957, 13.89, 6.181, 2.41, 0.1389, 1.585, 10.34, 0.3076, 2.888
)
ai <- res$ai
expect_equal(ai, expected_ai, tolerance = 1e-4)
expect_equal(res$dens_scale, 1.0)
expect_equal(res$method, "leopold")

res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "normlaplacian", verbose = FALSE, n_threads = 0, dens_scale = 0.5,
  ret_model = TRUE
)
expected_ai05 <- c(
  3.348, 2.027, 4.392, 2.93, 1.83, 0.4392, 1.484, 3.79, 0.6536, 2.003
)
expect_equal(res$ai, expected_ai05, tolerance = 1e-3)
expect_equal(res$dens_scale, 0.5)

ret_trans <- umap_transform(iris10, res)
expect_ok_matrix(res$embedding)

# 97: should be able to create a model without pre-computed nns and allow
# umap_transform to work with pre-computed nns
train_nn <- annoy_nn(
  X = iris10, k = 4, metric = "euclidean", n_threads = 0,
  ret_index = TRUE
)
set.seed(42)
umap_train_x_null <- umap(
  X = NULL, nn_method = train_nn, ret_model = TRUE,
  n_neighbors = 4
)
set.seed(42)
umap_train_x <- umap(X = iris10, ret_model = TRUE, n_neighbors = 4)

# make the test set a different size to the training set
iris9test <- x2m(iris[11:19, ])
query_ref_nn <- annoy_search(
  X = iris9test, k = 4, ann = train_nn$index,
  n_threads = 0
)
row.names(query_ref_nn$dist) <- row.names(iris9test)

# Success
set.seed(42)
umap_test_1 <- umap_transform(
  X = NULL, model = umap_train_x_null,
  nn_method = query_ref_nn
)

# This was throwing an error because umap_train_x doesn't have pre-computed
# neighbors (and there is no reason to insist that it should have just because
# the test data uses them)
set.seed(42)
umap_test_2 <- umap_transform(
  X = NULL, model = umap_train_x,
  nn_method = query_ref_nn
)
expect_equal(umap_test_1, umap_test_2)



res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5,
  init = "laplacian", verbose = FALSE, n_threads = 0,
  init_sdev = 0.1
)
expect_ok_matrix(res)

# 99 init_sdev = "range" range scales input data columns 0-10
res <- umap(iris10, n_neighbors = 4, init_sdev = "range", n_epochs = 0)
expect_equal(apply(res, 2, range), matrix(c(0, 10, 0, 10), ncol = 2))

# init_sdev = "range" should rescale with user-supplied input too
res <-
  umap(
    iris10,
    n_neighbors = 4,
    init = res,
    init_sdev = "range",
    n_epochs = 0
  )
expect_equal(apply(res, 2, range), matrix(c(0, 10, 0, 10), ncol = 2))

# 101 intersect and union
test_that("intersect and union", {
  # expected values are confirmed via the python implementation
  iris10_12 <-
    as(Matrix::drop0(matrix(
      c(
        0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
        1.0000000e+00, 1.0000000e+00, 0.0000000e+00, 1.0000000e+00,
        0.0000000e+00, 3.1662715e-09,
        0.0000000e+00, 0.0000000e+00, 5.2903235e-01, 4.7097012e-01,
        0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
        4.1861886e-01, 1.0000000e+00,
        0.0000000e+00, 5.2903235e-01, 0.0000000e+00, 1.0000000e+00,
        0.0000000e+00, 0.0000000e+00, 1.0000000e+00, 0.0000000e+00,
        5.8137214e-01, 7.9134369e-01,
        0.0000000e+00, 4.7097012e-01, 1.0000000e+00, 0.0000000e+00,
        0.0000000e+00, 0.0000000e+00, 8.1357503e-01, 0.0000000e+00,
        1.0000000e+00, 4.1732001e-01,
        1.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
        0.0000000e+00, 1.0000000e+00, 2.3949760e-01, 9.2372406e-01,
        0.0000000e+00, 0.0000000e+00,
        1.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
        1.0000000e+00, 0.0000000e+00, 0.0000000e+00, 1.5851905e-08,
        0.0000000e+00, 0.0000000e+00,
        0.0000000e+00, 0.0000000e+00, 1.0000000e+00, 8.1357503e-01,
        2.3949760e-01, 0.0000000e+00, 0.0000000e+00, 3.5861197e-01,
        0.0000000e+00, 0.0000000e+00,
        1.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
        9.2372406e-01, 1.5851905e-08, 3.5861197e-01, 0.0000000e+00,
        0.0000000e+00, 3.1848407e-01,
        0.0000000e+00, 4.1861886e-01, 5.8137214e-01, 1.0000000e+00,
        0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
        0.0000000e+00, 0.0000000e+00,
        3.1662715e-09, 1.0000000e+00, 7.9134369e-01, 4.1732001e-01,
        0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 3.1848407e-01,
        0.0000000e+00, 0.0000000e+00
      ),
      byrow = TRUE, nrow = 10
    )), "generalMatrix")

  iris10_34 <-
    as(Matrix::drop0(matrix(
      c(
        0.000000e+00, 1.000000e+00, 1.000000e+00, 9.995531e-01, 1.000000e+00,
        0.000000e+00, 1.000000e+00, 9.995531e-01, 1.000000e+00, 6.160182e-10,
        # 2
        1.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00, 1.000000e+00,
        0.000000e+00, 1.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00,
        # 3
        1.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00,
        0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
        # 4
        9.995531e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
        1.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00, 1.000000e+00,
        # 5
        1.000000e+00, 1.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00,
        0.000000e+00, 1.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00,
        # 6
        0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00,
        0.000000e+00, 3.771700e-08, 1.000000e+00, 0.000000e+00, 0.000000e+00,
        # 7
        1.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00,
        3.771700e-08, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
        # 8
        9.995531e-01, 0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00,
        1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00,
        # 9
        1.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00,
        0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
        # 10
        6.160182e-10, 0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00,
        0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00
      ),
      byrow = TRUE, nrow = 10
    )), "generalMatrix")

  expected_intersect <- matrix(
    c(
      0.0, 0.55551356, 0.6280843, 0.5780957, 1.0, 0.9556667, 0.65937775, 1.0, 0.8031627, 0.7241592,
      0.55551356, 0.0, 1.0, 0.6247057, 0.62943524, 0.0, 0.7119533, 0.0, 1.0, 0.78192693,
      0.6280843, 1.0, 0.0, 0.701477, 0.68993694, 0.0, 0.7589823, 0.0, 0.85413647, 0.8139497,
      0.5780957, 0.6247057, 0.701477, 0.0, 0.0, 0.96441525, 0.7222059, 0.7034965, 0.84200585, 1.0,
      1.0, 0.62943524, 0.68993694, 0.0, 0.0, 0.96303976, 0.99999994, 0.69032043, 0.8358984, 0.0,
      0.9556667, 0.0, 0.0, 0.96441525, 0.96303976, 0.0, 0.345057, 1.0, 0.0, 0.0,
      0.65937775, 0.7119533, 0.7589823, 0.7222059, 0.99999994, 0.345057, 0.0, 0.7411224, 0.0, 0.0,
      1.0, 0.0, 0.0, 0.7034965, 0.69032043, 1.0, 0.7411224, 0.0, 0.0, 0.999708,
      0.8031627, 1.0, 0.85413647, 0.84200585, 0.8358984, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.7241592, 0.78192693, 0.8139497, 1.0, 0.0, 0.0, 0.0, 0.999708, 0.0, 0.0
    ),
    byrow = TRUE, nrow = 10
  )

  expect_equal(as.matrix(simplicial_set_intersect(iris10_12, iris10_34)),
    expected_intersect,
    tol = 1e-7,
    check.attributes = FALSE
  )

  expected_union <- matrix(
    c(
      0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0002958188,
      1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0,
      1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.90659714,
      0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0,
      1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0,
      1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
      1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
      1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0002958188, 1.0, 0.90659714, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0
    ),
    byrow = TRUE, nrow = 10
  )

  expect_equal(as.matrix(simplicial_set_union(iris10_12, iris10_34)),
    expected_union,
    tol = 1e-7,
    check.attributes = FALSE
  )

  expected_intersect_weighted <- matrix(
    c(
      0.0, 0.40025446, 0.39215788, 0.36378226, 1.0, 0.99999994, 0.4569211, 1.0, 0.6892773, 0.51736516,
      0.40025446, 0.0, 1.0, 0.8227322, 0.4237253, 0.0, 0.4646226, 0.0, 1.0, 0.91525006,
      0.39215788, 1.0, 0.0, 0.8513419, 0.41594562, 0.0, 0.89491767, 0.0, 0.9377875, 0.906834,
      0.36378226, 0.8227322, 0.8513419, 0.0, 0.0, 0.7309169, 0.8804487, 0.4460863, 0.94368106, 1.0,
      1.0, 0.4237253, 0.41594562, 0.0, 0.0, 1.0, 1.0, 0.88101995, 0.70143735, 0.0,
      0.99999994, 0.0, 0.0, 0.7309169, 1.0, 0.0, 0.6468519, 0.7858214, 0.0, 0.0,
      0.4569211, 0.4646226, 0.89491767, 0.8804487, 1.0, 0.6468519, 0.0, 0.8816428, 0.0, 0.0,
      1.0, 0.0, 0.0, 0.4460863, 0.88101995, 0.7858214, 0.8816428, 0.0, 0.0, 0.99908066,
      0.6892773, 1.0, 0.9377875, 0.94368106, 0.70143735, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.51736516, 0.91525006, 0.906834, 1.0, 0.0, 0.0, 0.0, 0.99908066, 0.0, 0.0
    ),
    byrow = TRUE, nrow = 10
  )

  expect_equal(as.matrix(simplicial_set_intersect(iris10_12, iris10_34, weight = 0.25)),
    expected_intersect_weighted,
    tol = 1e-7,
    check.attributes = FALSE
  )
})

test_that("can set seed internally", {
  set.seed(42)
  res <-
    umap(
      iris10,
      n_neighbors = 4,
      n_epochs = 10,
      learning_rate = 0.5,
      n_sgd_threads = 1
    )
  expect_ok_matrix(res)

  # default doesn't reset the seed
  res2 <-
    umap(
      iris10,
      n_neighbors = 4,
      n_epochs = 10,
      learning_rate = 0.5,
      n_sgd_threads = 1
    )
  diff12 <- res - res2
  expect_gt(sqrt(sum(diff12 * diff12) / length(diff12)), 0.01)

  # set seed internally same as calling set.seed
  res3 <-
    umap(
      iris10,
      n_neighbors = 4,
      n_epochs = 10,
      learning_rate = 0.5,
      n_sgd_threads = 1,
      seed = 42
    )
  expect_equal(res, res3)

  # creating a model stores the seed but also forces annoy for nearest neighbors
  # which changes the RNG state more than when FNN can be used internally
  # 115: eh actually this is probably due more to irlba than RSpectra?
  res_model <-
    umap(
      iris10,
      n_neighbors = 4,
      n_epochs = 10,
      learning_rate = 0.5,
      n_sgd_threads = 1,
      seed = 42,
      ret_model = TRUE
    )
  expect_equal(res_model$seed, 42)
  diff1m <- res - res_model$embedding
  expect_gt(sqrt(sum(diff1m * diff1m) / length(diff1m)), 1e-6)

  # explicitly set annoy nn and things are reproducible again
  res4 <-
    umap(
      iris10,
      n_neighbors = 4,
      n_epochs = 10,
      learning_rate = 0.5,
      n_sgd_threads = 1,
      seed = 42,
      nn_method = "annoy"
    )
  expect_equal(res_model$embedding, res4)
})

test_that("can provide nn_args", {
  res <-
    umap(
      iris10,
      n_neighbors = 4,
      n_epochs = 10,
      learning_rate = 0.5,
      n_trees = 5,
      nn_args = list(n_trees = 10),
      ret_model = TRUE
    )
  expect_ok_matrix(res$embedding)
  expect_equal(res$nn_args$n_trees, 10)
})


test_that("deterministic negative sampling is reproducible", {
  res_seed42 <- umap(
    iris10,
    seed = 42,
    n_neighbors = 4,
    n_epochs = 2,
    learning_rate = 0.5,
    init = iris10_pca,
    verbose = FALSE,
    n_threads = 0,
    rng_type = "deterministic"
  )
  res_seed1337 <- umap(
    iris10,
    seed = 1337,
    n_neighbors = 4,
    n_epochs = 2,
    learning_rate = 0.5,
    init = iris10_pca,
    verbose = FALSE,
    n_threads = 0,
    rng_type = "deterministic"
  )
  expect_ok_matrix(res_seed42)
  expect_equal(res_seed42, res_seed1337)

  res_seed42_t2 <- umap(
    iris10,
    seed = 42,
    n_neighbors = 4,
    n_epochs = 2,
    learning_rate = 0.5,
    init = iris10_pca,
    verbose = FALSE,
    n_threads = 2,
    rng_type = "deterministic"
  )
  expect_equal(res_seed42, res_seed42_t2)
})
