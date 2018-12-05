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

# Results are repeatable with n_threads = 0 (or 1) and same seed
set.seed(1337)
res2 <- umap(iris10,
            n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
            init = "normlaplacian", verbose = FALSE, n_threads = 0
)
expect_equal(res2, res)

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

# init with matrix
iris10_pca <- prcomp(iris10, rank. = 2, retx = TRUE, center = TRUE,
                     scale. = FALSE)$x

res <- umap(iris10,
            n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
            init = iris10_pca, verbose = FALSE, n_threads = 0
)
expect_ok_matrix(res)
# Ensure that internal C++ code doesn't modify user-supplied initialization
expect_equal(iris10_pca, prcomp(iris10, rank. = 2, retx = TRUE, center = TRUE,
                            scale. = FALSE)$x)

# return nn
# reset seed here so we can compare output with next test result
set.seed(1337)
res <- umap(iris10,
            n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
            init = "spca", verbose = FALSE, n_threads = 0,
            ret_nn = TRUE)
expect_is(res, "list")
expect_ok_matrix(res$embedding)
expect_is(res$nn, "list")
expect_is(res$nn$euclidean, "list")
expect_ok_matrix(res$nn$euclidean$idx, nc = 4)
expect_ok_matrix(res$nn$euclidean$dist, nc = 4)

# Use pre-calculated nn: should be the same as previous result
set.seed(1337)
res_nn <- umap(iris10,
            nn_method = res$nn[[1]], n_epochs = 2, alpha = 0.5, min_dist = 0.001,
            init = "spca", verbose = FALSE, n_threads = 0)
expect_ok_matrix(res_nn)
expect_equal(res_nn, res$embedding)

# Passing nn list directly is also ok
set.seed(1337)
res_nnl <- umap(iris10,
               nn_method = res$nn, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
               init = "spca", verbose = FALSE, n_threads = 0,
               ret_nn = TRUE)
expect_ok_matrix(res_nnl$embedding)
expect_equal(res_nnl$embedding, res$embedding)
expect_equal(res_nnl$nn[[1]], res$nn[[1]])
expect_equal(names(res_nnl$nn), "precomputed")

# Use multiple nn data
res_nn2 <- umap(iris10,
               nn_method = list(nn, nn), n_epochs = 2, alpha = 0.5, 
               min_dist = 0.001,
               init = "spca", verbose = FALSE, n_threads = 0, ret_nn = TRUE)
expect_ok_matrix(res_nn2$embedding)
expect_equal(names(res_nn2$nn), c("precomputed", "precomputed"))


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

# return nn and a model
res <- umap(iris10,
            n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
            init = "rand", verbose = FALSE, n_threads = 1,
            ret_model = TRUE, ret_nn = TRUE)
expect_is(res, "list")
expect_ok_matrix(res$embedding)

expect_is(res$nn, "list")
expect_is_nn(res$nn[[1]], k = 4)
expect_equal(names(res$nn), "euclidean")
res_test <- umap_transform(iris10, res, n_threads = 0, verbose = FALSE)
expect_ok_matrix(res_test)


# https://github.com/jlmelville/uwot/issues/6
res <- umap(iris10, n_components = 1, n_neighbors = 4, n_epochs = 2,
            n_threads = 1, verbose = FALSE)
expect_ok_matrix(res, nc = 1)


# Supervised
set.seed(1337)
res_y <- umap(iris10, n_neighbors = 4, n_epochs = 2, alpha = 0.5,
               min_dist = 0.001, init = "spca", verbose = FALSE, n_threads = 0,
               y = 1 / (1:10) ^ 2, target_n_neighbors = 2)
expect_ok_matrix(res_y)

# Repeat using equivalent NN info for y
y_nn <- list(
  idx = matrix(c(
    1,    2,
    2,    3,
    3,    4,
    4,    5,
    5,    6,
    6,    7,
    7,    8,
    8,    9,
    9,   10,
    10,    9
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
res_ynn <- umap(iris10, n_neighbors = 4, n_epochs = 2, alpha = 0.5,
               min_dist = 0.001, init = "spca", verbose = FALSE, n_threads = 0,
               y = y_nn)
expect_ok_matrix(res_ynn)
# Should be the same result
expect_equal(res_ynn, res_y)

hamm <- structure(c(0L, 0L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
                    0L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L,
                    0L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L), .Dim = c(10L, 4L
                    ))
res <- umap(hamm, n_neighbors = 4, metric = "hamming", verbose = FALSE,
            n_threads = 1)
expect_ok_matrix(res)

# Multiple metrics
set.seed(1337)
res <- umap(iris10,
            n_neighbors = 4, n_epochs = 2, alpha = 0.5,
            init = "spca", verbose = FALSE, n_threads = 0,
            metric = list(euclidean = c(1, 2), euclidean = c(3, 4)),
            ret_model = TRUE)
res_trans <- umap_transform(iris10, model = res, verbose = FALSE, n_threads = 0,
                            n_epochs = 2)
expect_ok_matrix(res_trans)

