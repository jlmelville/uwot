library(uwot)
context("neighbors")

i10nn4dist <- matrix(c(
  0, 0.1414214, 0.1732051, 0.4690416,
  0, 0.1732051, 0.3000000, 0.3316625,
  0, 0.2449490, 0.2645751, 0.3000000,
  0, 0.2449490, 0.3000000, 0.3162278,
  0, 0.1414214, 0.2236068, 0.4582576,
  0, 0.6164414, 0.6164414, 0.7000000,
  0, 0.2645751, 0.3316625, 0.4242641,
  0, 0.1732051, 0.2236068, 0.3316625,
  0, 0.3000000, 0.4358899, 0.5099020,
  0, 0.1732051, 0.3162278, 0.3162278
), nrow = 10, byrow = TRUE)

i10nn4idx <- matrix(c(
  1, 5, 8, 10,
  2, 10, 3, 4,
  3, 4, 7, 2,
  4, 3, 9, 10,
  5, 1, 8, 7,
  6, 1, 5, 8,
  7, 3, 4, 8,
  8, 1, 5, 10,
  9, 4, 3, 2,
  10, 2, 3, 4
), nrow = 10, byrow = TRUE)

## Test specialized functions
res <- FNN_nn(iris10, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

res <- dist_nn(diris10, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

res <- sparse_nn(dmiris10z, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

# Test overall function
res <- find_nn(iris10, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

res <- find_nn(diris10, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

res <- find_nn(dmiris10z, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

# Test Annoy
# ten iris entries where the 4 nearest neighbors are distinct
uiris <- unique(iris)
uirism <- as.matrix(uiris[, -5])
ui10 <- uirism[6:15, ]

nn_index4 <- matrix(c(
  6, 10, 3, 7,
  7, 3, 5, 8,
  7, 5, 2, 8,
  9, 8, 2, 5,
  8, 3, 7, 2,
  1, 3, 10, 7,
  3, 2, 5, 8,
  5, 4, 7, 3,
  4, 8, 2, 5,
  6, 1, 3, 7
), nrow = 10, byrow = TRUE)

nn_dist4 <- matrix(c(
  0.3464102, 0.6782330, 0.7000000, 0.8124038,
  0.3000000, 0.4242641, 0.4795832, 0.4898979,
  0.2236068, 0.3316625, 0.4242641, 0.4690416,
  0.3464102, 0.4242641, 0.5477226, 0.5567764,
  0.1732051, 0.3316625, 0.3464102, 0.4795832,
  0.3464102, 0.5000000, 0.5830952, 0.6782330,
  0.2236068, 0.3000000, 0.3464102, 0.4582576,
  0.1732051, 0.4242641, 0.4582576, 0.4690416,
  0.3464102, 0.5830952, 0.6164414, 0.7280110,
  0.5830952, 0.6782330, 1.0440307, 1.2328828
), nrow = 10, byrow = TRUE)

self_nn_index4 <- matrix(c(
  1, 6, 10, 3,
  2, 7, 3, 5,
  3, 7, 5, 2,
  4, 9, 8, 2,
  5, 8, 3, 7,
  6, 1, 3, 10,
  7, 3, 2, 5,
  8, 5, 4, 7,
  9, 4, 8, 2,
  10, 6, 1, 3
), nrow = 10, byrow = TRUE)

self_nn_dist4 <- matrix(c(
  0, 0.3464102, 0.6782330, 0.7000000,
  0, 0.3000000, 0.4242641, 0.4795832,
  0, 0.2236068, 0.3316625, 0.4242641,
  0, 0.3464102, 0.4242641, 0.5477226,
  0, 0.1732051, 0.3316625, 0.3464102,
  0, 0.3464102, 0.5000000, 0.5830952,
  0, 0.2236068, 0.3000000, 0.3464102,
  0, 0.1732051, 0.4242641, 0.4582576,
  0, 0.3464102, 0.5830952, 0.6164414,
  0, 0.5830952, 0.6782330, 1.0440307
), nrow = 10, byrow = TRUE)

res <- annoy_nn(ui10, k = 4, n_threads = 0)
expect_equal(res$idx, self_nn_index4, check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4, check.attributes = FALSE, tol = 1e-6)

res <- annoy_nn(ui10, k = 4, n_threads = 0, ret_index = TRUE)
expect_equal(res$idx, self_nn_index4, check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4, check.attributes = FALSE, tol = 1e-6)
expect_true(!is.null(res$index))
expect_is(res$index, "Rcpp_AnnoyEuclidean")

res <- annoy_nn(ui10, k = 4, n_threads = 1)
expect_equal(res$idx, self_nn_index4, check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4, check.attributes = FALSE, tol = 1e-6)

res <- annoy_nn(ui10, k = 4, n_threads = 1, ret_index = TRUE)
expect_equal(res$idx, self_nn_index4, check.attributes = FALSE)
expect_equal(res$dist, self_nn_dist4, check.attributes = FALSE, tol = 1e-6)
expect_true(!is.null(res$index))
expect_is(res$index, "Rcpp_AnnoyEuclidean")

cos_index <- matrix(
  c(
    1, 2, 7, 3,
    2, 1, 7, 3,
    3, 6, 4, 7,
    4, 3, 5, 7,
    5, 8, 4, 3,
    6, 3, 9, 4,
    7, 3, 1, 4,
    8, 5, 4, 3,
    9, 6, 10, 3,
    10, 9, 6, 3
  ),
  byrow = TRUE, ncol = 4, nrow = 10
)


# Cosine distances from HNSW
cos_dist <- matrix(
  c(
    0, 0.000131368637084961, 0.00048297643661499, 0.000737011432647705,
    0, 0.000131368637084961, 0.000680804252624512, 0.000909507274627686,
    0, 0.000168740749359131, 0.000244021415710449, 0.000422179698944092,
    0, 0.000244021415710449, 0.000383198261260986, 0.000549376010894775,
    0, 7.09891319274902e-05, 0.000383198261260986, 0.000682294368743896,
    0, 0.000168740749359131, 0.000265955924987793, 0.000767052173614502,
    0, 0.000422179698944092, 0.00048297643661499, 0.000549376010894775,
    0, 7.09891319274902e-05, 0.000611364841461182, 0.000812351703643799,
    0, 0.000265955924987793, 0.00078284740447998, 0.000819146633148193,
    0, 0.00078284740447998, 0.00160372257232666, 0.00279802083969116
  ),
  byrow = TRUE, ncol = 4, nrow = 10
)

res <- annoy_nn(ui10, k = 4, n_threads = 0, ret_index = TRUE, metric = "cosine")
expect_equal(res$idx, cos_index, check.attributes = FALSE)
expect_equal(res$dist, cos_dist, check.attributes = FALSE, tol = 1e-6)
expect_true(!is.null(res$index))
expect_is(res$index, "Rcpp_AnnoyAngular")
