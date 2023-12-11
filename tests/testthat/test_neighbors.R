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


sparse_to_tri <- function(m, lower = TRUE) {
  sm <- summary(m)
  if (lower) {
    subtri <- subset(sm, i >= j)
  } else {
    subtri <- subset(sm, i <= j)
  }

  Matrix::sparseMatrix(i = subtri$i, j = subtri$j, x = subtri$x, dims = dim(m))
}

dmiris10zu <- sparse_to_tri(dmiris10z, lower = FALSE)
res <- sparse_tri_nn(dmiris10zu, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

dmiris10zl <- sparse_to_tri(dmiris10z, lower = TRUE)
res <- sparse_tri_nn(dmiris10zl, k = 4, include_self = TRUE)
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

res <- find_nn(dmiris10z, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

res <- find_nn(dmiris10zu, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])

res <- find_nn(dmiris10zl, k = 4, include_self = TRUE)
expect_equal(res$dist, i10nn4dist, tol = 1e-6)
expect_equal(res$idx[-6, ], i10nn4idx[-6, ])


# Test Annoy
res <- annoy_nn(ui10, k = 4, n_threads = 0)
expect_equal(res$idx, self_unn4$idx, check.attributes = FALSE)
expect_equal(res$dist, self_unn4$dist, check.attributes = FALSE, tol = 1e-6)

res <- annoy_nn(ui10, k = 4, n_threads = 0, ret_index = TRUE)
expect_equal(res$idx, self_unn4$idx, check.attributes = FALSE)
expect_equal(res$dist, self_unn4$dist, check.attributes = FALSE, tol = 1e-6)
expect_true(!is.null(res$index))
expect_is(res$index, "list")
expect_is(res$index$ann, "Rcpp_AnnoyEuclidean")
expect_equal(res$index$metric, "euclidean")


res <- annoy_nn(ui10, k = 4, n_threads = 1)
expect_equal(res$idx, self_unn4$idx, check.attributes = FALSE)
expect_equal(res$dist, self_unn4$dist, check.attributes = FALSE, tol = 1e-6)

res <- annoy_nn(ui10, k = 4, n_threads = 1, ret_index = TRUE)
expect_equal(res$idx, self_unn4$idx, check.attributes = FALSE)
expect_equal(res$dist, self_unn4$dist, check.attributes = FALSE, tol = 1e-6)
expect_true(!is.null(res$index))
expect_is(res$index, "list")
expect_is(res$index$ann, "Rcpp_AnnoyEuclidean")
expect_equal(res$index$metric, "euclidean")

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
expect_is(res$index, "list")
expect_is(res$index$ann, "Rcpp_AnnoyAngular")
expect_equal(res$index$metric, "cosine")



# Correlation distances from sklearn
cor_index <- matrix(
  c(
    1, 3, 6, 8, 5, 7, 4, 9, 10, 2,
    2, 10, 9, 4, 8, 3, 1, 6, 5, 7,
    3, 1, 8, 6, 5, 7, 9, 4, 10, 2,
    4, 9, 8, 10, 3, 1, 6, 2, 5, 7,
    5, 7, 6, 1, 3, 8, 4, 9, 10, 2,
    6, 5, 7, 1, 3, 8, 4, 9, 10, 2,
    7, 5, 6, 1, 3, 8, 4, 9, 10, 2,
    8, 3, 4, 1, 9, 6, 5, 10, 7, 2,
    9, 4, 10, 8, 2, 3, 1, 6, 5, 7,
    10, 9, 4, 2, 8, 3, 1, 6, 5, 7
  ),
  byrow = TRUE, ncol = 10, nrow = 10
)
cor_dist <- matrix(
  c(
    0.00000000e+00, 2.60889537e-05, 4.13946853e-04, 4.61852891e-04,
    6.52668500e-04, 1.18880691e-03, 1.83154822e-03, 1.92335771e-03,
    3.44803882e-03, 4.00133876e-03,
    0.00000000e+00, 9.67144550e-04, 1.45365862e-03, 2.60336976e-03,
    2.88194472e-03, 3.39291433e-03, 4.00133876e-03, 6.40810993e-03,
    7.76732119e-03, 9.27944638e-03,
    0.00000000e+00, 2.60889537e-05, 3.95491414e-04, 6.22703492e-04,
    9.38868047e-04, 1.56232107e-03, 1.64390757e-03, 1.66652079e-03,
    3.01440442e-03, 3.39291433e-03,
    0.00000000e+00, 1.66690505e-04, 4.54440488e-04, 6.93152364e-04,
    1.66652079e-03, 1.83154822e-03, 2.16742062e-03, 2.60336976e-03,
    3.28117923e-03, 3.86063303e-03,
    0.00000000e+00, 8.60443273e-05, 1.16680316e-04, 6.52668500e-04,
    9.38868047e-04, 1.49684289e-03, 3.28117923e-03, 3.96913822e-03,
    6.23882524e-03, 7.76732119e-03,
    0.00000000e+00, 1.16680316e-04, 2.77394147e-04, 4.13946853e-04,
    6.22703492e-04, 8.21174669e-04, 2.16742062e-03, 2.78434617e-03,
    4.73942264e-03, 6.40810993e-03,
    1.11022302e-16, 8.60443273e-05, 2.77394147e-04, 1.18880691e-03,
    1.56232107e-03, 2.04787836e-03, 3.86063303e-03, 4.78602797e-03,
    7.27277830e-03, 9.27944638e-03,
    0.00000000e+00, 3.95491414e-04, 4.54440488e-04, 4.61852891e-04,
    5.93789371e-04, 8.21174669e-04, 1.49684289e-03, 1.62634825e-03,
    2.04787836e-03, 2.88194472e-03,
    0.00000000e+00, 1.66690505e-04, 2.60225275e-04, 5.93789371e-04,
    1.45365862e-03, 1.64390757e-03, 1.92335771e-03, 2.78434617e-03,
    3.96913822e-03, 4.78602797e-03,
    0.00000000e+00, 2.60225275e-04, 6.93152364e-04, 9.67144550e-04,
    1.62634825e-03, 3.01440442e-03, 3.44803882e-03, 4.73942264e-03,
    6.23882524e-03, 7.27277830e-03
  ),
  byrow = TRUE, ncol = 10, nrow = 10
)

res <- annoy_nn(iris10, k = 10, n_threads = 0, ret_index = TRUE, metric = "correlation")
expect_equal(res$idx, cor_index, check.attributes = FALSE)
expect_equal(res$dist, cor_dist, check.attributes = FALSE, tol = 1e-6)
expect_true(!is.null(res$index))
expect_is(res$index, "list")
expect_is(res$index$ann, "Rcpp_AnnoyAngular")
expect_equal(res$index$metric, "correlation")
