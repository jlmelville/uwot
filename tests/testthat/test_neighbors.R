library(uwot)
context("neighbors")

i10nn4dist <-  matrix(c(
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

i10nn4idx <-  matrix(c(
  1,    5,    8,   10,
  2,   10,    3,    4,
  3,    4,    7,    2,
  4,    3,    9,   10,
  5,    1,    8,    7,
  6,    1,    5,    8,
  7,    3,    4,    8,
  8,    1,    5,   10,
  9,    4,    3,    2,
  10,   2,    3,    4
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

