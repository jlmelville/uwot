# Small -ve distances are possible
dist2 <- function(X) {
  D2 <- rowSums(X * X)
  D2 + sweep(X %*% t(X) * -2, 2, t(D2), `+`)
}

# Squared Euclidean distances, ensuring no small -ve distances can occur
safe_dist2 <- function(X) {
  D2 <- dist2(X)
  D2[D2 < 0] <- 0
  D2
}

# convert dataframe to distance matrix
x2d <- function(X) {
  sqrt(safe_dist2(x2m(X)))
}

# Covert a vector into a 2D matrix for generating Y output
c2y <- function(...) {
  matrix(unlist(list(...)), ncol = 2)
}

iris10 <- NULL
iris10_Y <- NULL
diris10 <- NULL
dmiris10 <- NULL
dmiris10z <- NULL
ycat <- NULL
ycat2 <- NULL
ynum <- NULL
ynum2 <- NULL
nn <- NULL

ui10 <- NULL
unn4 <- NULL
self_unn4 <- NULL

create_data <- function() {
  iris10 <<- x2m(iris[1:10, ])
  iris10_Y <<- pca_init(iris10, ndim = 2)
  diris10 <<- stats::dist(iris10)

  # Sparse iris10 dist
  dmiris10 <<- as.matrix(diris10)

  dmiris10zl <- dmiris10
  dmiris10zl[dmiris10zl > 0.71] <- 0
  dmiris10z <<- as(Matrix::drop0(dmiris10zl), "generalMatrix")

  # some Y data
  ycat <<- as.factor(c(levels(iris$Species)[rep(1:3, each = 3)], NA))
  ycat2 <<- as.factor(c(NA, levels(iris$Species)[rep(1:3, times = 3)]))
  ynum <<- (1:10) / 10
  ynum2 <<- seq(from = 10, to = -10, length.out = 10) / 100

  nnl <- find_nn(iris10,
    k = 4, method = "fnn", metric = "euclidean",
    n_threads = 0, verbose = FALSE
  )
  row.names(nnl$idx) <- row.names(iris10)
  row.names(nnl$dist) <- row.names(iris10)
  nn <<- nnl


  # ten iris entries where the 4 nearest neighbors are distinct
  uiris <- unique(iris)
  uirism <- as.matrix(uiris[, -5])
  ui10 <<- uirism[6:15, ]
  unn4 <<- list(
    idx = matrix(c(
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
    ), nrow = 10, byrow = TRUE),
    dist = matrix(c(
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
  )

  self_unn4 <<- list(
    idx = matrix(c(
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
    ), nrow = 10, byrow = TRUE),
    dist = matrix(c(
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
  )
}

# Just test that res is a matrix with valid numbers
expect_ok_matrix <- function(res, nr = 10, nc = 2) {
  expect_is(res, "matrix")
  expect_equal(nrow(res), nr)
  expect_equal(ncol(res), nc)
  expect_false(any(is.infinite(res)))
}

expect_is_nn <- function(res, nr = 10, k = 4) {
  expect_is(res, "list")
  expect_is_nn_matrix(res$dist, nr, k)
  expect_is_nn_matrix(res$idx, nr, k)
}

expect_is_nn_matrix <- function(res, nr = 10, k = 4) {
  expect_is(res, "matrix")
  expect_equal(nrow(res), nr)
  expect_equal(ncol(res), k)
}

create_data()
