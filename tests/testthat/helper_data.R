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

iris10 <- x2m(iris[1:10, ])
iris10_Y <- pca_scores(iris10, ncol = 2)

diris10 <- dist(iris10)

# Sparse iris10 dist
dmiris10 <- as.matrix(diris10)
dmiris10z <- dmiris10
dmiris10z[dmiris10z > 0.71] <- 0
dmiris10z <- Matrix::drop0(dmiris10z)

# some Y data
ycat <- as.factor(c(levels(iris$Species)[rep(1:3, each = 3)], NA))
ycat2 <- as.factor(c(NA, levels(iris$Species)[rep(1:3, times = 3)]))
ynum <- (1:10) / 10
ynum2 <- seq(from = 10, to = -10, length.out = 10) / 100


nn <- find_nn(iris10,
  k = 4, method = "fnn", metric = "euclidean",
  n_threads = 0, verbose = FALSE
)

# Just test that res is a matrix with valid numbers
expect_ok_matrix <- function(res, nr = nrow(iris10), nc = 2) {
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
