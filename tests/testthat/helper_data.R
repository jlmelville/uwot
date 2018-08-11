# convert data frame to matrix using numeric columns
x2m <- function(X) {
  if (!methods::is(X, "matrix")) {
    m <- as.matrix(X[, which(vapply(X, is.numeric, logical(1)))])
    attr(m, "dimnames") <- NULL
  }
  else {
    m <- X
  }
  m
}

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

nn <- find_nn(iris10,
  k = 4, method = "fnn", metric = "euclidean",
  n_threads = 0, verbose = FALSE
)

# Just test that res is a matrix with valid numbers
expect_ok_matrix <- function(res) {
  expect_is(res, "matrix")
  expect_false(any(is.infinite(res)))
}
