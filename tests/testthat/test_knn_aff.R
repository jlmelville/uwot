library(uwot)
context("knn affinity")

expected_sparse <- matrix(0, nrow = 10, ncol = 10)
for (i in seq_len(nrow(nn$idx))) {
  for (j in seq_len(ncol(nn$idx))) {
    expected_sparse[i, nn$idx[i, j]] <- 2
  }
}
expected_sparse <- Matrix::drop0(expected_sparse)

res <- nn_to_sparse(nn$idx, val = 2)
expect_equal(res, expected_sparse)

v <- 1
expected_sparse_mv <- matrix(0, nrow = 10, ncol = 10)
for (i in seq_len(nrow(nn$idx))) {
  nnr <- sort(nn$idx[i, ])
  for (j in seq_len(ncol(nn$idx))) {
    expected_sparse_mv[i, nnr[j]] <- v
    v <- v + 1
  }
}
expect_equal(nn_to_sparse(nn$idx, matrix(1:40, nrow = 10, byrow = TRUE)),
  Matrix::drop0(expected_sparse_mv),
  check.attributes = FALSE
)

res <- perplexity_similarities(iris10, 4, kernel = "knn", nn = nn)
expected_sym_nn_graph <- matrix(0, nrow = 10, ncol = 10)
o3 <- 1 / 3
o6 <- 1 / 6
expected_sym_nn_graph[1, c(5, 6, 8, 10)] <- c(o3, o6, o3, o6)
expected_sym_nn_graph[2, c(3, 4, 9, 10)] <- c(o3, o6, o6, o3)
expected_sym_nn_graph[3, c(2, 4, 7, 9, 10)] <- c(o3, o3, o3, o6, o6)
expected_sym_nn_graph[4, c(2, 3, 7, 9, 10)] <- c(o6, o3, o6, o3, o3)
expected_sym_nn_graph[5, c(1, 6, 7, 8)] <- c(o3, o6, o6, o3)
expected_sym_nn_graph[6, c(1, 5, 8)] <- c(o6, o6, o6)
expected_sym_nn_graph[7, c(3, 4, 5, 8)] <- c(o3, o6, o6, o6)
expected_sym_nn_graph[8, c(1, 5, 6, 7, 10)] <- c(o3, o3, o6, o6, o6)
expected_sym_nn_graph[9, c(2, 3, 4)] <- c(o6, o6, o3)
expected_sym_nn_graph[10, c(1, 2, 3, 4, 8)] <- c(o6, o3, o6, o3, o6)

expect_equal(sum(res), 10)
expect_true(Matrix::isSymmetric(res))
expect_equal(as.matrix(res), expected_sym_nn_graph,
  check.attributes = FALSE,
  tol = 1e-7
)
