library(uwot)

test_that("approximate-power optimization updates the precomputed layout", {
  init <- cbind(seq_len(nrow(iris10)), rev(seq_len(nrow(iris10))))
  layout <- function(approx_pow) {
    umap(
      NULL,
      nn_method = nn,
      init = init,
      n_epochs = 10,
      seed = 42,
      n_threads = 0,
      batch = TRUE,
      approx_pow = approx_pow
    )
  }
  approx <- layout(TRUE)
  expect_ok_matrix(approx)
  expect_false(isTRUE(all.equal(
    as.numeric(dist(approx)),
    as.numeric(dist(init))
  )))
  expect_false(isTRUE(all.equal(approx, layout(FALSE))))
})

test_that("batch optimization visits every epoch across worker windows", {
  # More than two polling windows of edges, including with two workers.
  n <- 370L
  graph <- Matrix::Matrix(matrix(1, n, n), sparse = TRUE)
  diag(graph) <- 0
  set.seed(19)
  init <- matrix(rnorm(2 * n), n, 2)
  layout <- function(threads) {
    epochs <- integer()
    result <- optimize_graph_layout(
      graph,
      init = init,
      n_epochs = 3,
      batch = TRUE,
      n_sgd_threads = threads,
      rng_type = "deterministic",
      epoch_callback = function(epoch, n_epochs, coords) {
        epochs <<- c(epochs, epoch)
        expect_ok_matrix(coords, nr = n)
      }
    )
    expect_equal(epochs, seq_len(3))
    result
  }
  one <- layout(1)
  expect_false(isTRUE(all.equal(as.vector(dist(one)), as.vector(dist(init)))))
  expect_equal(layout(2), one)
})

test_that("batch Tausworthe sampling is reproducible with two workers", {
  layout <- function() {
    umap(
      iris10,
      n_neighbors = 4,
      init = "random",
      n_epochs = 4,
      batch = TRUE,
      n_threads = 0,
      n_sgd_threads = 2,
      rng_type = "tausworthe",
      seed = 43
    )
  }
  embedding <- layout()
  expect_ok_matrix(embedding)
  expect_equal(layout(), embedding)
})
