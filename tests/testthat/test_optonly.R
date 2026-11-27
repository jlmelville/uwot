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
