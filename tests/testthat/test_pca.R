library(uwot)
context("PCA")

test_that("PCA initialization", {
  iris10prcomp <- prcomp(iris10, retx = TRUE, center = TRUE, scale. = FALSE)$x
  iris10_pca_scores <- pca_scores(iris10, ncol = 2)
  suppressWarnings(iris10_irlba_scores <- irlba_scores(iris10, ncol = 2))

  expect_equal(abs(iris10prcomp[, 1:2]), abs(iris10_pca_scores),
    check.attributes = FALSE
  )
  expect_equal(abs(iris10prcomp[, 1:2]), abs(iris10_irlba_scores),
    check.attributes = FALSE
  )
})

test_that("1 component initialization works", {
  expect_ok_matrix(pca_init(iris10, ndim = 1))
  expect_ok_matrix(scaled_pca(iris10, ndim = 1))
})
