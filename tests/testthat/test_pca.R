library(uwot)
context("PCA")

iris10prcomp <- prcomp(iris10, retx = TRUE, center = TRUE, scale. = FALSE)

test_that("PCA initialization", {
  iris10_pca_scores <- pca_init(iris10, ndim = 2)
  suppressWarnings(iris10_irlba_scores <- irlba_scores(iris10, ncol = 2))

  expect_equal(abs(iris10prcomp$x[, 1:2]), abs(iris10_pca_scores),
    check.attributes = FALSE
  )
  expect_equal(abs(iris10prcomp$x[, 1:2]), abs(iris10_irlba_scores),
    check.attributes = FALSE
  )

  suppressWarnings(iris10_svdr_scores <- irlba_svdr_scores(iris10, ncol = 2))
  expect_equal(abs(iris10prcomp$x[, 1:2]), abs(iris10_svdr_scores),
    check.attributes = FALSE
  )
})

test_that("1 component initialization works", {
  expect_ok_matrix(pca_init(iris10, ndim = 1), nc = 1)
})

test_that("PCA returns model data", {
  iris10_pca_scores <- pca_init(iris10, ndim = 2, ret_extra = TRUE)
  expect_equal(abs(iris10prcomp$x[, 1:2]),
    abs(iris10_pca_scores$scores),
    check.attributes = FALSE
  )
  expect_equal(abs(iris10prcomp$rotation[, 1:2]),
    abs(iris10_pca_scores$rotation),
    check.attributes = FALSE
  )
  expect_equal(abs(iris10prcomp$center),
    abs(iris10_pca_scores$center),
    check.attributes = FALSE
  )

  suppressWarnings(iris10_irlba_scores <- irlba_scores(iris10,
    ncol = 2,
    ret_extra = TRUE
  ))
  expect_equal(abs(iris10prcomp$x[, 1:2]),
    abs(iris10_irlba_scores$scores),
    check.attributes = FALSE
  )
  expect_equal(abs(iris10prcomp$rotation[, 1:2]),
    abs(iris10_irlba_scores$rotation),
    check.attributes = FALSE
  )
  expect_equal(abs(iris10prcomp$center),
    abs(iris10_irlba_scores$center),
    check.attributes = FALSE
  )
})

test_that("logical pca ok", {
  set.seed(1337)
  random_logical <- matrix(rnorm(1000), nrow = 100) > 0.5
  random_int <- random_logical * 1

  expect_equal(abs(irlba_scores(random_logical, ncol = 2)),
               abs(irlba_scores(random_int, ncol = 2)))
})
