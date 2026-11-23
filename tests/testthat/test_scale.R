library(uwot)

test_that("standardization centers and rescales each feature", {
  iris10_sd <- apply(iris10, 2, sd)
  iris10_mean <- apply(iris10, 2, mean)

  iris10_none <- scale_input(iris10, scale_type = FALSE)
  expect_equal(apply(iris10_none, 2, sd), iris10_sd)
  expect_equal(apply(iris10_none, 2, mean), iris10_mean)

  iris10_scale <- scale_input(iris10, scale_type = TRUE)
  expect_equal(unname(apply(iris10_scale, 2, sd)), rep(1, 4))
  expect_equal(
    unname(apply(iris10_scale, 2, mean)),
    rep(0, 4)
  )

  # "scale" and "z" and TRUE are synonyms
  expect_equal(scale_input(iris10, scale_type = "scale"), iris10_scale)
  expect_equal(scale_input(iris10, scale_type = "Z"), iris10_scale)
})

test_that("range and maximum-absolute scaling use the requested bounds", {
  iris10_maxabs <- scale_input(iris10, scale_type = "maxabs")
  expect_equal(
    unname(apply(iris10_maxabs, 2, mean)),
    rep(0, 4)
  )
  expect_equal(max(abs(iris10_maxabs)), 1)

  iris10_range <- scale_input(iris10, scale_type = "range")
  expect_equal(max(iris10_range), 1)
  expect_equal(min(iris10_range), 0)

  iris10_colrange <- scale_input(iris10, scale_type = "colrange")
  expect_equal(
    unname(apply(iris10_colrange, 2, max)),
    rep(1, 4)
  )
  expect_equal(
    unname(apply(iris10_colrange, 2, min)),
    rep(0, 4)
  )
})

test_that("scaling applied outside umap should not appear in the model", {
  iris10s <- scale(iris10)
  iris10s_umap0 <-
    umap(
      iris10s,
      n_neighbors = 4,
      n_epochs = 0,
      init = "rand",
      ret_model = TRUE,
      scale = FALSE
    )
  expect_null(iris10s_umap0$scale_info)

  iris10s_umap0s <-
    umap(
      iris10s,
      n_neighbors = 4,
      n_epochs = 0,
      init = "rand",
      ret_model = TRUE,
      scale = TRUE
    )
  expect_equal(
    names(iris10s_umap0s$scale_info),
    c("scaled:center", "scaled:scale", "scaled:nzvcols")
  )
})
