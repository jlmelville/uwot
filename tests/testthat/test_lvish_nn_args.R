test_that("lvish uses the Annoy tree count supplied in nn_args", {
  messages <- testthat::capture_messages(
    invisible(lvish(
      iris10,
      n_neighbors = 4,
      perplexity = 2,
      nn_method = "annoy",
      nn_args = list(n_trees = 3),
      n_threads = 0,
      n_epochs = 0,
      verbose = TRUE
    ))
  )

  expect_true(any(grepl("n_trees = 3$", trimws(messages))))
})
