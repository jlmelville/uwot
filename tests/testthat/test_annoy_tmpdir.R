expect_annoy_tmpdir <- function(fit) {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  messages <- testthat::capture_messages(
    invisible(fit(
      iris10,
      nn_method = "annoy",
      n_neighbors = 4,
      n_threads = 2,
      tmpdir = dir,
      verbose = TRUE
    ))
  )
  expect_true(any(grepl(file.path(dir, ""), messages, fixed = TRUE)))
}

test_that("umap writes its Annoy index in the requested directory", {
  expect_annoy_tmpdir(umap)
})

test_that("umap2 writes its Annoy index in the requested directory", {
  expect_annoy_tmpdir(umap2)
})

test_that("similarity_graph writes its Annoy index in the requested directory", {
  expect_annoy_tmpdir(similarity_graph)
})
