expect_annoy_tmpdir <- function(fit) {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  saw_index <- FALSE
  withCallingHandlers(
    fit(
      iris10,
      nn_method = "annoy",
      n_neighbors = 4,
      n_threads = 2,
      tmpdir = dir,
      verbose = TRUE
    ),
    message = function(m) {
      files <- list.files(dir, full.names = TRUE)
      if (any(file.size(files) > 0)) {
        saw_index <<- TRUE
      }
      invokeRestart("muffleMessage")
    }
  )
  expect_true(saw_index)
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
