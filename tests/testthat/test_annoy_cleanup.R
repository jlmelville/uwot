test_that("an interrupted Annoy search removes its temporary index", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  expect_error(
    withCallingHandlers(
      tumap(
        iris10,
        nn_method = "annoy",
        n_neighbors = 4,
        n_threads = 2,
        tmpdir = dir,
        verbose = TRUE
      ),
      message = function(m) {
        if (length(list.files(dir)) > 0L) {
          stop("interrupted after writing the index")
        }
        invokeRestart("muffleMessage")
      }
    ),
    "interrupted after writing the index"
  )
  expect_length(list.files(dir), 0)
})

test_that("a failed threaded transform preserves its loaded Annoy model", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  archive <- tempfile()
  on.exit(unlink(archive), add = TRUE)
  model <- umap(
    iris10,
    nn_method = "annoy",
    n_neighbors = 4,
    init = "random",
    n_epochs = 3,
    n_threads = 0,
    ret_model = TRUE,
    verbose = FALSE
  )
  save_uwot(model, archive, unload = TRUE)
  model <- load_uwot(archive)
  on.exit(unload_uwot(model), add = TRUE, after = FALSE)
  project <- function(threads, verbose) {
    set.seed(42)
    umap_transform(
      iris10,
      model,
      n_threads = threads,
      n_sgd_threads = 1,
      tmpdir = dir,
      verbose = verbose
    )
  }
  expected <- project(0, FALSE)
  expect_error(
    withCallingHandlers(project(2, TRUE), message = function(m) {
      if (
        grepl("Searching Annoy index using", conditionMessage(m), fixed = TRUE)
      ) {
        stop("interrupted threaded transform")
      }
      invokeRestart("muffleMessage")
    }),
    "interrupted threaded transform"
  )
  expect_length(list.files(dir), 0)
  expect_equal(project(0, FALSE), expected)
})
