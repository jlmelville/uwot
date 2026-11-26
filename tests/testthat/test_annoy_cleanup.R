# Interrupt after the search index is available and before the native query.
interrupt_annoy_search <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (
      grepl("Searching Annoy index using", conditionMessage(m), fixed = TRUE)
    ) {
      stop("interrupted Annoy search")
    }
    invokeRestart("muffleMessage")
  })
}

make_annoy_model <- function(threads = 0, dir = tempdir()) {
  umap(
    iris10,
    nn_method = "annoy",
    n_neighbors = 4,
    init = "random",
    n_epochs = 3,
    n_threads = threads,
    tmpdir = dir,
    ret_model = TRUE,
    verbose = FALSE
  )
}

test_that("an interrupted Annoy search removes its temporary index", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  expect_error(
    interrupt_annoy_search(tumap(
      iris10,
      nn_method = "annoy",
      n_neighbors = 4,
      n_threads = 2,
      tmpdir = dir,
      verbose = TRUE
    )),
    "interrupted Annoy search"
  )
  expect_length(list.files(dir), 0)
})

test_that("a failed threaded transform preserves its loaded Annoy model", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  archive <- tempfile()
  on.exit(unlink(archive), add = TRUE)
  model <- make_annoy_model()
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
    interrupt_annoy_search(project(2, TRUE)),
    "interrupted Annoy search"
  )
  expect_length(list.files(dir), 0)
  expect_equal(project(0, FALSE), expected)
})

test_that("a live Annoy model reuses its backing file until unload", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  model <- make_annoy_model()
  query <- function(threads, verbose = FALSE) {
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
  expected <- query(0)
  expect_error(
    interrupt_annoy_search(query(2, TRUE)),
    "interrupted Annoy search"
  )
  backing <- list.files(dir)
  expect_length(backing, as.integer(.Platform$OS.type == "windows"))
  expect_equal(query(2), expected)
  expect_identical(list.files(dir), backing)
  unload_uwot(model)
  expect_length(list.files(dir), 0)
})

test_that("collecting an Annoy model releases its temporary backing file", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  old_wd <- setwd(dirname(dir))
  on.exit(setwd(old_wd), add = TRUE)
  model <- make_annoy_model(2, basename(dir))
  setwd(old_wd)
  expect_length(list.files(dir), as.integer(.Platform$OS.type == "windows"))
  alias <- model
  rm(model)
  invisible(gc())
  expect_length(list.files(dir), as.integer(.Platform$OS.type == "windows"))
  rm(alias)
  invisible(gc())
  expect_length(list.files(dir), 0)
})


test_that("saving a live Annoy model releases its old temporary backing", {
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  model <- make_annoy_model(2, dir)
  archive <- file.path(dir, "model.tar")
  saved <- save_uwot(model, archive)
  on.exit(unload_uwot(saved), add = TRUE, after = FALSE)

  expect_identical(list.files(dir), "model.tar")
  expect_true(all(is.finite(umap_transform(iris10, saved, n_threads = 0))))
})

test_that("archive failures preserve the caller's live Annoy index", {
  model <- make_annoy_model()
  on.exit(unload_uwot(model, cleanup = FALSE))
  expected <- query_neighbors(model)
  wd <- getwd()
  original_dirs <- list.dirs(tempdir(), recursive = FALSE)
  missing <- file.path(tempfile(), "model.tar")
  expect_error(suppressWarnings(save_uwot(model, missing, unload = TRUE)))
  expect_false(file.exists(missing))
  expect_equal(getwd(), wd)
  expect_equal(query_neighbors(model), expected)

  # An external archiver can return a failure status without raising an R error.
  old_tar <- Sys.getenv("tar", unset = NA_character_)
  on.exit(
    if (is.na(old_tar)) Sys.unsetenv("tar") else Sys.setenv(tar = old_tar),
    add = TRUE
  )
  stub <- tempfile(fileext = ".R")
  on.exit(unlink(stub), add = TRUE)
  writeLines("quit(status = 1L)", stub)
  Sys.setenv(
    tar = paste(
      shQuote(file.path(R.home("bin"), "Rscript")),
      "--vanilla",
      shQuote(stub)
    )
  )
  archive <- tempfile()
  expect_error(
    suppressWarnings(save_uwot(model, archive, unload = TRUE)),
    "archive"
  )
  expect_false(file.exists(archive))
  expect_equal(query_neighbors(model), expected)
  unload_uwot(model, cleanup = FALSE)
  expect_equal(list.dirs(tempdir(), recursive = FALSE), original_dirs)
})
