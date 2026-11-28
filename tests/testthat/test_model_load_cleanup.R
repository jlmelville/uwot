test_that("an incomplete model archive leaves no extraction directory", {
  archive <- tempfile(fileext = ".tar")
  on.exit(unlink(archive))
  utils::tar(archive, files = character(), tar = "internal")

  before <- list.dirs(tempdir(), recursive = FALSE)
  expect_error(load_uwot(archive), "Can't find model")
  expect_equal(list.dirs(tempdir(), recursive = FALSE), before)
})

test_that("a failed load releases Annoy indexes opened earlier in the archive", {
  fixture <- tempfile()
  dir.create(fixture)
  old_wd <- setwd(fixture)
  on.exit({
    setwd(old_wd)
    unlink(fixture, recursive = TRUE)
  })

  dir.create("uwot")
  saveRDS(list(metric = list(euclidean = 1:4, manhattan = 1:4)), "uwot/model")
  ann <- methods::new(RcppAnnoy::AnnoyEuclidean, 4L)
  ann$addItem(0L, iris10[1, ])
  ann$build(1L)
  ann$save("uwot/nn1")
  ann$unload()
  utils::tar("incomplete.tar", files = "uwot", tar = "internal")

  before <- list.dirs(tempdir(), recursive = FALSE)
  expect_error(load_uwot("incomplete.tar"), "Can't find nearest neighbor index")
  expect_equal(list.dirs(tempdir(), recursive = FALSE), before)
})

test_that("an extraction failure is reported before reading model files", {
  archive <- tempfile(fileext = ".tar")
  on.exit(unlink(archive), add = TRUE)
  utils::tar(archive, files = character(), tar = "internal")

  # Use a quiet external command to exercise tar's nonzero exit status.
  old_tar <- Sys.getenv("TAR", unset = NA_character_)
  on.exit(
    if (is.na(old_tar)) Sys.unsetenv("TAR") else Sys.setenv(TAR = old_tar),
    add = TRUE
  )
  Sys.setenv(
    TAR = paste(
      shQuote(file.path(R.home("bin"), "Rscript")),
      "--vanilla -e",
      shQuote("quit(status=1L)"),
      "--args"
    )
  )
  before <- list.dirs(tempdir(), recursive = FALSE)
  expect_error(
    suppressWarnings(load_uwot(archive)),
    "Could not extract model archive .*tar exited with status 1"
  )
  expect_equal(list.dirs(tempdir(), recursive = FALSE), before)
})
