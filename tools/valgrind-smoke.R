library(uwot)

# Supply the library containing the build being checked as the first argument.
lib <- normalizePath(commandArgs(trailingOnly = TRUE)[1L], mustWork = TRUE)
loaded <- normalizePath(find.package("uwot"))
stopifnot(identical(loaded, file.path(lib, "uwot")))
message("Testing uwot ", packageVersion("uwot"), " at ", loaded)

local({
  x <- as.matrix(iris[, 1:4])
  file <- tempfile(fileext = ".tar")
  on.exit(unlink(file), add = TRUE)
  model <- umap(
    x,
    nn_method = "annoy",
    n_neighbors = 5,
    n_epochs = 3,
    n_threads = 2,
    n_sgd_threads = 1,
    ret_model = TRUE,
    seed = 42,
    verbose = FALSE
  )
  model <- save_uwot(model, file)
  project <- function(model) {
    umap_transform(
      x[1:4, , drop = FALSE],
      model,
      n_epochs = 3,
      n_threads = 2,
      n_sgd_threads = 1,
      seed = 42,
      verbose = FALSE
    )
  }
  expected <- project(model)
  unload_uwot(model)
  restored <- load_uwot(file)
  on.exit(unload_uwot(restored), add = TRUE)
  actual <- project(restored)
  stopifnot(
    identical(dim(actual), c(4L, 2L)),
    all(is.finite(actual)),
    isTRUE(all.equal(actual, expected))
  )
})
