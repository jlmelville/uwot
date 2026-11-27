library(uwot)
library(RSpectra)

test_that("can save and load simple model", {
  set.seed(1337)
  model <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = "euclidean",
    verbose = FALSE,
    n_threads = 0,
    ret_model = TRUE
  )

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname)
  expect_true(file.exists(mod_fname))

  # Can use model after saving
  set.seed(1337)
  res_trans <- umap_transform(iris10, model)
  expect_ok_matrix(res_trans)

  # Clean up temp dir from saving
  expect_true(file.exists(model$mod_dir))
  unload_uwot(model)
  expect_false(file.exists(model$mod_dir))
  # Can't use transform now model is unloaded
  expect_error(umap_transform(iris10, model), "is unloaded")

  modelload <- load_uwot(file = mod_fname)
  set.seed(1337)
  resload_trans <- umap_transform(iris10, modelload)
  expect_ok_matrix(resload_trans)

  expect_equal(resload_trans, res_trans)
  if (file.exists(mod_fname)) {
    unlink(mod_fname)
  }

  # Clean up temp dir from loading
  expect_true(file.exists(modelload$mod_dir))
  unload_uwot(modelload)
  expect_false(file.exists(modelload$mod_dir))
})


test_that("can save and load mixed distance model", {
  set.seed(1337)
  jiris10 <- jitter(iris10)
  metric2 <- list(
    "euclidean" = c(1, 2),
    "cosine" = c("Petal.Length", "Petal.Width")
  )
  model <- umap(
    jiris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = metric2,
    verbose = FALSE,
    n_threads = 0,
    ret_nn = TRUE,
    ret_model = TRUE
  )

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname)
  expect_true(file.exists(mod_fname))

  # Can use model after saving
  set.seed(1337)
  res_trans <- umap_transform(jiris10, model)
  expect_ok_matrix(res_trans)

  # Clean up temp dir from saving
  expect_true(file.exists(model$mod_dir))
  unload_uwot(model)
  expect_false(file.exists(model$mod_dir))
  # Can't use transform now model is unloaded
  expect_error(umap_transform(iris10, model), "is unloaded")

  modelload <- load_uwot(file = mod_fname)
  set.seed(1337)
  resload_trans <- umap_transform(jiris10, modelload)
  expect_ok_matrix(resload_trans)

  expect_equal(resload_trans, res_trans)
  if (file.exists(mod_fname)) {
    unlink(mod_fname)
  }

  # Clean up temp dir from loading
  expect_true(file.exists(modelload$mod_dir))
  unload_uwot(modelload)
  expect_false(file.exists(modelload$mod_dir))
})

test_that("unloading a model on save", {
  set.seed(1337)
  model <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = "euclidean",
    verbose = FALSE,
    n_threads = 0,
    ret_model = TRUE
  )

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname, unload = TRUE)
  expect_false(file.exists(model$mod_dir))

  # Trying to transform with a model that got unloaded won't work
  expect_error(umap_transform(iris10, model), "is unloaded")

  modelload <- load_uwot(file = mod_fname)
  # Clean up temp dir from loading
  expect_true(file.exists(modelload$mod_dir))
  # Can avoid cleaning up if you really want that
  unload_uwot(modelload, cleanup = FALSE)
  expect_true(file.exists(modelload$mod_dir))
  # Can unload multiple times
  unload_uwot(modelload, cleanup = TRUE)
  expect_false(file.exists(modelload$mod_dir))
})

# #88
test_that("save-load-save", {
  set.seed(1337)
  X <- matrix(rnorm(100), 10, 10)

  model <- uwot::umap(X, n_neighbors = 4, ret_model = TRUE)
  model_file <- tempfile(tmpdir = tempdir())
  model <- uwot::save_uwot(model, file = model_file)
  model2 <- uwot::load_uwot(file = model_file)
  new_file <- tempfile(tmpdir = tempdir())
  uwot::save_uwot(model2, file = new_file)
  expect_true(file.exists(new_file))

  modelm <- uwot::umap(
    X,
    n_neighbors = 4,
    metric = list("euclidean" = 1:5, "euclidean" = 6:10),
    ret_model = TRUE
  )
  modelm_file <- tempfile(tmpdir = tempdir())
  modelm <- uwot::save_uwot(modelm, file = modelm_file)
  modelm2 <- uwot::load_uwot(file = modelm_file)
  new_filem <- tempfile(tmpdir = tempdir())
  uwot::save_uwot(modelm2, file = new_filem)
  expect_true(file.exists(new_filem))
})

# #117 correlation metric not correctly restored
test_that("reload-correlation", {
  set.seed(1337)
  model <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = "correlation",
    verbose = FALSE,
    n_threads = 0,
    ret_model = TRUE,
    ret_extra = c("nn")
  )
  expect_equal(names(model$metric), "correlation")
  expect_equal(model$nn_index$metric, "correlation")

  set.seed(1337)
  transformed_before_reload <-
    umap_transform(iris10, model, n_epochs = 2, ret_extra = c("nn"))

  expect_angular_distances(
    transformed_before_reload$nn$correlation$dist,
    model$nn$correlation$dist,
    tolerance = 1e-7
  )

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname, unload = TRUE)

  modelload <- load_uwot(file = mod_fname)

  expect_equal(names(modelload$metric), "correlation")
  expect_equal(modelload$nn_index$metric, "correlation")

  set.seed(1337)
  transformed_after_reload <-
    umap_transform(iris10, modelload, n_epochs = 2, ret_extra = c("nn"))
  expect_angular_distances(
    transformed_after_reload$nn$correlation$dist,
    model$nn$correlation$dist,
    tolerance = 1e-7
  )
  expect_identical(transformed_after_reload$nn, transformed_before_reload$nn)

  if (file.exists(mod_fname)) {
    unlink(mod_fname)
  }
  expect_true(file.exists(modelload$mod_dir))
  unload_uwot(modelload)
  expect_false(file.exists(modelload$mod_dir))
})

test_that("save-load hnsw", {
  testthat::skip_if_not_installed("RcppHNSW")
  set.seed(1337)
  model <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = "euclidean",
    verbose = FALSE,
    n_threads = 0,
    ret_model = TRUE,
    ret_extra = c("nn"),
    nn_method = "hnsw"
  )
  expect_equal(model$nn_method, "hnsw")

  set.seed(1337)
  transformed_before_reload <-
    umap_transform(iris10, model, n_epochs = 2, ret_extra = c("nn"))

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname, unload = TRUE)

  modelload <- load_uwot(file = mod_fname)

  expect_equal(modelload$nn_method, "hnsw")

  set.seed(1337)
  transformed_after_reload <-
    umap_transform(iris10, modelload, n_epochs = 2, ret_extra = c("nn"))

  if (file.exists(mod_fname)) {
    unlink(mod_fname)
  }
  expect_true(file.exists(modelload$mod_dir))
  unload_uwot(modelload)
  expect_false(file.exists(modelload$mod_dir))

  expect_equal(model$nn$euclidean$idx, modelload$nn$euclidean$idx)
  expect_equal(model$nn$euclidean$dist, modelload$nn$euclidean$dist)

  expect_equal(
    transformed_before_reload$nn$euclidean$idx,
    transformed_after_reload$nn$euclidean$idx
  )
  expect_equal(
    transformed_before_reload$nn$euclidean$dist,
    transformed_after_reload$nn$euclidean$dist,
    tolerance = 1e-7
  )

  expect_equal(
    transformed_before_reload$embedding,
    transformed_after_reload$embedding
  )
})

test_that("save-load nndescent", {
  testthat::skip_if_not_installed("rnndescent")
  set.seed(1337)
  model <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = "euclidean",
    verbose = FALSE,
    n_threads = 0,
    ret_model = TRUE,
    ret_extra = c("nn"),
    nn_method = "nndescent"
  )
  expect_equal(model$nn_method, "nndescent")

  set.seed(1337)
  transformed_before_reload <-
    umap_transform(iris10, model, n_epochs = 2, ret_extra = c("nn"))

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname, unload = TRUE)

  modelload <- load_uwot(file = mod_fname)

  expect_equal(modelload$nn_method, "nndescent")

  set.seed(1337)
  transformed_after_reload <-
    umap_transform(iris10, modelload, n_epochs = 2, ret_extra = c("nn"))

  if (file.exists(mod_fname)) {
    unlink(mod_fname)
  }
  expect_true(file.exists(modelload$mod_dir))
  unload_uwot(modelload)
  expect_false(file.exists(modelload$mod_dir))

  expect_equal(model$nn$euclidean$idx, modelload$nn$euclidean$idx)
  expect_equal(model$nn$euclidean$dist, modelload$nn$euclidean$dist)

  expect_equal(
    transformed_before_reload$nn$euclidean$idx,
    transformed_after_reload$nn$euclidean$idx
  )
  expect_equal(
    transformed_before_reload$nn$euclidean$dist,
    transformed_after_reload$nn$euclidean$dist,
    tolerance = 1e-7
  )

  expect_equal(
    transformed_before_reload$embedding,
    transformed_after_reload$embedding
  )

  mod_fname2 <- tempfile(tmpdir = tempdir())
  saveRDS(modelload, mod_fname2)
  modelload2 <- readRDS(mod_fname2)
  expect_equal(modelload2$nn_method, "nndescent")
  set.seed(1337)
  transformed_after_reload2 <-
    umap_transform(iris10, modelload2, n_epochs = 2, ret_extra = c("nn"))
  expect_equal(
    transformed_after_reload$nn$euclidean$idx,
    transformed_after_reload2$nn$euclidean$idx
  )
  expect_equal(
    transformed_after_reload$nn$euclidean$dist,
    transformed_after_reload2$nn$euclidean$dist,
    tolerance = 1e-7
  )
  expect_equal(
    transformed_after_reload$embedding,
    transformed_after_reload2$embedding
  )
  if (file.exists(mod_fname2)) {
    unlink(mod_fname2)
  }
})

# 131: can't use a relative path for saving
test_that("save-load relative path", {
  set.seed(1337)
  model <- umap(
    iris10,
    n_neighbors = 4,
    n_epochs = 2,
    init = "spca",
    metric = "euclidean",
    verbose = FALSE,
    n_threads = 0,
    ret_model = TRUE
  )
  # remember to go back to the original working directory
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  # move to a temp directory for this test
  test_dir <- tempfile("test_relative_path_")
  dir.create(test_dir)
  setwd(test_dir)

  # Test 1: relative path no sub-folders
  rel_path <- "model.uwot"
  model <- save_uwot(model, file = rel_path)
  expect_true(file.exists(rel_path))

  modelload <- load_uwot(file = rel_path)
  set.seed(1337)
  resload_trans <- umap_transform(iris10, modelload)
  expect_ok_matrix(resload_trans)

  # Test 2: Relative path within a sub-folder
  dir.create("my_folder")
  rel_sub_path <- file.path("my_folder", "model.uwot")
  model <- save_uwot(model, file = rel_sub_path)
  expect_true(file.exists(rel_sub_path))

  modelload <- load_uwot(file = rel_sub_path)
  set.seed(1337)
  resload_trans <- umap_transform(iris10, modelload)
  expect_ok_matrix(resload_trans)
})

for (backend in c("annoy", "hnsw", "nndescent")) {
  test_that(paste("repeated saves replace staging files for", backend), {
    if (backend == "hnsw") {
      skip_if_not_installed("RcppHNSW")
    }
    if (backend == "nndescent") {
      skip_if_not_installed("rnndescent")
    }
    archives <- c(tempfile(), tempfile())
    on.exit(unlink(archives))
    blocks <- list(euclidean = 1:2, euclidean = 3:4)
    m <- fit_small_model(nn_method = backend, metric = blocks)
    query <- function(m) {
      umap_transform(
        iris10,
        m,
        ret_extra = "nn",
        n_epochs = 0,
        n_threads = 0
      )$nn
    }
    expected <- query(m)
    m <- save_uwot(m, archives[1])
    previous <- m$mod_dir
    m <- save_uwot(m, archives[2])
    expect_false(dir.exists(previous))
    expect_true(all(file.exists(archives)))
    metadata <- readRDS(file.path(m$mod_dir, "uwot", "model"))
    expect_null(metadata$nn_index)
    expect_equal(query(m), expected)
    unload_uwot(m)
    expect_false(dir.exists(m$mod_dir))
    for (archive in archives) {
      restored <- load_uwot(archive)
      expect_equal(query(restored), expected)
      unload_uwot(restored)
    }
  })
}

test_that("a native index write failure preserves the archive and live model", {
  archive <- tempfile()
  on.exit(unlink(archive))
  for (backend in c("annoy", "hnsw")) {
    if (backend == "hnsw" && !requireNamespace("RcppHNSW", quietly = TRUE)) {
      next
    }
    m <- fit_small_model(nn_method = backend)
    expected <- query_neighbors(m)
    # Direct the real backend to a missing parent inside the staging folder.
    failed <- m
    native <- m$nn_index$ann
    failed$nn_index$ann <- list(
      save = function(file) {
        native$save(file.path(dirname(file), "absent", "nn1"))
      }
    )
    if (backend == "annoy") {
      failed$nn_index$ann$getNTrees <- native$getNTrees
    }
    expect_error(
      save_uwot(failed, archive, unload = TRUE),
      "Could not save nearest neighbor index|[Uu]nable to open|[Cc]annot open"
    )
    expect_false(file.exists(archive))
    expect_equal(query_neighbors(m), expected)
    unload_uwot(m, cleanup = FALSE)
  }
})

test_that("saving inside the old working directory preserves the new archive", {
  first <- tempfile()
  on.exit(unlink(first))
  m <- save_uwot(fit_small_model(), first)
  old_dir <- m$mod_dir
  on.exit(unlink(old_dir, recursive = TRUE), add = TRUE)
  second <- file.path(old_dir, "backup.tar")
  saved <- save_uwot(m, second, unload = TRUE)
  expect_false(dir.exists(saved$mod_dir))
  expect_identical(list.files(old_dir), "backup.tar")
  restored <- load_uwot(second)
  expect_equal(query_neighbors(restored)$idx[1, 1], 1L)
  unload_uwot(restored)
})

test_that("an Annoy write diagnostic rejects a nonempty staged file", {
  m <- fit_small_model()
  expected <- query_neighbors(m)
  archive <- tempfile()
  on.exit(unload_uwot(m, cleanup = FALSE))
  failed <- m
  proxy <- new.env(parent = emptyenv())
  proxy$getNTrees <- m$nn_index$ann$getNTrees
  proxy$save <- function(path) {
    writeBin(as.raw(1:12), path)
    # Use the real RcppAnnoy diagnostic while a partial staged file exists.
    m$nn_index$ann$save(file.path(dirname(path), "absent", "nn1"))
  }
  failed$nn_index$ann <- proxy
  outer <- utils::capture.output(
    {
      expect_error(save_uwot(failed, archive, unload = TRUE), "Unable to open")
      message("outer sink restored")
    },
    type = "message"
  )
  expect_true("outer sink restored" %in% outer)
  expect_false(file.exists(archive))
  expect_equal(query_neighbors(m), expected)
})

test_that("HNSW saves reject a partial header, base table or final record", {
  skip_if_not_installed("RcppHNSW")
  m <- fit_small_model(nn_method = "hnsw")
  expected <- query_neighbors(m)
  archive <- tempfile()
  m <- save_uwot(m, archive)
  on.exit({
    unload_uwot(m)
    unlink(archive)
  })
  path <- file.path(m$mod_dir, "uwot", "nn1")
  index <- readBin(path, "raw", n = file.size(path))
  expect_true(hnsw_index_complete(path))

  failed <- m
  failed$nn_index$ann <- list(save = function(path) {
    writeBin(index[seq_len(n)], path)
  })
  for (n in c(1L, length(index) %/% 2L, length(index) - 1L)) {
    dest <- tempfile()
    expect_error(
      save_uwot(failed, dest, unload = TRUE),
      "Could not save nearest neighbor index"
    )
    expect_false(file.exists(dest))
    expect_true(file.exists(path))
  }
  expect_equal(query_neighbors(m), expected)
})
