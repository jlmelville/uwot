library(uwot)
library(RSpectra)
context("load/save model")

test_that("can save and load simple model", {
  set.seed(1337)
  model <- umap(iris10,
    n_neighbors = 4, n_epochs = 2, init = "spca",
    metric = "euclidean", verbose = FALSE, n_threads = 0,
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
  model <- umap(jiris10,
    n_neighbors = 4, n_epochs = 2, init = "spca",
    metric = metric2,
    verbose = FALSE, n_threads = 0,
    ret_nn = TRUE, ret_model = TRUE
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
  model <- umap(iris10,
    n_neighbors = 4, n_epochs = 2, init = "spca",
    metric = "euclidean", verbose = FALSE, n_threads = 0,
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

  modelm <- uwot::umap(X, n_neighbors = 4, metric = list("euclidean" = 1:5, "euclidean" = 6:10), ret_model = TRUE)
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
  model <- umap(iris10,
                n_neighbors = 4, n_epochs = 2, init = "spca",
                metric = "correlation", verbose = FALSE, n_threads = 0,
                ret_model = TRUE, ret_extra = c("nn")
  )
  expect_equal(names(model$metric), "correlation")
  expect_equal(model$nn_index$metric, "correlation")

  set.seed(1337)
  transformed_before_reload <-
    umap_transform(iris10,
                   model,
                   n_epochs = 2,
                   ret_extra = c("nn"))

  expect_equal(
    transformed_before_reload$nn$correlation$dist,
    model$nn$correlation$dist,
    check.attributes = FALSE,
    tol = 1e-7
  )

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname, unload = TRUE)

  modelload <- load_uwot(file = mod_fname)

  expect_equal(names(model$metric), "correlation")
  expect_equal(model$nn_index$metric, "correlation")

  set.seed(1337)
  transformed_after_reload <-
    umap_transform(iris10,
                   modelload,
                   n_epochs = 2,
                   ret_extra = c("nn"))
  expect_equal(
    transformed_after_reload$nn$correlation$dist,
    model$nn$correlation$dist,
    check.attributes = FALSE,
    tol = 1e-7
  )

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
  model <- umap(iris10,
    n_neighbors = 4, n_epochs = 2, init = "spca",
    metric = "euclidean", verbose = FALSE, n_threads = 0,
    ret_model = TRUE, ret_extra = c("nn"), nn_method = "hnsw"
  )
  expect_equal(model$nn_method, "hnsw")

  set.seed(1337)
  transformed_before_reload <-
    umap_transform(iris10,
      model,
      n_epochs = 2,
      ret_extra = c("nn")
    )

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname, unload = TRUE)

  modelload <- load_uwot(file = mod_fname)

  expect_equal(modelload$nn_method, "hnsw")

  set.seed(1337)
  transformed_after_reload <-
    umap_transform(iris10,
      modelload,
      n_epochs = 2,
      ret_extra = c("nn")
    )

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
    transformed_after_reload$nn$euclidean$idx,
  )
  expect_equal(
    transformed_before_reload$nn$euclidean$dist,
    transformed_after_reload$nn$euclidean$dist,
    check.attributes = FALSE,
    tol = 1e-7
  )

  expect_equal(
    transformed_before_reload$embedding,
    transformed_after_reload$embedding
  )
})

test_that("save-load nndescent", {
  testthat::skip_if_not_installed("rnndescent")
  set.seed(1337)
  model <- umap(iris10,
                n_neighbors = 4, n_epochs = 2, init = "spca",
                metric = "euclidean", verbose = FALSE, n_threads = 0,
                ret_model = TRUE, ret_extra = c("nn"), nn_method = "nndescent"
  )
  expect_equal(model$nn_method, "nndescent")

  set.seed(1337)
  transformed_before_reload <-
    umap_transform(iris10,
                   model,
                   n_epochs = 2,
                   ret_extra = c("nn")
    )

  mod_fname <- tempfile(tmpdir = tempdir())
  model <- save_uwot(model, file = mod_fname, unload = TRUE)

  modelload <- load_uwot(file = mod_fname)

  expect_equal(modelload$nn_method, "nndescent")

  set.seed(1337)
  transformed_after_reload <-
    umap_transform(iris10,
                   modelload,
                   n_epochs = 2,
                   ret_extra = c("nn")
    )

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
    transformed_after_reload$nn$euclidean$idx,
  )
  expect_equal(
    transformed_before_reload$nn$euclidean$dist,
    transformed_after_reload$nn$euclidean$dist,
    check.attributes = FALSE,
    tol = 1e-7
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
    umap_transform(iris10,
                   modelload2,
                   n_epochs = 2,
                   ret_extra = c("nn")
    )
  expect_equal(
    transformed_after_reload$nn$euclidean$idx,
    transformed_after_reload2$nn$euclidean$idx,
  )
  expect_equal(
    transformed_after_reload$nn$euclidean$dist,
    transformed_after_reload2$nn$euclidean$dist,
    check.attributes = FALSE,
    tol = 1e-7
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
  model <- umap(iris10,
                n_neighbors = 4, n_epochs = 2, init = "spca",
                metric = "euclidean", verbose = FALSE, n_threads = 0,
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
