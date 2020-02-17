library(uwot)
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
