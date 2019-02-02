library(uwot)
context("load/save model")

test_that("can save and load simple model", {
  set.seed(1337)
  model <- umap(iris10, n_neighbors = 4, n_epochs = 2, init = "spca",
              metric = "euclidean", verbose = FALSE, n_threads = 0, 
              ret_model = TRUE)
  
  mod_fname <- tempfile()
  save_uwot(model, file = mod_fname)
  expect_true(file.exists(mod_fname))
  
  # Can use model after saving
  set.seed(1337)
  res_trans <- umap_transform(iris10, model)
  expect_ok_matrix(res_trans)
  
  modelload <- load_uwot(file = mod_fname)
  set.seed(1337)
  resload_trans <- umap_transform(iris10, modelload)
  expect_ok_matrix(resload_trans)
  
  expect_equal(resload_trans, res_trans)
  if (file.exists(mod_fname)) {
    unlink(mod_fname)
  }
})


test_that("can save and load mixed distance model", {
  set.seed(1337)
  jiris10 <- jitter(iris10)
  metric2 <- list("euclidean" = c(1, 2),
                  "cosine" = c("Petal.Length", "Petal.Width"))
  model <- umap(jiris10, n_neighbors = 4, n_epochs = 2, init = "spca",
                  metric = metric2,
                  verbose = FALSE, n_threads = 0,
                  ret_nn = TRUE, ret_model = TRUE)
  
  mod_fname <- tempfile()
  save_uwot(model, file = mod_fname)
  expect_true(file.exists(mod_fname))
  
  # Can use model after saving
  set.seed(1337)
  res_trans <- umap_transform(jiris10, model)
  expect_ok_matrix(res_trans)
  
  modelload <- load_uwot(file = mod_fname)
  set.seed(1337)
  resload_trans <- umap_transform(jiris10, modelload)
  expect_ok_matrix(resload_trans)
  
  expect_equal(resload_trans, res_trans)
  if (file.exists(mod_fname)) {
    unlink(mod_fname)
  }
})
