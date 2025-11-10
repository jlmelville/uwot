# Save or Load a Model

Functions to write a UMAP model to a file, and to restore.

## Usage

``` r
save_uwot(model, file, unload = FALSE, verbose = FALSE)
```

## Arguments

- model:

  a UMAP model create by
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.md).

- file:

  name of the file where the model is to be saved or read from.

- unload:

  if `TRUE`, unload all nearest neighbor indexes for the model. The
  `model` will no longer be valid for use in
  [`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md)
  and the temporary working directory used during model saving will be
  deleted. You will need to reload the model with `load_uwot` to use the
  model. If `FALSE`, then the model can be re-used without reloading,
  but you must manually unload the NN index when you are finished using
  it if you want to delete the temporary working directory. To unload
  manually, use
  [`unload_uwot`](https://jlmelville.github.io/uwot/reference/unload_uwot.md).
  The absolute path of the working directory is found in the `mod_dir`
  item of the return value.

- verbose:

  if `TRUE`, log information to the console.

## Value

`model` with one extra item: `mod_dir`, which contains the path to the
working directory. If `unload = FALSE` then this directory still exists
after this function returns, and can be cleaned up with
[`unload_uwot`](https://jlmelville.github.io/uwot/reference/unload_uwot.md).
If you don't care about cleaning up this directory, or `unload = TRUE`,
then you can ignore the return value.

## See also

[`load_uwot`](https://jlmelville.github.io/uwot/reference/load_uwot.md),
[`unload_uwot`](https://jlmelville.github.io/uwot/reference/unload_uwot.md)

## Examples

``` r
iris_train <- iris[c(1:10, 51:60), ]
iris_test <- iris[100:110, ]

# create model
model <- umap(iris_train, ret_model = TRUE, n_epochs = 20)

# save without unloading: this leaves behind a temporary working directory
model_file <- tempfile("iris_umap")
model <- save_uwot(model, file = model_file)

# The model can continue to be used
test_embedding <- umap_transform(iris_test, model)

# To manually unload the model from memory when finished and to clean up
# the working directory (this doesn't touch your model file)
unload_uwot(model)

# At this point, model cannot be used with umap_transform, this would fail:
# test_embedding2 <- umap_transform(iris_test, model)

# restore the model: this also creates a temporary working directory
model2 <- load_uwot(file = model_file)
test_embedding2 <- umap_transform(iris_test, model2)

# Unload and clean up the loaded model temp directory
unload_uwot(model2)

# clean up the model file
unlink(model_file)

# save with unloading: this deletes the temporary working directory but
# doesn't allow the model to be re-used
model3 <- umap(iris_train, ret_model = TRUE, n_epochs = 20)
model_file3 <- tempfile("iris_umap")
model3 <- save_uwot(model3, file = model_file3, unload = TRUE)
```
