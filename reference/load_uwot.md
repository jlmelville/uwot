# Save or Load a Model

Functions to write a UMAP model to a file, and to restore.

## Usage

``` r
load_uwot(file, verbose = FALSE)
```

## Arguments

- file:

  name of the file where the model is to be saved or read from.

- verbose:

  if `TRUE`, log information to the console.

## Value

The model saved at `file`, for use with
[`umap_transform`](https://jlmelville.github.io/uwot/reference/umap_transform.md).
Additionally, it contains an extra item: `mod_dir`, which contains the
path to the temporary working directory used during loading of the
model. This directory cannot be removed until this model has been
unloaded by using
[`unload_uwot`](https://jlmelville.github.io/uwot/reference/unload_uwot.md).

## See also

[`save_uwot`](https://jlmelville.github.io/uwot/reference/save_uwot.md),
[`unload_uwot`](https://jlmelville.github.io/uwot/reference/unload_uwot.md)

## Examples

``` r
library(RSpectra)

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
