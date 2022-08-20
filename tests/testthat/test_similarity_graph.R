library(uwot)
context("similarity graph")

# #96: more convenient way to just get the high dimensional similarity graph

# Hard way first (by using umap function)
# allow for no initialization if n_epochs = 0
# Allowable but returns nothing of value
expect_warning(res <- umap(iris10, n_neighbors = 4, init = NULL, n_epochs = 0),
               "will be returned")
expect_null(res)

# More sensibly, return high-dimensional data
res <- umap(iris10, n_neighbors = 4, init = NULL, n_epochs = 0,
            ret_extra = c("fgraph"))
expect_is(res, "list")
expect_null(res$embedding)
expect_is(res$fgraph, "sparseMatrix")

# also a bad idea but maybe ok to extract something from the return value
# manually
expect_warning(res_with_model <- umap(iris10, n_neighbors = 4, init = NULL, n_epochs = 0,
                           ret_model = TRUE),
               "will not be valid for transforming")
expect_is(res_with_model, "list")
expect_null(res_with_model$embedding)
# but you cannot use umap_transform with the model
expect_error(umap_transform(iris10, res_with_model), "(?i)invalid embedding coordinates")
