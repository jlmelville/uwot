library(uwot)
library(RSpectra)
context("similarity graph")

# #96: more convenient way to just get the high dimensional similarity graph

# Hard way first (by using umap function)
# allow for no initialization if n_epochs = 0
# Allowable but returns nothing of value
expect_warning(
  res <- umap(iris10, n_neighbors = 4, init = NULL, n_epochs = 0),
  "will be returned"
)
expect_null(res)

# More sensibly, return high-dimensional data
res <- umap(iris10,
  n_neighbors = 4, init = NULL, n_epochs = 0,
  ret_extra = c("fgraph")
)
expect_is(res, "list")
expect_null(res$embedding)
expect_is(res$fgraph, "sparseMatrix")

# also a bad idea but maybe ok to extract something from the return value
# manually
expect_warning(
  res_with_model <- umap(iris10,
    n_neighbors = 4, init = NULL, n_epochs = 0,
    ret_model = TRUE
  ),
  "will not be valid for transforming"
)
expect_is(res_with_model, "list")
expect_null(res_with_model$embedding)
# but you cannot use umap_transform with the model
expect_error(umap_transform(iris10, res_with_model), "(?i)invalid embedding coordinates")

# Simpler way (by using similarity_graph function)
set.seed(42)
sim_graph <- similarity_graph(iris10, n_neighbors = 4)
expect_equal(res$fgraph, sim_graph)

# can return extra data
set.seed(42)
sim_graph_extra <- similarity_graph(iris10,
  n_neighbors = 4,
  ret_extra = c("sigma", "nn")
)
expect_is(sim_graph_extra, "list")
expect_equal(sim_graph, sim_graph_extra$similarity_graph)
expect_equal(length(sim_graph_extra$sigma), 10)
expect_equal(length(sim_graph_extra$rho), 10)
expect_is(sim_graph_extra$nn$euclidean, "list")

# can use pre-computed nn instead of data
sim_graph_nn <- similarity_graph(
  nn_method = sim_graph_extra$nn,
  ret_extra = c("sigma")
)
expect_equal(sim_graph, sim_graph_nn$similarity_graph)
expect_equal(sim_graph_extra$sigma, sim_graph_nn$sigma)

# can use largevis for t-SNE-like graph
sim_graph_largevis <- similarity_graph(iris10,
  method = "largevis",
  perplexity = 5, ret_extra = c("sigma")
)
expect_is(sim_graph_largevis, "list")
expect_is(sim_graph_largevis$similarity_graph, "sparseMatrix")
expect_equal(dim(sim_graph_largevis$similarity_graph), c(10, 10))
expect_equal(length(sim_graph_extra$sigma), 10)

# specific use case of bbknnR
fss <- fuzzy_simplicial_set(
  nn = sim_graph_extra$nn$euclidean,
  set_op_mix_ratio = 0.5,
  local_connectivity = 2
)

sim_graph_bbknnR <- similarity_graph(
  nn_method = sim_graph_extra$nn$euclidean,
  set_op_mix_ratio = 0.5,
  local_connectivity = 2
)
expect_equal(sim_graph_bbknnR, fss)

# supervised
sim_graphy <- similarity_graph(iris10, n_neighbors = 4, y = ynum)
expect_is(sim_graphy, "sparseMatrix")

# binary edge weights
sim_graphb <- similarity_graph(iris10, n_neighbors = 4, binary_edge_weights = TRUE)
expect_true(all(sim_graphb@x == 1))


test_that("optimize graph layout", {
  iris30 <- iris[c(1:10, 51:60, 101:110), ]
  iris30_sim_graph <- similarity_graph(iris30, n_neighbors = 10)
  set.seed(42)
  iris30_opt <- optimize_graph_layout(iris30_sim_graph, X = iris30)

  set.seed(42)
  iris30_umap <- umap(iris30, n_neighbors = 10)
  expect_equal(iris30_opt, iris30_umap)
})
