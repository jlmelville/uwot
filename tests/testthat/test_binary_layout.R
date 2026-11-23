test_that("layout with binary weights matches a manually binarized graph", {
  graph <- similarity_graph(iris[1:20, 1:4], n_neighbors = 4, n_threads = 0)
  binary <- graph
  binary@x[] <- 1

  layout <- function(graph, binary_edge_weights) {
    set.seed(42)
    coords <- optimize_graph_layout(
      graph,
      binary_edge_weights = binary_edge_weights,
      init = "random",
      n_epochs = 20,
      n_sgd_threads = 1,
      verbose = FALSE
    )
    as.vector(coords)
  }

  expected <- layout(binary, FALSE)
  expect_false(identical(layout(graph, FALSE), expected))
  expect_equal(layout(graph, TRUE), expected)
})
