# Merge Similarity Graph by Simplicial Set Union

Combine two similarity graphs by treating them as fuzzy topological sets
and forming the union.

## Usage

``` r
simplicial_set_union(x, y, n_threads = NULL, verbose = FALSE)
```

## Arguments

- x:

  A sparse matrix representing the first similarity graph in the union
  operation.

- y:

  A sparse matrix representing the second similarity graph in the union
  operation.

- n_threads:

  Number of threads to use when resetting the local metric. Default is
  half the number of concurrent threads supported by the system.

- verbose:

  If `TRUE`, log progress to the console.

## Value

A sparse matrix containing the union of `x` and `y`.

## Examples

``` r
# Form two different "views" of the same data
iris30 <- iris[c(1:10, 51:60, 101:110), ]
iris_sg12 <- similarity_graph(iris30[, 1:2], n_neighbors = 5)
iris_sg34 <- similarity_graph(iris30[, 3:4], n_neighbors = 5)

# Combine the two representations into one
iris_combined <- simplicial_set_union(iris_sg12, iris_sg34)

# Optimize the layout based on the combined view
iris_combined_umap <- optimize_graph_layout(iris_combined, n_epochs = 100)
```
