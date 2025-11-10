# Merge Similarity Graph by Simplicial Set Intersection

Combine two similarity graphs by treating them as fuzzy topological sets
and forming the intersection.

## Usage

``` r
simplicial_set_intersect(x, y, weight = 0.5, n_threads = NULL, verbose = FALSE)
```

## Arguments

- x:

  A sparse matrix representing the first similarity graph in the
  intersection operation.

- y:

  A sparse matrix representing the second similarity graph in the
  intersection operation.

- weight:

  A value between `0 - 1`, controlling the relative influence of `x` and
  `y` in the intersection. Default (`0.5`) gives equal influence. Values
  smaller than `0.5` put more weight on `x`. Values greater than `0.5`
  put more weight on `y`.

- n_threads:

  Number of threads to use when resetting the local metric. Default is
  half the number of concurrent threads supported by the system.

- verbose:

  If `TRUE`, log progress to the console.

## Value

A sparse matrix containing the intersection of `x` and `y`.

## Examples

``` r
# Form two different "views" of the same data
iris30 <- iris[c(1:10, 51:60, 101:110), ]
iris_sg12 <- similarity_graph(iris30[, 1:2], n_neighbors = 5)
iris_sg34 <- similarity_graph(iris30[, 3:4], n_neighbors = 5)

# Combine the two representations into one
iris_combined <- simplicial_set_intersect(iris_sg12, iris_sg34)

# Optimize the layout based on the combined view
iris_combined_umap <- optimize_graph_layout(iris_combined, n_epochs = 100)
```
