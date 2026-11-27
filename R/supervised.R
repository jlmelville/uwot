# Combine a fuzzy simplicial set with another fuzzy simplicial set
# generated from categorical data using categorical distances. The target
# data is assumed to be categorical label data (a vector of labels),
# and this will update the fuzzy simplicial set to respect that label data.
# TODO: optional category cardinality based weighting of distance
# simplicial_set The input fuzzy simplicial set.
# target The categorical labels to use in the intersection.
# unknown_dist The distance an unknown label (-1) is assumed to be from any point.
# far_dist The distance between unmatched labels.
# Return The resulting intersected fuzzy simplicial set.
categorical_simplicial_set_intersection <- function(
  simplicial_set,
  target,
  unknown_dist = 1.0,
  far_dist = 5.0,
  verbose = FALSE
) {
  # Convert to dgTMatrix to get to the j indices
  simplicial_set <- methods::as(simplicial_set, "TsparseMatrix")
  simplicial_set@x <- fast_intersection_cpp(
    simplicial_set@i,
    simplicial_set@j,
    simplicial_set@x,
    target,
    unknown_dist,
    far_dist
  )

  # drop0 converts back to dgCMatrix
  reset_local_connectivity(Matrix::drop0(simplicial_set))
}

# Reset the local connectivity requirement -- each data sample should
# have complete confidence in at least one 1-simplex in the simplicial set.
# We can enforce this by locally rescaling confidences, and then remerging the
# different local simplicial sets together.
reset_local_connectivity <-
  function(
    simplicial_set,
    reset_local_metric = FALSE,
    num_local_metric_neighbors = 15,
    n_threads = NULL,
    verbose = FALSE
  ) {
    # Python UMAP stores graph as CSR uwot uses CSC so need to be careful about
    # which axis to normalize
    simplicial_set <- col_max_normalize(simplicial_set)
    if (reset_local_metric) {
      if (is.null(n_threads)) {
        n_threads <- default_num_threads()
      }
      tsmessage(
        "Resetting local metric",
        pluralize("thread", n_threads, " using")
      )
      metric_res <-
        reset_local_metrics_parallel(
          simplicial_set@p,
          simplicial_set@x,
          num_local_metric_neighbors = num_local_metric_neighbors,
          n_threads = n_threads
        )
      simplicial_set@x <- metric_res$values
      # TODO: at least some failures are very typical and it doesn't seem to
      # affect results, so not worth reporting this for now.
      # if (metric_res$n_failures > 0) {
      #   tsmessage(metric_res$n_failures, " local metric reset failures")
      # }
    }
    fuzzy_set_union(simplicial_set)
  }

general_simplicial_set_intersection <- function(left, right, weight) {
  left <- as_neighbor_graph(left)
  right <- as_neighbor_graph(right)
  result <- methods::as(left + right, "TsparseMatrix")

  result@x <- general_sset_intersection_cpp(
    left@p,
    left@i,
    left@x,
    right@p,
    right@i,
    right@x,
    result@i,
    result@j,
    result@x,
    weight
  )

  result
}

# Sparse Matrix functions -------------------------------------------------

# normalize each column of a dgCMatrix by its maximum
# https://stackoverflow.com/questions/39284774/column-rescaling-for-a-very-large-sparse-matrix-in-r
col_max_normalize <- function(X) {
  X@x <- X@x / rep.int(colMaxs(X), diff(X@p))
  X
}

# column maximums of a dgCMatrix
colMaxs <- function(X) {
  ptr <- X@p
  xs <- X@x
  vapply(
    1:ncol(X),
    function(i) {
      if (ptr[i + 1] > ptr[i]) {
        max(xs[(ptr[i] + 1):ptr[i + 1]])
      } else {
        0
      }
    },
    numeric(1)
  )
}
