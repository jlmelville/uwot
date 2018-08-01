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
  simplicial_set, target, unknown_dist = 1.0, far_dist = 5.0, verbose = FALSE) {

  # Convert to dgTMatrix to get to the j indices
  simplicial_set <- methods::as(simplicial_set, "dgTMatrix")
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
reset_local_connectivity <- function(simplicial_set) {
  fuzzy_set_union(row_max_normalize(simplicial_set))
}

# Under the assumption of categorical distance for the intersecting
# simplicial set perform a fast intersection.
# This is not at all fast in R, use fast_intersection_cpp instead
fast_intersection <- function(rows, cols, values, target, unknown_dist = 1.0,
                              far_dist = 5.0) {
  ex_unknown <- exp(-unknown_dist)
  ex_far <- exp(-far_dist)

  for (nz in 1:length(values)) {
    i <- rows[nz]
    j <- cols[nz]
    if (is.na(target[i]) || is.na(target[j])) {
      values[nz] <- values[nz] * ex_unknown
    }
    else if (target[i] != target[j]) {
      values[nz] <- values[nz] * ex_far
    }
  }

  values
}

# column maximums of a dgCMatrix
colMaxs <- function(X) {
  nr <- nrow(X)
  result <- rep(0, nr)

  dX <- diff(X@p)

  for (i in 1:nr) {
    if (dX[i] > 0) {
      result[i] <- max(X@x[(X@p[i] + 1):X@p[i + 1]])
    }
  }

  result
}

# row maximums of a dgCMatrix
rowMaxs <- function(X) {
  colMaxs(Matrix::t(X))
}
