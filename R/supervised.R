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
                                                    simplicial_set, target,
                                                    unknown_dist = 1.0,
                                                    far_dist = 5.0,
                                                    verbose = FALSE) {

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

  for (nz in seq_len(length(values))) {
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

general_simplicial_set_intersection <- function(left, right, weight) {
  result <- methods::as(left + right, "dgTMatrix")

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

# An R translation of the Python function. Not very fast,
# so use the C++ version instead
general_sset_intersection <- function(indptr1,
                                      indices1,
                                      data1,
                                      indptr2,
                                      indices2,
                                      data2,
                                      result_row,
                                      result_col,
                                      result_val,
                                      mix_weight = 0.5) {
  left_min <- max(min(data1) / 2.0, 1.0e-8)
  right_min <- max(min(data2) / 2.0, 1.0e-8)

  for (idx in seq_len(length(result_row))) {
    i <- result_col[idx] + 1
    j <- result_row[idx]

    left_val <- left_min
    for (k in (indptr1[i]):(indptr1[i + 1] - 1)) {
      if (indices1[k + 1] == j) {
        left_val <- data1[k + 1]
      }
    }

    right_val <- right_min
    for (k in (indptr2[i]):(indptr2[i + 1] - 1)) {
      if (indices2[k + 1] == j) {
        right_val <- data2[k + 1]
      }
    }

    if (left_val > left_min || right_val > right_min) {
      if (mix_weight < 0.5) {
        result_val[idx] <- left_val *
          right_val^(mix_weight / (1.0 - mix_weight))
      }
      else {
        result_val[idx] <- right_val *
          left_val^(((1.0 - mix_weight) / mix_weight))
      }
    }
  }
  result_val
}


# Sparse Matrix functions -------------------------------------------------

# normalize each column of a dgCMatrix by its maximum
# https://stackoverflow.com/questions/39284774/column-rescaling-for-a-very-large-sparse-matrix-in-r
col_max_normalize <- function(X) {
  X@x <- X@x / rep.int(colMaxs(X), diff(X@p))
  X
}

# normalize each row of a dgCMatrix by its maximum
row_max_normalize <- function(X) {
  Matrix::t(col_max_normalize(Matrix::t(X)))
}

col_sum_normalize <- function(X) {
  X@x <- X@x / rep.int(Matrix::colSums(X), diff(X@p))
  X
}

row_sum_normalize <- function(X) {
  Matrix::t(col_sum_normalize(Matrix::t(X)))
}

# column maximums of a dgCMatrix
colMaxs <- function(X) {
  nc <- ncol(X)
  result <- rep(0, nc)

  dX <- diff(X@p)

  for (i in 1:nc) {
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
