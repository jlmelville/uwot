//  UWOT -- An R package for dimensionality reduction using UMAP
//
//  Copyright (C) 2018 James Melville
//
//  This file is part of UWOT
//
//  UWOT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  UWOT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with UWOT.  If not, see <http://www.gnu.org/licenses/>.

#include <algorithm>
#include <vector>

#include <Rcpp.h>

#include "pforr.h"
#include "uwot/transform.h"

using namespace Rcpp;

// Initialize embedding as a weighted average of nearest neighbors of each point
// train_embedding: dim x n_train matrix of final embedding coordinates
// nn_index: n_nbrs x n_test matrix of indexes of neighbors in X_train that are
//   nearest neighbors of X_test
// weights: n_nbrs x n_test weight matrix
// Returns the dim x n_test matrix of initialized coordinates.
// [[Rcpp::export]]
NumericMatrix
init_transform_parallel(NumericMatrix train_embedding, IntegerVector nn_index,
                        IntegerVector nn_ptr, std::size_t n_test_vertices,
                        Nullable<NumericVector> nn_weights,
                        std::size_t n_threads = 0, std::size_t grain_size = 1) {

  if (nn_weights.isNotNull() &&
      NumericVector(nn_weights).size() != nn_index.size()) {
    stop("Expected one weight per neighbor");
  }
  std::size_t ndim = train_embedding.nrow();
  auto ptr = as<std::vector<std::size_t>>(nn_ptr);
  if (ptr.size() == 1) {
    if (ptr[0] == 0 || nn_index.size() / ptr[0] != n_test_vertices ||
        nn_index.size() % ptr[0] != 0) {
      stop("Invalid fixed neighbor count");
    }
  } else if (ptr.size() != n_test_vertices + 1 || ptr.front() != 0 ||
             ptr.back() != static_cast<std::size_t>(nn_index.size()) ||
             !std::is_sorted(ptr.begin(), ptr.end()) ||
             std::adjacent_find(ptr.begin(), ptr.end()) != ptr.end()) {
    stop("Invalid neighbor column pointers");
  }

  auto train_embeddingv = as<std::vector<float>>(train_embedding);
  auto nn_indexv = as<std::vector<int>>(nn_index);
  // Convert to zero-indexing
  for (int &i : nn_indexv) {
    if (i < 1 || i > train_embedding.ncol())
      stop("Invalid training neighbor index");
    --i;
  }
  std::vector<float> embedding(n_test_vertices * ndim);

  std::vector<float> nn_weightsv(0);
  if (nn_weights.isNotNull()) {
    nn_weightsv = as<std::vector<float>>(nn_weights);
  }
  auto worker = [&](std::size_t begin, std::size_t end) {
    uwot::init_by_mean(begin, end, ndim, ptr, nn_indexv, nn_weightsv,
                       train_embeddingv, embedding);
  };
  pforr::parallel_for(n_test_vertices, worker, n_threads, grain_size);

  return NumericMatrix(ndim, n_test_vertices, embedding.begin());
}
