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

#include <vector>

#include <Rcpp.h>

#include "RcppPerpendicular.h"
#include "uwot/transform.h"

using namespace Rcpp;

// Initialize embedding as a weighted average of nearest neighbors of each point
// train_embedding: dim x n_train matrix of final embedding coordinates
// nn_index: n_nbrs x n_test matrix of indexes of neighbors in X_train that are
//   nearest neighbors of X_test
// weights: n_nbrs x n_test weight matrix
// Returns the dim x n_test matrix of initialized coordinates.
// [[Rcpp::export]]
NumericMatrix init_transform_parallel(NumericMatrix train_embedding,
                                      IntegerVector nn_index,
                                      std::size_t n_test_vertices,
                                      Nullable<NumericVector> nn_weights,
                                      std::size_t n_threads = 0,
                                      std::size_t grain_size = 1) {

  std::size_t n_train_vertices = train_embedding.ncol();
  std::size_t ndim = train_embedding.nrow();
  std::size_t n_neighbors = nn_index.size() / n_test_vertices;

  auto train_embeddingv = as<std::vector<float>>(train_embedding);
  auto nn_indexv = as<std::vector<int>>(nn_index);
  // Convert to zero-indexing
  for (int &i : nn_indexv) {
    --i;
  }
  std::vector<float> embedding(n_test_vertices * ndim);

  std::vector<float> nn_weightsv(0);
  if (nn_weights.isNotNull()) {
    nn_weightsv = as<std::vector<float>>(nn_weights);
  }
  auto worker = [&](std::size_t begin, std::size_t end) {
    uwot::init_by_mean(begin, end, ndim, n_neighbors, nn_indexv, nn_weightsv,
                       n_test_vertices, train_embeddingv, n_train_vertices,
                       embedding);
  };
  RcppPerpendicular::parallel_for(n_test_vertices, worker, n_threads,
                                  grain_size);

  return NumericMatrix(ndim, n_test_vertices, embedding.begin());
}
