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

// [[Rcpp::export]]
NumericMatrix init_transform_av_parallel(NumericMatrix train_embedding,
                                         IntegerMatrix nn_index,
                                         std::size_t n_threads = 0,
                                         std::size_t grain_size = 1) {

  std::size_t n_train_vertices = train_embedding.nrow();
  std::size_t ndim = train_embedding.ncol();
  std::size_t n_test_vertices = nn_index.nrow();

  auto train_embeddingv = as<std::vector<float>>(train_embedding);
  auto nn_indexv = as<std::vector<int>>(nn_index);
  // Convert to zero-indexing
  for (int &i : nn_indexv) {
    --i;
  }

  uwot::AverageWorker worker(train_embeddingv, n_train_vertices, nn_indexv,
                             n_test_vertices);

  if (n_threads > 0) {
    RcppPerpendicular::parallel_for(0, n_test_vertices, worker, n_threads,
                                    grain_size);
  } else {
    worker(0, n_test_vertices);
  }

  return NumericMatrix(n_test_vertices, ndim, worker.embedding.begin());
}

// Initialize embedding as a weighted average of nearest neighbors of each point
// train_embedding: n_train x dim matrix of final embedding coordinates
// nn_index: n_test x n_nbrs matrix of indexes of neighbors in X_train that are
//   nearest neighbors of X_test
// weights: n_test x n_nbrs weight matrix
// Returns the n_test x dim matrix of initialized coordinates.
// [[Rcpp::export]]
NumericMatrix init_transform_parallel(NumericMatrix train_embedding,
                                      IntegerMatrix nn_index,
                                      NumericMatrix nn_weights,
                                      std::size_t n_threads = 0,
                                      std::size_t grain_size = 1) {

  std::size_t n_train_vertices = train_embedding.nrow();
  std::size_t ndim = train_embedding.ncol();
  std::size_t n_test_vertices = nn_index.nrow();

  auto train_embeddingv = as<std::vector<float>>(train_embedding);
  auto nn_indexv = as<std::vector<int>>(nn_index);
  // Convert to zero-indexing
  for (int &i : nn_indexv) {
    --i;
  }
  auto nn_weightsv = as<std::vector<float>>(nn_weights);

  uwot::WeightedAverageWorker worker(train_embeddingv, n_train_vertices,
                                     nn_indexv, nn_weightsv, n_test_vertices);

  if (n_threads > 0) {
    RcppPerpendicular::parallel_for(0, n_test_vertices, worker, n_threads,
                                    grain_size);
  } else {
    worker(0, n_test_vertices);
  }

  return NumericMatrix(n_test_vertices, ndim, worker.embedding.begin());
}
