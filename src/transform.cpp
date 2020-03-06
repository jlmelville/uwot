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

struct AverageWorker {

  const std::vector<float> &train_embedding;
  std::size_t n_train_vertices;

  const std::vector<int> &nn_index;
  std::size_t n_test_vertices;

  std::size_t ndim;
  std::size_t n_neighbors;
  std::vector<float> embedding;

  AverageWorker(const std::vector<float> &train_embedding,
                std::size_t n_train_vertices, const std::vector<int> &nn_index,
                std::size_t n_test_vertices)
      : train_embedding(train_embedding), n_train_vertices(n_train_vertices),
        nn_index(nn_index), n_test_vertices(n_test_vertices),
        ndim(train_embedding.size() / n_train_vertices),
        n_neighbors(nn_index.size() / n_test_vertices),
        embedding(n_test_vertices * n_neighbors) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> sumc(ndim);
    for (std::size_t i = begin; i < end; i++) {
      std::fill(sumc.begin(), sumc.end(), 0.0);

      for (std::size_t j = 0; j < n_neighbors; j++) {
        std::size_t nbr = nn_index[i + j * n_test_vertices];
        for (std::size_t k = 0; k < ndim; k++) {
          sumc[k] += train_embedding[nbr + k * n_train_vertices];
        }
      }

      for (std::size_t k = 0; k < ndim; k++) {
        embedding[i + k * n_test_vertices] = sumc[k] / n_neighbors;
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix init_transform_av_parallel(
    Rcpp::NumericMatrix train_embedding, Rcpp::IntegerMatrix nn_index,
    std::size_t n_threads = 0, std::size_t grain_size = 1) {

  std::size_t n_train_vertices = train_embedding.nrow();
  std::size_t ndim = train_embedding.ncol();
  std::size_t n_test_vertices = nn_index.nrow();

  auto train_embeddingv = Rcpp::as<std::vector<float>>(train_embedding);
  auto nn_indexv = Rcpp::as<std::vector<int>>(nn_index);
  // Convert to zero-indexing
  for (int &i : nn_indexv) {
    --i;
  }

  AverageWorker worker(train_embeddingv, n_train_vertices, nn_indexv,
                       n_test_vertices);

  if (n_threads > 0) {
    RcppPerpendicular::parallel_for(0, n_test_vertices, worker, n_threads,
                                    grain_size);
  } else {
    worker(0, n_test_vertices);
  }

  return Rcpp::NumericMatrix(n_test_vertices, ndim, worker.embedding.begin());
}

struct WeightedAverageWorker {

  const std::vector<float> &train_embedding;
  std::size_t n_train_vertices;

  const std::vector<int> &nn_index;
  const std::vector<float> &nn_weights;
  std::size_t n_test_vertices;

  std::size_t ndim;
  std::size_t n_neighbors;
  std::vector<float> embedding;

  WeightedAverageWorker(const std::vector<float> &train_embedding,
                        std::size_t n_train_vertices,
                        const std::vector<int> &nn_index,
                        const std::vector<float> &nn_weights,
                        std::size_t n_test_vertices)
      : train_embedding(train_embedding), n_train_vertices(n_train_vertices),
        nn_index(nn_index), nn_weights(nn_weights),
        n_test_vertices(n_test_vertices),
        ndim(train_embedding.size() / n_train_vertices),
        n_neighbors(nn_index.size() / n_test_vertices),
        embedding(n_test_vertices * n_neighbors) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> sumc(ndim);
    for (std::size_t i = begin; i < end; i++) {
      std::fill(sumc.begin(), sumc.end(), 0.0);

      double sumw = 0.0;

      for (std::size_t j = 0; j < n_neighbors; j++) {
        std::size_t nbr = nn_index[i + j * n_test_vertices];
        float w = nn_weights[i + j * n_test_vertices];
        sumw += w;
        for (std::size_t k = 0; k < ndim; k++) {
          sumc[k] += train_embedding[nbr + k * n_train_vertices] * w;
        }
      }

      for (std::size_t k = 0; k < ndim; k++) {
        embedding[i + k * n_test_vertices] = sumc[k] / sumw;
      }
    }
  }
};

// Initialize embedding as a weighted average of nearest neighbors of each point
// train_embedding: n_train x dim matrix of final embedding coordinates
// nn_index: n_test x n_nbrs matrix of indexes of neighbors in X_train that are
//   nearest neighbors of X_test
// weights: n_test x n_nbrs weight matrix
// Returns the n_test x dim matrix of initialized coordinates.
// [[Rcpp::export]]
Rcpp::NumericMatrix init_transform_parallel(Rcpp::NumericMatrix train_embedding,
                                            Rcpp::IntegerMatrix nn_index,
                                            Rcpp::NumericMatrix nn_weights,
                                            std::size_t n_threads = 0,
                                            std::size_t grain_size = 1) {

  std::size_t n_train_vertices = train_embedding.nrow();
  std::size_t ndim = train_embedding.ncol();
  std::size_t n_test_vertices = nn_index.nrow();

  auto train_embeddingv = Rcpp::as<std::vector<float>>(train_embedding);
  auto nn_indexv = Rcpp::as<std::vector<int>>(nn_index);
  // Convert to zero-indexing
  for (int &i : nn_indexv) {
    --i;
  }
  auto nn_weightsv = Rcpp::as<std::vector<float>>(nn_weights);

  WeightedAverageWorker worker(train_embeddingv, n_train_vertices, nn_indexv,
                               nn_weightsv, n_test_vertices);

  if (n_threads > 0) {
    RcppPerpendicular::parallel_for(0, n_test_vertices, worker, n_threads,
                                    grain_size);
  } else {
    worker(0, n_test_vertices);
  }

  return Rcpp::NumericMatrix(n_test_vertices, ndim, worker.embedding.begin());
}
