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

#include <Rcpp.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]




struct AverageWorker : public RcppParallel::Worker {

  const RcppParallel::RMatrix<double> train_embedding;
  const RcppParallel::RMatrix<int> nn_index;
  RcppParallel::RMatrix<double> embedding;
  const size_t nc;
  const size_t nnbrs;
  const double one_over_n;

  AverageWorker(Rcpp::NumericMatrix train_embedding, Rcpp::IntegerMatrix nn_index,
                Rcpp::NumericMatrix embedding
  ) :
    train_embedding(train_embedding), nn_index(nn_index),
    embedding(embedding),
    nc(train_embedding.ncol()), nnbrs(nn_index.ncol()),
    one_over_n(1.0 / nnbrs)
  {  }

  void operator()(std::size_t begin, std::size_t end) {
    double sumc[nc];
    for (std::size_t i = begin; i < end; i++) {
      for (size_t k = 0; k < nc; k++) {
        sumc[k] = 0.0;
      }

      for (size_t j = 0; j < nnbrs; j++) {
        auto nbr = nn_index(i, j) - 1;
        for (size_t k = 0; k < nc; k++) {
          sumc[k] += train_embedding(nbr, k);
        }
      }

      for (size_t k = 0; k < nc; k++) {
        embedding(i, k) = sumc[k] * one_over_n;
      }
    }
  }
};


// [[Rcpp::export]]
Rcpp::NumericMatrix init_transform_av_parallel(Rcpp::NumericMatrix train_embedding,
                                               Rcpp::IntegerMatrix nn_index,
                                               bool parallelize = true,
                                               const size_t grain_size = 1) {
  Rcpp::NumericMatrix embedding(nn_index.nrow(), train_embedding.ncol());

  AverageWorker worker(train_embedding, nn_index, embedding);

  if (parallelize) {
    RcppParallel::parallelFor(0, nn_index.nrow(), worker, grain_size);
  }
  else {
    worker(0, nn_index.nrow());
  }

  return embedding;
}

struct WeightedAverageWorker : public RcppParallel::Worker {

  const RcppParallel::RMatrix<double> train_embedding;
  const RcppParallel::RMatrix<int> nn_index;
  const RcppParallel::RMatrix<double> nn_weights;
  RcppParallel::RMatrix<double> embedding;
  const size_t nc;
  const size_t nnbrs;

  WeightedAverageWorker(Rcpp::NumericMatrix train_embedding, Rcpp::IntegerMatrix nn_index,
                        const Rcpp::NumericMatrix& nn_weights, Rcpp::NumericMatrix embedding
  ) :
    train_embedding(train_embedding), nn_index(nn_index), nn_weights(nn_weights),
    embedding(embedding),
    nc(train_embedding.ncol()), nnbrs(nn_index.ncol())
  {  }

  void operator()(std::size_t begin, std::size_t end) {
    double sumc[nc];
    for (std::size_t i = begin; i < end; i++) {
      for (size_t k = 0; k < nc; k++) {
        sumc[k] = 0.0;
      }
      double sumw = 0.0;

      for (size_t j = 0; j < nnbrs; j++) {
        auto nbr = nn_index(i, j) - 1;
        double w = nn_weights(i, j);
        sumw += w;
        for (size_t k = 0; k < nc; k++) {
          sumc[k] += train_embedding(nbr, k) * w;
        }
      }

      for (size_t k = 0; k < nc; k++) {
        embedding(i, k) = sumc[k] / sumw;
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
                                            const size_t grain_size = 1,
                                            bool parallelize = true) {
  Rcpp::NumericMatrix embedding(nn_index.nrow(), train_embedding.ncol());

  WeightedAverageWorker worker(train_embedding, nn_index, nn_weights, embedding);

  if (parallelize) {
    RcppParallel::parallelFor(0, nn_index.nrow(), worker, grain_size);
  }
  else {
    worker(0, nn_index.nrow());
  }

  return embedding;
}
