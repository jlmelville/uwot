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

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Initialize embedding as a weighted average of nearest neighbors of each point
// train_embedding: n_train x dim matrix of final embedding coordinates
// nn_index: n_test x n_nbrs matrix of indexes of neighbors in X_train that are 
//   nearest neighbors of X_test
// weights: n_train x n_test sparse weight matrix, with n_nbrs weights per columns
//  corresponding to the indexes in nn_index. Some of these may be zero. There
//  should not be any non-zero weights in the other locations of the matrix.
//  TODO: this input type is used because it's the output of smooth_knn_distances,
//  but it would make more sense to use a dense output matrix like nn_index here.
//  The sparse matrix only needs forming prior to symmetrization or forming the
//  graph.
// Returns the n_test x dim matrix of initialized coordinates.
// [[Rcpp::export]]
Rcpp::NumericMatrix init_transform_cpp(Rcpp::NumericMatrix train_embedding, 
                                       Rcpp::IntegerMatrix nn_index,
                                       arma::sp_mat weights) {
  auto nr = nn_index.nrow();
  auto nc = train_embedding.ncol();
  auto nnbrs = nn_index.ncol();
  
  Rcpp::NumericMatrix embedding(nr, nc);

  double sumc[nc];

  for (auto i = 0; i < nr; i++) {
    for (auto k = 0; k < nc; k++) {
      sumc[k] = 0.0;
    }
    double sumw = 0.0;
    
    for (auto j = 0; j < nnbrs; j++) {
      auto nbr = nn_index.at(i, j) - 1;
      double w = weights(nbr, i);
      sumw += w;
      for (auto k = 0; k < nc; k++) {
        sumc[k] += train_embedding.at(nbr, k) * w;
      }
    }
    
    for (auto k = 0; k < nc; k++) {
      embedding.at(i, k) = sumc[k] / sumw;
    }
  }
  
  return embedding;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix init_transform_av_cpp(Rcpp::NumericMatrix train_embedding, 
                                          Rcpp::IntegerMatrix nn_index) {
  auto nr = nn_index.nrow();
  auto nc = train_embedding.ncol();
  auto nnbrs = nn_index.ncol();
  
  Rcpp::NumericMatrix embedding(nr, nc);
  
  double sumc[nc];
  double one_over_n = 1.0 / nnbrs;

  for (auto i = 0; i < nr; i++) {
    for (auto k = 0; k < nc; k++) {
      sumc[k] = 0.0;
    }

    for (auto j = 0; j < nnbrs; j++) {
      auto nbr = nn_index.at(i, j) - 1;
      for (auto k = 0; k < nc; k++) {
        sumc[k] += train_embedding.at(nbr, k);
      }
    }
    
    for (auto k = 0; k < nc; k++) {
      embedding.at(i, k) = sumc[k] * one_over_n;
    }
  }
  
  return embedding;
}