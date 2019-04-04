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
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct PerplexityWorker : public RcppParallel::Worker {
  RcppParallel::RMatrix<double> res;
  const RcppParallel::RMatrix<double> nn_dist;
  const RcppParallel::RMatrix<int> nn_idx;
  const unsigned int n_vertices;
  const unsigned int n_neighbors;

  const double target;
  const unsigned int n_iter;
  const double tol;
  const double double_max = std::numeric_limits<double>::max();
  
  PerplexityWorker(
    Rcpp::NumericMatrix res,
    const Rcpp::NumericMatrix nn_dist,
    const Rcpp::IntegerMatrix nn_idx,
    const double perplexity, 
    const unsigned int n_iter,
    const double tol
  ) :
    res(res), nn_dist(nn_dist), nn_idx(nn_idx), n_vertices(nn_dist.nrow()), 
    n_neighbors(nn_dist.ncol()),target(std::log(perplexity)), n_iter(n_iter), 
    tol(tol)
  {  }
  
  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> d2(n_neighbors - 1, 0.0);
    
    for (std::size_t i = begin; i < end; i++) {
      double beta = 1.0;
      double lo = 0.0;
      double hi = double_max;
      
      for (unsigned int k = 1; k < n_neighbors; k++) {
        d2[k - 1] = nn_dist(i, k) * nn_dist(i, k);
      }
      
      for (unsigned int iter = 0; iter < n_iter; iter++) {
        double Z = 0.0;
        double H = 0.0;
        double sum_D2_W = 0.0;
        for (unsigned int k = 0; k < n_neighbors - 1; k++) {
          double W = std::exp(-d2[k] * beta);
          Z += W;
          sum_D2_W += d2[k] * W;
        }
        if (Z > 0) {
          H = std::log(Z) + beta * sum_D2_W / Z;
        }

        if (std::abs(H - target) < tol) {
          break;
        }
        
        if (H < target) {
          hi = beta;
          beta = 0.5 * (lo + hi);
        }
        else {
          lo = beta;
          if (hi == double_max) {
            beta *= 2.0;
          }
          else {
            beta = 0.5 * (lo + hi);
          }
        }
      }
      
      double Z = 0.0;
      for (unsigned int k = 0; k < n_neighbors - 1; k++) {
        double W = std::exp(-d2[k] * beta);
        Z += W;
        // no longer need d2 at this point, store final W there
        d2[k] = W;
      }
      
      // This will index over d2, skipping when i == j
      std::size_t widx = 0;
      for (unsigned int k = 0; k < n_neighbors; k++) {
        unsigned int j = nn_idx(i, k) - 1;
        if (i != j) {
          res(i, k) = d2[widx] / Z;
          ++widx;
        }
        else {
          res(i, k) = 0.0;
        }
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix calc_row_probabilities_parallel(
    const Rcpp::NumericMatrix nn_dist, 
    const Rcpp::IntegerMatrix nn_idx,
    const double perplexity,
    const unsigned int n_iter = 200,
    const double tol = 1e-5,
    const bool parallelize = true,
    const std::size_t grain_size = 1,
    const bool verbose = false) {
  Rcpp::NumericMatrix res = Rcpp::NumericMatrix(nn_dist.nrow(), nn_dist.ncol());
  const unsigned int n_vertices = nn_dist.nrow();
  PerplexityWorker worker(res, nn_dist, nn_idx, perplexity, n_iter, tol);
  
  if (parallelize) {
    RcppParallel::parallelFor(0, n_vertices, worker, grain_size);
  }
  else {
    worker(0, n_vertices);
  }
  
  return res;
}
