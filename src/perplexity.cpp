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

#include <limits>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
// // [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>

struct PerplexityWorker : public RcppParallel::Worker {
  const RcppParallel::RMatrix<double> nn_dist;
  const RcppParallel::RMatrix<int> nn_idx;
  const unsigned int n_vertices;
  const unsigned int n_neighbors;

  arma::umat locations;
  arma::vec values;

  const double target;
  const unsigned int n_iter;
  const double tol;
  const double double_max = std::numeric_limits<double>::max();

  // Progress progress;
  // tthread::mutex mutex;

  PerplexityWorker(const Rcpp::NumericMatrix& nn_dist,
                   const Rcpp::IntegerMatrix&  nn_idx,
                   const double perplexity, const unsigned int n_iter,
                   const double tol
                     // , Progress& progress
  ) :
    nn_dist(nn_dist), nn_idx(nn_idx), n_vertices(nn_dist.nrow()), n_neighbors(nn_dist.ncol()),
    locations(2, n_vertices * n_neighbors), values(n_vertices * n_neighbors),
    target(std::log(perplexity)), n_iter(n_iter), tol(tol)
    // , progress(progress)
  {  }

  void operator()(std::size_t begin, std::size_t end) {
    double d2[n_neighbors - 1];

    // first vertex: guess initial beta as 1.0
    // subsequent vertices, guess the last optimized beta
    double beta = 1.0;

    for (std::size_t i = begin; i < end; i++) {
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

      unsigned int loc = i * n_neighbors;
      // loc is incremented in the loop
      for (unsigned int k = 0; k < n_neighbors; k++, loc++) {
        unsigned int j = nn_idx(i, k) - 1;

        locations(0, loc) = i;
        locations(1, loc) = j;
        if (i != j) {
          values(loc) = d2[k - 1] / Z;
        }
        else {
          values(loc) = 0.0;
        }
      }

      // {
      //   tthread::lock_guard<tthread::mutex> guard(mutex);
      //   progress.increment();
      //   if (Progress::check_abort()) {
      //     return;
      //   }
      // }
    }
  }
};

// [[Rcpp::export]]
arma::sp_mat calc_row_probabilities_parallel(const Rcpp::NumericMatrix& nn_dist, const Rcpp::IntegerMatrix& nn_idx,
                                             const double perplexity,
                                             const unsigned int n_iter = 200,
                                             const double tol = 1e-5,
                                             const bool parallelize = true,
                                             const std::size_t grain_size = 1,
                                             const bool verbose = false) {
  const unsigned int n_vertices = nn_dist.nrow();
  // Progress progress(n_vertices, verbose);
  PerplexityWorker worker(nn_dist, nn_idx, perplexity, n_iter, tol
                            // , progress
  );

  if (parallelize) {
    RcppParallel::parallelFor(0, n_vertices, worker, grain_size);
  }
  else {
    worker(0, n_vertices);
  }

  return arma::sp_mat(
    false, // add_values
    worker.locations,
    worker.values,
    n_vertices, n_vertices
  );
}
