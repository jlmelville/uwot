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
  
  tthread::mutex mutex;
  std::size_t n_search_fails;
  
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
    tol(tol), n_search_fails(0)
  {  }
  
  void operator()(std::size_t begin, std::size_t end) {
    // number of binary search failures in this window
    std::size_t n_window_search_fails = 0;
    std::vector<double> d2(n_neighbors - 1, 0.0);
    
    for (std::size_t i = begin; i < end; i++) {
      double beta = 1.0;
      double lo = 0.0;
      double hi = double_max;

      // best value seen is used only if binary search fails
      // (usually only happens if there are multiple degenerate distances)
      double beta_best = beta;
      double adiff_min = double_max;
      bool converged = false;
      
      // calculate squared distances and remember the minimum
      double dmin = double_max;
      double dtmp;
      for (unsigned int k = 1; k < n_neighbors; k++) {
        dtmp = nn_dist(i, k) * nn_dist(i, k);
        d2[k - 1] = dtmp;
        if (dtmp < dmin) {
          dmin = dtmp;
        }
      }
      // shift distances by minimum: this implements the log-sum-exp trick
      // D2, W and Z are their shifted versions
      // but P (and hence Shannon entropy) is unchanged
      for (unsigned int k = 1; k < n_neighbors; k++) {
        d2[k - 1] -= dmin;
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

        const double adiff = std::abs(H - target);
        if (adiff < tol) {
          converged = true;
          break;
        }
        
        // store best beta in case binary search fails
        if (adiff < adiff_min) {
          adiff_min = adiff;
          beta_best = beta;
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
      if (!converged) {
        ++n_window_search_fails;
        beta = beta_best;
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
    
    // Update global count of failures
    {
      tthread::lock_guard<tthread::mutex> guard(mutex);
      n_search_fails += n_window_search_fails;
    }
  }
};

// [[Rcpp::export]]
Rcpp::List calc_row_probabilities_parallel(
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

  return Rcpp::List::create(
    Rcpp::Named("matrix") = res,
    Rcpp::Named("n_failures") = worker.n_search_fails
  );
}
