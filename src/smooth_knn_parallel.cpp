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
#include <numeric>
#include <limits>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
// // [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>

struct SmoothKnnWorker : public RcppParallel::Worker {
  const RcppParallel::RMatrix<double> nn_dist;
  const RcppParallel::RMatrix<int> nn_idx;
  const unsigned int n_vertices;
  const unsigned int n_neighbors;
  
  arma::umat locations;
  arma::vec values;
  
  const double target;
  const unsigned int n_iter;
  const double local_connectivity;
  const double bandwidth;
  const double tol;
  const double min_k_dist_scale;
  const double mean_distances;
  const double double_max = std::numeric_limits<double>::max();
  const bool self_nbr;
  
  // Progress progress;
  // tthread::mutex mutex;
  
  SmoothKnnWorker(const Rcpp::NumericMatrix& nn_dist, const Rcpp::IntegerMatrix&  nn_idx,
                  const unsigned int n_iter, const double local_connectivity,
                  const double bandwidth, const double tol, const double min_k_dist_scale,
                  const bool self_nbr
                    // , Progress& progress
  ) :
    nn_dist(nn_dist), nn_idx(nn_idx), n_vertices(nn_dist.nrow()), n_neighbors(nn_dist.ncol()),
    locations(2, n_vertices * n_neighbors), values(n_vertices * n_neighbors),
    target(std::log2(n_neighbors)),
    n_iter(n_iter), local_connectivity(local_connectivity), bandwidth(bandwidth),
    tol(tol), min_k_dist_scale(min_k_dist_scale),
    mean_distances(mean(nn_dist)),
    self_nbr(self_nbr)
    // , progress(progress)
  {  }
  
  
  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> non_zero_distances(n_neighbors);
    
    double sigma = 1.0;
    for (std::size_t i = begin; i < end; i++) {
      non_zero_distances.clear();
      double lo = 0.0;
      double hi = double_max;
      
      auto ith_distances = nn_dist.row(i);
      for (size_t k = 0; k < ith_distances.size(); k++) {
        if (ith_distances[k] > 0.0) {
          non_zero_distances.push_back(ith_distances[k]);
        }
      }
      
      // Find rho, the distance to the nearest neighbor (excluding zero distance neighbors)
      double rho = 0.0;
      if (non_zero_distances.size() >= local_connectivity) {
        int index = static_cast<int>(std::floor(local_connectivity));
        double interpolation = local_connectivity - index;
        if (index > 0) {
          rho = non_zero_distances[index - 1];
          if (interpolation >= tol) {
            rho += interpolation * (non_zero_distances[index] - non_zero_distances[index - 1]);
          }
        }
        else {
          rho = interpolation * non_zero_distances[0];
        }
      }
      else if (non_zero_distances.size() > 0) {
        rho = *std::max_element(non_zero_distances.begin(), non_zero_distances.end());
      }
      
      for (unsigned int iter = 0; iter < n_iter; iter++) {
        double val = 0.0;
        // NB we iterate from 1, not 0: don't use the self-distance.
        // Makes using Rcpp sugar sufficiently awkward so do the explicit loop
        for (unsigned int k = 1; k < n_neighbors; k++) {
          double dist = std::max(0.0, ith_distances[k] - rho);
          val += std::exp(-dist / sigma);
        }
        
        if (std::abs(val - target) < tol) {
          break;
        }
        
        if (val > target) {
          hi = sigma;
          sigma = 0.5 * (lo + hi);
        }
        else {
          lo = sigma;
          if (hi == double_max) {
            sigma *= 2;
          }
          else {
            sigma = 0.5 * (lo + hi);
          }
        }
      }
      
      if (rho > 0.0) {
        double mean = std::accumulate(ith_distances.begin(), ith_distances.end(), 0.0) / ith_distances.size();
        sigma = std::max(min_k_dist_scale * mean, sigma);
      }
      else {
        sigma = std::max(min_k_dist_scale * mean_distances, sigma);
      }
      
      double res[n_neighbors];
      for (size_t k = 0; k < n_neighbors; k++) {
        double rk = ith_distances[k] - rho;
        if (rk <= 0) {
          res[k] = 1.0;
        }
        else {
          res[k] = std::exp(-rk / (sigma * bandwidth));
        }
      }
      
      unsigned int loc = i * n_neighbors;
      // loc is incremented in the loop
      for (unsigned int k = 0; k < n_neighbors; k++, loc++) {
        unsigned int j = nn_idx(i, k) - 1;
        
        locations(0, loc) = j;
        locations(1, loc) = i;
        if (i != j || !self_nbr) {
          values(loc) = res[k];
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
arma::sp_mat smooth_knn_distances_parallel(const Rcpp::NumericMatrix& nn_dist, 
                                           const Rcpp::IntegerMatrix& nn_idx,
                                           const unsigned int n_iter,
                                           const double local_connectivity,
                                           const double bandwidth,
                                           const double tol,
                                           const double min_k_dist_scale,
                                           const size_t grain_size = 1,
                                           const bool self_nbr = true,
                                           const unsigned int n_reference_vertices = 0,
                                           const bool verbose = false) {
  const unsigned int n_vertices = nn_dist.nrow();
  const unsigned int n_ref_vertices = n_reference_vertices > 0 ? n_reference_vertices : n_vertices;
  
  // Progress progress(n_vertices, verbose);
  SmoothKnnWorker worker(nn_dist, nn_idx, n_iter, local_connectivity,
                         bandwidth, tol, min_k_dist_scale,
                         self_nbr
                           // , progress
  );
  
  RcppParallel::parallelFor(0, n_vertices, worker, grain_size);
  
  return arma::sp_mat(
    false, // add_values
    worker.locations,
    worker.values,
    n_ref_vertices, n_vertices
  );
}
