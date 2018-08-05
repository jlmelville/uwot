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
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

// [[Rcpp::export]]
Rcpp::NumericMatrix smooth_knn_distances_cpp(
                                      const Rcpp::NumericMatrix& nn_dist,
                                      const Rcpp::IntegerMatrix& nn_idx,
                                      const unsigned int n_iter = 64,
                                      const double local_connectivity = 1.0,
                                      const double bandwidth = 1.0,
                                      const double tol = 1e-5,
                                      const double min_k_dist_scale = 1e-3,
                                      const bool verbose = false) {
  const unsigned int n_vertices = nn_dist.nrow();
  const unsigned int n_neighbors = nn_dist.ncol();

  // bandwidth is not used here on purpose
  const double target = log2(n_neighbors);
  const double double_max = std::numeric_limits<double>::max();

  double mean_distances = mean(nn_dist);

  Rcpp::NumericMatrix values(n_vertices, n_neighbors);
  Progress progress(n_vertices, verbose);

  // Initial guess for i = 0
  // On subsequent iterations we shall start from the previous optimized
  // value: seems to save around 20% on total number of iterations
  double sigma = 1.0;
  for (unsigned int i = 0; i < n_vertices; i++) {
    double lo = 0.0;
    double hi = double_max;

    Rcpp::NumericVector ith_distances = nn_dist.row(i);
    Rcpp::NumericVector non_zero_distances = ith_distances[ith_distances > 0.0];

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
      rho = Rcpp::max(non_zero_distances);
    }

    for (unsigned int iter = 0; iter < n_iter; iter++) {
      double val = 0.0;
      // NB we iterate from 1, not 0: don't use the self-distance.
      // Makes using Rcpp sugar sufficiently awkward so do the explicit loop
      for (unsigned int j = 1; j < n_neighbors; j++) {
        double dist = std::max(0.0, ith_distances[j] - rho);
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
      sigma = std::max(min_k_dist_scale * mean(ith_distances), sigma);
    }
    else {
      sigma = std::max(min_k_dist_scale * mean_distances, sigma);
    }

    Rcpp::NumericVector res = Rcpp::exp(-(ith_distances - rho) / (sigma * bandwidth));
    res[ith_distances - rho <= 0] = 1.0;

    for (unsigned int k = 0; k < n_neighbors; k++) {
      values.at(i, k) = res(k);
    }

    if (progress.check_abort()) {
      Rcpp::stop("Progress aborted by user");
    }
    progress.increment();
  }

  return values;
}
