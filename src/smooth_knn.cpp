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
#include <limits>
#include <mutex>
#include <numeric>
#include <vector>

#include "RcppPerpendicular.h"
#include <Rcpp.h>

#include "matrix.h"

// As in R's internals and Rcpp, use Kahan Summation to compensate for
// numerical error
// https://stackoverflow.com/q/17866149/4096483
auto mean_average(std::vector<double> v) -> double {
  std::size_t n = v.size();
  long double s = std::accumulate(v.begin(), v.end(), 0.0L);
  s /= n;

  long double t = 0.0;
  for (std::size_t i = 0; i < n; i++) {
    t += v[i] - s;
  }
  s += t / n;

  return static_cast<double>(s);
}

struct SmoothKnnWorker : public RcppPerpendicular::Worker {
  const std::vector<double> &nn_dist;

  const unsigned int n_vertices;
  const unsigned int n_neighbors;

  const double target;
  const unsigned int n_iter;
  const double local_connectivity;
  const double bandwidth;
  const double tol;
  const double min_k_dist_scale;
  const double mean_distances;
  const double double_max = (std::numeric_limits<double>::max)();

  std::vector<double> nn_weights;

  std::mutex mutex;
  std::size_t n_search_fails;

  SmoothKnnWorker(const std::vector<double> &nn_dist,
                  const unsigned int n_vertices, const unsigned int n_iter,
                  const double local_connectivity, const double bandwidth,
                  const double tol, const double min_k_dist_scale)
      : nn_dist(nn_dist), n_vertices(n_vertices),
        n_neighbors(nn_dist.size() / n_vertices),
        target(std::log2(n_neighbors)), n_iter(n_iter),
        local_connectivity(local_connectivity), bandwidth(bandwidth), tol(tol),
        min_k_dist_scale(min_k_dist_scale),
        mean_distances(mean_average(nn_dist)),
        nn_weights(n_vertices * n_neighbors), n_search_fails(0) {}

  void operator()(std::size_t begin, std::size_t end) override {
    // number of binary search failures in this window
    std::size_t n_window_search_fails = 0;
    std::vector<double> non_zero_distances;
    non_zero_distances.reserve(n_neighbors);

    for (std::size_t i = begin; i < end; i++) {
      double sigma = 1.0;
      non_zero_distances.clear();
      double lo = 0.0;
      double hi = double_max;
      // best value seen is used only if binary search fails
      // NB there is already a safeguard against sigma getting too large
      // so this is less of a problem than with the perplexity search
      double sigma_best = sigma;
      double adiff_min = double_max;

      std::vector<double> ith_distances(n_neighbors);
      get_row(nn_dist, n_vertices, n_neighbors, i, ith_distances);
      for (double ith_distance : ith_distances) {
        if (ith_distance > 0.0) {
          non_zero_distances.push_back(ith_distance);
        }
      }

      // Find rho, the distance to the nearest neighbor (excluding zero distance
      // neighbors)
      double rho = 0.0;
      if (non_zero_distances.size() >= local_connectivity) {
        int index = static_cast<int>(std::floor(local_connectivity));
        double interpolation = local_connectivity - index;
        if (index > 0) {
          rho = non_zero_distances[index - 1];
          if (interpolation >= tol) {
            rho += interpolation *
                   (non_zero_distances[index] - non_zero_distances[index - 1]);
          }
        } else if (non_zero_distances.size() > 0) {
          rho = interpolation * non_zero_distances[0];
        }
      } else if (non_zero_distances.size() > 0) {
        rho = *std::max_element(non_zero_distances.begin(),
                                non_zero_distances.end());
      }

      bool converged = false;
      for (unsigned int iter = 0; iter < n_iter; iter++) {
        double val = 0.0;
        // NB we iterate from 1, not 0: don't use the self-distance.
        for (unsigned int k = 1; k < n_neighbors; k++) {
          double dist = (std::max)(0.0, ith_distances[k] - rho);
          val += std::exp(-dist / sigma);
        }

        const double adiff = std::abs(val - target);
        if (adiff < tol) {
          converged = true;
          break;
        }

        // store best sigma in case binary search fails (usually in the presence
        // of multiple degenerate distances)
        if (adiff < adiff_min) {
          adiff_min = adiff;
          sigma_best = sigma;
        }

        if (val > target) {
          hi = sigma;
          sigma = 0.5 * (lo + hi);
        } else {
          lo = sigma;
          if (hi == double_max) {
            sigma *= 2;
          } else {
            sigma = 0.5 * (lo + hi);
          }
        }
      }
      if (!converged) {
        ++n_window_search_fails;
        sigma = sigma_best;
      }

      if (rho > 0.0) {
        sigma =
            (std::max)(min_k_dist_scale * mean_average(ith_distances), sigma);
      } else {
        sigma = (std::max)(min_k_dist_scale * mean_distances, sigma);
      }

      std::vector<double> res(n_neighbors, 0.0);
      for (std::size_t k = 0; k < n_neighbors; k++) {
        double rk = ith_distances[k] - rho;
        if (rk <= 0) {
          res[k] = 1.0;
        } else {
          res[k] = std::exp(-rk / (sigma * bandwidth));
        }
      }

      set_row(nn_weights, n_vertices, n_neighbors, i, res);
    }
    // Update global count of failures
    {
      std::lock_guard<std::mutex> guard(mutex);
      n_search_fails += n_window_search_fails;
    }
  }
};

// [[Rcpp::export]]
Rcpp::List smooth_knn_distances_parallel(
    const Rcpp::NumericMatrix &nn_dist, const unsigned int n_iter = 64,
    const double local_connectivity = 1.0, const double bandwidth = 1.0,
    const double tol = 1e-5, const double min_k_dist_scale = 1e-3,
    const bool parallelize = true, const std::size_t grain_size = 1,
    const bool verbose = false) {

  const unsigned int n_vertices = nn_dist.nrow();
  const unsigned int n_neighbors = nn_dist.ncol();

  auto nn_distv = Rcpp::as<std::vector<double>>(nn_dist);
  SmoothKnnWorker worker(nn_distv, n_vertices, n_iter, local_connectivity,
                         bandwidth, tol, min_k_dist_scale);

  if (parallelize) {
    RcppPerpendicular::parallelFor(0, n_vertices, worker, grain_size);
  } else {
    worker(0, n_vertices);
  }

  return Rcpp::List::create(Rcpp::Named("matrix") =
                                Rcpp::NumericMatrix(n_vertices, n_neighbors,
                                                    worker.nn_weights.begin()),
                            Rcpp::Named("n_failures") = worker.n_search_fails);
}
