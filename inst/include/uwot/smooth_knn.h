// BSD 2-Clause License
//
// Copyright 2020 James Melville
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// OF SUCH DAMAGE.

#ifndef UWOT_SMOOTH_KNN_H
#define UWOT_SMOOTH_KNN_H

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "RcppPerpendicular.h"
#include "uwot/matrix.h"

namespace uwot {

// Welford-style mean calculation
auto mean_average(std::vector<double> v) -> double {
  long double mean = 0.0;
  std::size_t n = v.size();

  for (std::size_t i = 0; i < n; ++i) {
    mean += (v[i] - mean) / (i + 1);
  }
  return static_cast<double>(mean);
}

struct SmoothKnnWorker {
  const std::vector<double> &nn_dist;

  std::size_t n_vertices;
  std::size_t n_neighbors;

  double target;
  std::size_t n_iter;
  double local_connectivity;
  double bandwidth;
  double tol;
  double min_k_dist_scale;
  double mean_distances;
  double double_max = (std::numeric_limits<double>::max)();

  std::vector<double> nn_weights;

  std::atomic_size_t n_search_fails;

  SmoothKnnWorker(const std::vector<double> &nn_dist, std::size_t n_vertices,
                  std::size_t n_iter, double local_connectivity,
                  double bandwidth, double tol, double min_k_dist_scale)
      : nn_dist(nn_dist), n_vertices(n_vertices),
        n_neighbors(nn_dist.size() / n_vertices),
        target(std::log2(n_neighbors)), n_iter(n_iter),
        local_connectivity(local_connectivity), bandwidth(bandwidth), tol(tol),
        min_k_dist_scale(min_k_dist_scale),
        mean_distances(mean_average(nn_dist)),
        nn_weights(n_vertices * n_neighbors), n_search_fails(0) {}

  void operator()(std::size_t begin, std::size_t end) {
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
      uwot::get_row(nn_dist, n_vertices, n_neighbors, i, ith_distances);
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
      for (std::size_t iter = 0; iter < n_iter; iter++) {
        double val = 0.0;
        // NB we iterate from 1, not 0: don't use the self-distance.
        for (std::size_t k = 1; k < n_neighbors; k++) {
          double dist = (std::max)(0.0, ith_distances[k] - rho);
          val += std::exp(-dist / sigma);
        }

        double adiff = std::abs(val - target);
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

      uwot::set_row(nn_weights, n_vertices, n_neighbors, i, res);
    }

    // Update global count of failures
    n_search_fails += n_window_search_fails;
  }
};
} // namespace uwot

#endif // UWOT_SMOOTH_KNN_H
