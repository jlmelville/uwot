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

namespace uwot {

// Welford-style mean calculation
auto mean_average(const std::vector<double> &v) -> double {
  long double mean = 0.0;
  const std::size_t n = v.size();

  for (std::size_t i = 0; i < n; ++i) {
    mean += (v[i] - mean) / (i + 1);
  }
  return static_cast<double>(mean);
}

// Find rho, the distance to the local_connectivity-th nearest neighbor
// (excluding zero distance neighbors)
auto find_rho(std::size_t i, const std::vector<double> &non_zero_distances,
              double local_connectivity, double tol) -> double {
  double rho = 0.0;
  std::size_t nnzero = non_zero_distances.size();
  if (nnzero >= local_connectivity) {
    int index = static_cast<int>(std::floor(local_connectivity));
    double interpolation = local_connectivity - index;
    if (index > 0) {
      rho = non_zero_distances[index - 1];
      if (interpolation >= tol) {
        rho += interpolation *
               (non_zero_distances[index] - non_zero_distances[index - 1]);
      }
    } else if (nnzero > 0) {
      rho = interpolation * non_zero_distances[0];
    }
  } else if (nnzero > 0) {
    rho =
        *std::max_element(non_zero_distances.begin(), non_zero_distances.end());
  }

  return rho;
}

// Find the normalization factor for the smoothed distances
auto find_sigma(const std::vector<double> &ith_distances,
                std::size_t n_neighbors, double target, double rho, double tol,
                std::size_t n_iter, double min_k_dist_scale,
                double mean_distances, std::size_t &n_window_search_fails)
    -> double {
  constexpr auto double_max = (std::numeric_limits<double>::max)();

  // best value seen is used only if binary search fails
  // NB there is already a safeguard against sigma getting too large
  // so this is less of a problem than with the perplexity search
  double sigma = 1.0;
  double sigma_best = sigma;
  double adiff_min = double_max;

  double lo = 0.0;
  double hi = double_max;

  bool converged = false;
  for (std::size_t iter = 0; iter < n_iter; iter++) {
    double val = 0.0;
    // NB we iterate from 1, not 0: don't use the self-distance.
    for (std::size_t k = 1; k < n_neighbors; k++) {
      auto rk = ith_distances[k] - rho;
      val += rk <= 0.0 ? 1.0 : std::exp(-rk / sigma);
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
    sigma = (std::max)(min_k_dist_scale * mean_average(ith_distances), sigma);
  } else {
    sigma = (std::max)(min_k_dist_scale * mean_distances, sigma);
  }

  return sigma;
}

void smooth_knn(std::size_t i, const std::vector<double> &nn_dist,
                std::size_t n_vertices, std::size_t n_neighbors, double target,
                double local_connectivity, double tol, std::size_t n_iter,
                double bandwidth, double min_k_dist_scale,
                double mean_distances, bool save_sigmas,
                std::vector<double> &nn_weights,
                std::vector<double> &non_zero_distances,
                std::vector<double> &sigmas,
                std::size_t &n_window_search_fails) {
  non_zero_distances.clear();

  std::vector<double> ith_distances(n_neighbors);
  for (std::size_t j = 0; j < n_neighbors; j++) {
    auto ith_distance = nn_dist[i + j * n_vertices];
    ith_distances[j] = ith_distance;
    if (ith_distance > 0.0) {
      non_zero_distances.push_back(ith_distance);
    }
  }

  auto rho = find_rho(i, non_zero_distances, local_connectivity, tol);
  auto sigma =
      find_sigma(ith_distances, n_neighbors, target, rho, tol, n_iter,
                 min_k_dist_scale, mean_distances, n_window_search_fails);

  auto sigma_b = sigma * bandwidth;
  for (std::size_t k = 0; k < n_neighbors; k++) {
    auto rk = ith_distances[k] - rho;
    nn_weights[i + k * n_vertices] = rk <= 0.0 ? 1.0 : std::exp(-rk / sigma_b);
  }

  if (save_sigmas) {
    sigmas[i] = sigma;
  }
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

  std::vector<double> nn_weights;

  bool save_sigmas;
  std::vector<double> sigmas;

  std::atomic_size_t n_search_fails;

  SmoothKnnWorker(const std::vector<double> &nn_dist, std::size_t n_vertices,
                  std::size_t n_iter, double local_connectivity,
                  double bandwidth, double tol, double min_k_dist_scale,
                  bool save_sigmas, double target)
      : nn_dist(nn_dist), n_vertices(n_vertices),
        n_neighbors(nn_dist.size() / n_vertices), target(target),
        n_iter(n_iter), local_connectivity(local_connectivity),
        bandwidth(bandwidth), tol(tol), min_k_dist_scale(min_k_dist_scale),
        mean_distances(mean_average(nn_dist)),
        nn_weights(n_vertices * n_neighbors), save_sigmas(save_sigmas),
        sigmas(save_sigmas ? n_vertices : 0), n_search_fails(0) {}

  void operator()(std::size_t begin, std::size_t end) {
    // number of binary search failures in this window
    std::size_t n_window_search_fails = 0;
    std::vector<double> non_zero_distances;
    non_zero_distances.reserve(n_neighbors);

    for (std::size_t i = begin; i < end; i++) {
      smooth_knn(i, nn_dist, n_vertices, n_neighbors, target,
                 local_connectivity, tol, n_iter, bandwidth, min_k_dist_scale,
                 mean_distances, save_sigmas, nn_weights, non_zero_distances,
                 sigmas, n_window_search_fails);
    }

    // Update global count of failures
    n_search_fails += n_window_search_fails;
  }
};
} // namespace uwot

#endif // UWOT_SMOOTH_KNN_H
