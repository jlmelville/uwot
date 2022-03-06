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

#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "RcppPerpendicular.h"

namespace uwot {

// Welford-style mean calculation
auto mean_average(const std::vector<double> &v, std::size_t begin,
                  std::size_t end) -> double {
  long double mean = 0.0;
  auto b1 = 1 - begin;
  for (auto i = begin; i < end; ++i) {
    mean += (v[i] - mean) / (i + b1);
  }
  return static_cast<double>(mean);
}

auto mean_average(const std::vector<double> &v) -> double {
  return mean_average(v, 0, v.size());
}

// nn_dist is sorted non-decreasing nearest neighbor distances
// nzero_begin points to the index of the first non-zero distance
// nzero_end points to one past the index of the last non-zero distance
// (n_neighbors + 1)
auto find_rho(const std::vector<double> &nn_dist, std::size_t nzero_begin,
              std::size_t nzero_end, double local_connectivity, double tol)
    -> double {
  double rho = 0.0;
  auto nnzero = nzero_end - nzero_begin;
  if (nnzero >= local_connectivity) {
    auto index = static_cast<int>(std::floor(local_connectivity));
    double interpolation = local_connectivity - index;
    if (index > 0) {
      rho = nn_dist[nzero_begin + index - 1];
      if (interpolation >= tol) {
        // rho = (1 - interp) * rho + interp * d_{i+1}
        rho += interpolation * (nn_dist[nzero_begin + index] - rho);
      }
    } else if (nnzero > 0) {
      // typical code-path: rho is the smallest non-zero distance
      rho = interpolation * nn_dist[nzero_begin];
    }
  } else if (nnzero > 0) {
    // not enough non-zero distances, return the largest non-zero distance
    rho = nn_dist[nzero_end - 1];
  }

  return rho;
}

// Find the normalization factor for the smoothed distances
auto find_sigma(const std::vector<double> &nn_dist, std::size_t i_begin,
                std::size_t i_end, double target, double rho, double tol,
                std::size_t n_iter, std::size_t &n_window_search_fails)
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
    // NB i_begin should point to the first non-self neighbor
    for (auto j = i_begin; j < i_end; j++) {
      auto rj = nn_dist[j] - rho;
      val += rj <= 0.0 ? 1.0 : std::exp(-rj / sigma);
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

  return sigma;
}

// NB nn_dist must be in sorted non-decreasing order
void smooth_knn(std::size_t i, const std::vector<double> &nn_dist,
                const std::vector<std::size_t> &nn_ptr, bool skip_first,
                const std::vector<double> &target, double local_connectivity,
                double tol, std::size_t n_iter, double min_k_dist_scale,
                double mean_distances, bool save_sigmas,
                std::vector<double> &nn_weights, std::vector<double> &sigmas,
                std::vector<double> &rhos, std::size_t &n_window_search_fails) {

  // i_begin points to start of ith distances
  // i_end points to one past end of ith distances
  auto i_begin = 0;
  auto i_end = 0;
  auto n_neighbors = 0;
  // Space optimization for kNN (typical case): store the number of neighbors
  // as the only entry in nn_ptr
  if (nn_ptr.size() == 1) {
    n_neighbors = nn_ptr[0];
    i_begin = n_neighbors * i;
    i_end = i_begin + n_neighbors;
  } else {
    i_begin = nn_ptr[i];
    i_end = nn_ptr[i + 1];
    n_neighbors = i_end - i_begin;
  }

  // nzero_begin points to start of ith non-zero distances
  auto nzero_begin = i_end;
  for (auto j = i_begin; j < i_end; j++) {
    if (nn_dist[j] > 0.0) {
      nzero_begin = j;
      break;
    }
  }

  auto rho = find_rho(nn_dist, nzero_begin, i_end, local_connectivity, tol);
  double targeti = target.size() == 1 ? target[0] : target[i];
  // in case where self-distance (0) is passed as the nearest neighbor, skip
  // first item in neighbors when calculating sigma
  auto sigma = find_sigma(nn_dist, i_begin + (skip_first ? 1 : 0), i_end,
                          targeti, rho, tol, n_iter, n_window_search_fails);
  // safeguard sigma
  if (rho > 0.0) {
    sigma = (std::max)(min_k_dist_scale * mean_average(nn_dist, i_begin, i_end),
                       sigma);
  } else {
    sigma = (std::max)(min_k_dist_scale * mean_distances, sigma);
  }

  // create the final membership strengths
  for (auto j = i_begin; j < i_end; j++) {
    auto rj = nn_dist[j] - rho;
    nn_weights[j] = rj <= 0.0 ? 1.0 : std::exp(-rj / sigma);
  }

  if (save_sigmas) {
    sigmas[i] = sigma;
    rhos[i] = rho;
  }
}

void smooth_knn(std::size_t begin, std::size_t end,
                const std::vector<double> &nn_dist,
                const std::vector<std::size_t> &nn_ptr, bool skip_first,
                const std::vector<double> &target, double local_connectivity,
                double tol, std::size_t n_iter, double min_k_dist_scale,
                double mean_distances, bool save_sigmas,
                std::vector<double> &nn_weights, std::vector<double> &sigmas,
                std::vector<double> &rhos, std::atomic_size_t &n_search_fails) {
  // number of binary search failures in this window
  std::size_t n_window_search_fails = 0;

  for (std::size_t i = begin; i < end; i++) {
    smooth_knn(i, nn_dist, nn_ptr, skip_first, target, local_connectivity, tol,
               n_iter, min_k_dist_scale, mean_distances, save_sigmas,
               nn_weights, sigmas, rhos, n_window_search_fails);
  }

  // Update global count of failures
  n_search_fails += n_window_search_fails;
}

} // namespace uwot

#endif // UWOT_SMOOTH_KNN_H
