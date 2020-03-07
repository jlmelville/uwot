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

#ifndef UWOT_PERPLEXITY_H
#define UWOT_PERPLEXITY_H

#include <atomic>
#include <limits>
#include <vector>

namespace uwot {

struct PerplexityWorker {
  const std::vector<double> &nn_dist;
  const std::vector<int> &nn_idx;
  std::size_t n_vertices;
  std::size_t n_neighbors;

  double target;
  std::size_t n_iter;
  double tol;
  double double_max = (std::numeric_limits<double>::max)();

  std::vector<double> res;
  std::atomic_size_t n_search_fails;

  PerplexityWorker(const std::vector<double> &nn_dist,
                   const std::vector<int> &nn_idx, std::size_t n_vertices,
                   double perplexity, std::size_t n_iter, double tol)
      : nn_dist(nn_dist), nn_idx(nn_idx), n_vertices(n_vertices),
        n_neighbors(nn_dist.size() / n_vertices), target(std::log(perplexity)),
        n_iter(n_iter), tol(tol), res(n_vertices * n_neighbors),
        n_search_fails(0) {}

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
      for (std::size_t k = 1; k < n_neighbors; k++) {
        dtmp = nn_dist[i + k * n_vertices] * nn_dist[i + k * n_vertices];
        d2[k - 1] = dtmp;
        if (dtmp < dmin) {
          dmin = dtmp;
        }
      }
      // shift distances by minimum: this implements the log-sum-exp trick
      // D2, W and Z are their shifted versions
      // but P (and hence Shannon entropy) is unchanged
      for (std::size_t k = 1; k < n_neighbors; k++) {
        d2[k - 1] -= dmin;
      }

      for (std::size_t iter = 0; iter < n_iter; iter++) {
        double Z = 0.0;
        double H = 0.0;
        double sum_D2_W = 0.0;
        for (std::size_t k = 0; k < n_neighbors - 1; k++) {
          double W = std::exp(-d2[k] * beta);
          Z += W;
          sum_D2_W += d2[k] * W;
        }
        if (Z > 0) {
          H = std::log(Z) + beta * sum_D2_W / Z;
        }

        double adiff = std::abs(H - target);
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
        } else {
          lo = beta;
          if (hi == double_max) {
            beta *= 2.0;
          } else {
            beta = 0.5 * (lo + hi);
          }
        }
      }
      if (!converged) {
        ++n_window_search_fails;
        beta = beta_best;
      }

      double Z = 0.0;
      for (std::size_t k = 0; k < n_neighbors - 1; k++) {
        double W = std::exp(-d2[k] * beta);
        Z += W;
        // no longer need d2 at this point, store final W there
        d2[k] = W;
      }

      // This will index over d2, skipping when i == j
      std::size_t widx = 0;
      for (std::size_t k = 0; k < n_neighbors; k++) {
        std::size_t j = nn_idx[i + k * n_vertices] - 1;
        if (i != j) {
          res[i + k * n_vertices] = d2[widx] / Z;
          ++widx;
        }
      }
    }

    // Update global count of failures
    n_search_fails += n_window_search_fails;
  }
};

} // namespace uwot

#endif // UWOT_PERPLEXITY_H
