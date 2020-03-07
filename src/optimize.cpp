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

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "RcppPerpendicular.h"
#include "uwot/gradient.h"
#include "uwot/optimize.h"
#include "uwot/sampler.h"

#include "rng.h"

using namespace Rcpp;

template <typename T, bool DoMove = true, typename RandFactory = pcg_factory>
auto optimize_layout(const T &gradient, std::vector<float> &head_embedding,
                     std::vector<float> &tail_embedding,
                     const std::vector<unsigned int> &positive_head,
                     const std::vector<unsigned int> &positive_tail,
                     unsigned int n_epochs, unsigned int n_vertices,
                     const std::vector<float> &epochs_per_sample,
                     float initial_alpha, float negative_sample_rate,
                     std::size_t n_threads = 0, std::size_t grain_size = 1,
                     bool verbose = false) -> std::vector<float> {
  uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);

  uwot::SgdWorker<T, DoMove, RandFactory> worker(
      gradient, positive_head, positive_tail, sampler, head_embedding,
      tail_embedding, head_embedding.size() / n_vertices);

  Progress progress(n_epochs, verbose);
  const auto n_epochs_per_sample = epochs_per_sample.size();
  float alpha = initial_alpha;

  for (auto n = 0U; n < n_epochs; n++) {
    worker.set_alpha(alpha);
    worker.set_n(n);
    worker.reseed();
    if (n_threads > 0) {
      RcppPerpendicular::parallel_for(0, n_epochs_per_sample, worker, n_threads,
                                      grain_size);
    } else {
      worker(0, n_epochs_per_sample);
    }
    alpha = initial_alpha * (1.0 - (float(n) / float(n_epochs)));

    if (Progress::check_abort()) {
      progress.cleanup();
      return head_embedding;
    }
    if (verbose) {
      progress.increment();
    }
  }
  return head_embedding;
}

// [[Rcpp::export]]
NumericMatrix optimize_layout_umap(
    NumericMatrix head_embedding, Nullable<NumericMatrix> tail_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float a, float b, float gamma, float initial_alpha,
    float negative_sample_rate, bool approx_pow, bool pcg_rand = true,
    std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool move_other = true, bool verbose = false) {
  // For normal UMAP, tail_embedding is NULL and we want to pass
  // a shallow copy of head_embedding as tail_embedding.
  // When updating new values, tail_embedding is the new coordinate to optimize
  // and gets passed as normal.
  auto head_vec = as<std::vector<float>>(head_embedding);
  std::vector<float> *tail_vec_ptr = nullptr;
  bool delete_tail_ptr = false;
  if (tail_embedding.isNull()) {
    tail_vec_ptr = &head_vec;
  } else {
    tail_vec_ptr =
        new std::vector<float>(as<std::vector<float>>(tail_embedding));
    delete_tail_ptr = true;
  }

  std::vector<float> result;
  if (approx_pow) {
    const uwot::apumap_gradient gradient(a, b, gamma);
    if (move_other) {
      if (pcg_rand) {
        result = optimize_layout<uwot::apumap_gradient, true, pcg_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      } else {
        result = optimize_layout<uwot::apumap_gradient, true, tau_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      }
    } else {
      if (pcg_rand) {
        result = optimize_layout<uwot::apumap_gradient, false, pcg_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      } else {
        result = optimize_layout<uwot::apumap_gradient, false, tau_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      }
    }
  } else {
    const uwot::umap_gradient gradient(a, b, gamma);
    if (move_other) {
      if (pcg_rand) {
        result = optimize_layout<uwot::umap_gradient, true, pcg_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      } else {
        result = optimize_layout<uwot::umap_gradient, true, tau_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      }
    } else {
      if (pcg_rand) {
        result = optimize_layout<uwot::umap_gradient, false, pcg_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      } else {
        result = optimize_layout<uwot::umap_gradient, false, tau_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      }
    }
  }

  if (delete_tail_ptr) {
    delete (tail_vec_ptr);
  }

  return NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                       result.begin());
}

// [[Rcpp::export]]
NumericMatrix optimize_layout_tumap(
    NumericMatrix head_embedding, Nullable<NumericMatrix> tail_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float initial_alpha, float negative_sample_rate, bool pcg_rand = true,
    std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool move_other = true, bool verbose = false) {
  const uwot::tumap_gradient gradient;
  auto head_vec = as<std::vector<float>>(head_embedding);
  std::vector<float> *tail_vec_ptr = nullptr;
  bool delete_tail_ptr = false;
  if (tail_embedding.isNull()) {
    tail_vec_ptr = &head_vec;
  } else {
    tail_vec_ptr =
        new std::vector<float>(as<std::vector<float>>(tail_embedding));
    delete_tail_ptr = true;
  }

  std::vector<float> result;

  if (move_other) {
    if (pcg_rand) {
      result = optimize_layout<uwot::tumap_gradient, true, pcg_factory>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size, verbose);
    } else {
      result = optimize_layout<uwot::tumap_gradient, true, tau_factory>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size, verbose);
    }
  } else {
    if (pcg_rand) {
      result = optimize_layout<uwot::tumap_gradient, false, pcg_factory>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size, verbose);
    } else {
      result = optimize_layout<uwot::tumap_gradient, false, tau_factory>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size, verbose);
    }
  }

  if (delete_tail_ptr) {
    delete (tail_vec_ptr);
  }

  return NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                       result.begin());
}

// [[Rcpp::export]]
NumericMatrix optimize_layout_largevis(
    NumericMatrix head_embedding, const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float gamma, float initial_alpha, float negative_sample_rate,
    bool pcg_rand = true, std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool verbose = false) {
  // We don't support adding extra points for LargeVis, so this is much simpler
  // than the UMAP case
  const uwot::largevis_gradient gradient(gamma);
  auto head_vec = as<std::vector<float>>(head_embedding);

  std::vector<float> result;

  if (pcg_rand) {
    result = optimize_layout<uwot::largevis_gradient, true, pcg_factory>(
        gradient, head_vec, head_vec, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        n_threads, grain_size, verbose);
  } else {
    result = optimize_layout<uwot::largevis_gradient, true, tau_factory>(
        gradient, head_vec, head_vec, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        n_threads, grain_size, verbose);
  }

  return NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                       result.begin());
}
