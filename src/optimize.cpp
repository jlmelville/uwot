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
void optimize_layout(const T &gradient, std::vector<float> &head_embedding,
                     std::vector<float> &tail_embedding,
                     const std::vector<unsigned int> &positive_head,
                     const std::vector<unsigned int> &positive_tail,
                     unsigned int n_epochs, unsigned int n_vertices,
                     const std::vector<float> &epochs_per_sample,
                     float initial_alpha, float negative_sample_rate,
                     std::size_t n_threads = 0, std::size_t grain_size = 1,
                     bool verbose = false) {
  const std::size_t ndim = head_embedding.size() / n_vertices;

  uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);
  uwot::InPlaceUpdate<DoMove> update(head_embedding, tail_embedding,
                                     initial_alpha);
  // uwot::BatchUpdate<DoMove> update(head_embedding, tail_embedding,
  //                                  initial_alpha, n_threads);
  uwot::SgdWorker<T, decltype(update), RandFactory> worker(
      gradient, update, positive_head, positive_tail, sampler, ndim,
      tail_embedding.size() / ndim);

  Progress progress(n_epochs, verbose);
  const auto n_edges = epochs_per_sample.size();

  for (auto n = 0U; n < n_epochs; n++) {
    worker.epoch_begin(n, n_epochs);

    RcppPerpendicular::pfor(0, n_edges, worker, n_threads, grain_size);

    worker.epoch_end(n, n_epochs);

    if (Progress::check_abort()) {
      progress.cleanup();
      break;
    }
    if (verbose) {
      progress.increment();
    }
  }
}

struct UmapFactory {
  bool move_other;
  bool pcg_rand;
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  const std::vector<unsigned int> &positive_head;
  const std::vector<unsigned int> &positive_tail;
  unsigned int n_epochs;
  unsigned int n_vertices;
  const std::vector<float> &epochs_per_sample;
  float initial_alpha;
  float negative_sample_rate;
  std::size_t n_threads;
  std::size_t grain_size;
  bool verbose;

  UmapFactory(bool move_other, bool pcg_rand,
              std::vector<float> &head_embedding,
              std::vector<float> &tail_embedding,
              const std::vector<unsigned int> &positive_head,
              const std::vector<unsigned int> &positive_tail,
              unsigned int n_epochs, unsigned int n_vertices,
              const std::vector<float> &epochs_per_sample, float initial_alpha,
              float negative_sample_rate, std::size_t n_threads,
              std::size_t grain_size, bool verbose)
      : move_other(move_other), pcg_rand(pcg_rand),
        head_embedding(head_embedding), tail_embedding(tail_embedding),
        positive_head(positive_head), positive_tail(positive_tail),
        n_epochs(n_epochs), n_vertices(n_vertices),
        epochs_per_sample(epochs_per_sample), initial_alpha(initial_alpha),
        negative_sample_rate(negative_sample_rate), n_threads(n_threads),
        grain_size(grain_size), verbose(verbose) {}

  template <typename T> void create(const T &gradient) {
    if (move_other) {
      if (pcg_rand) {
        optimize_layout<T, true, pcg_factory>(
            gradient, head_embedding, tail_embedding, positive_head,
            positive_tail, n_epochs, n_vertices, epochs_per_sample,
            initial_alpha, negative_sample_rate, n_threads, grain_size,
            verbose);
      } else {
        optimize_layout<T, true, tau_factory>(
            gradient, head_embedding, tail_embedding, positive_head,
            positive_tail, n_epochs, n_vertices, epochs_per_sample,
            initial_alpha, negative_sample_rate, n_threads, grain_size,
            verbose);
      }
    } else {
      if (pcg_rand) {
        optimize_layout<T, false, pcg_factory>(
            gradient, head_embedding, tail_embedding, positive_head,
            positive_tail, n_epochs, n_vertices, epochs_per_sample,
            initial_alpha, negative_sample_rate, n_threads, grain_size,
            verbose);
      } else {
        optimize_layout<T, false, tau_factory>(
            gradient, head_embedding, tail_embedding, positive_head,
            positive_tail, n_epochs, n_vertices, epochs_per_sample,
            initial_alpha, negative_sample_rate, n_threads, grain_size,
            verbose);
      }
    }
  }
};

auto r_to_coords(NumericMatrix head_embedding,
                 Nullable<NumericMatrix> tail_embedding) -> uwot::Coords {
  auto head_vec = as<std::vector<float>>(head_embedding);
  if (tail_embedding.isNull()) {
    return uwot::Coords(head_vec);
  } else {
    auto tail_vec = as<std::vector<float>>(tail_embedding);
    return uwot::Coords(head_vec, tail_vec);
  }
}

auto r_to_coords(NumericMatrix head_embedding) -> uwot::Coords {
  auto head_vec = as<std::vector<float>>(head_embedding);
  return uwot::Coords(head_vec);
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

  auto coords = r_to_coords(head_embedding, tail_embedding);

  UmapFactory umap_factory(
      move_other, pcg_rand, coords.get_head_embedding(),
      coords.get_tail_embedding(), positive_head, positive_tail, n_epochs,
      n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
      n_threads, grain_size, verbose);
  if (approx_pow) {
    const uwot::apumap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  } else {
    const uwot::umap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  }

  return NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                       coords.get_head_embedding().begin());
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

  auto coords = r_to_coords(head_embedding, tail_embedding);

  const uwot::tumap_gradient gradient;
  UmapFactory umap_factory(
      move_other, pcg_rand, coords.get_head_embedding(),
      coords.get_tail_embedding(), positive_head, positive_tail, n_epochs,
      n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
      n_threads, grain_size, verbose);
  umap_factory.create(gradient);

  return NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                       coords.get_head_embedding().begin());
}

// [[Rcpp::export]]
NumericMatrix optimize_layout_largevis(
    NumericMatrix head_embedding, const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float gamma, float initial_alpha, float negative_sample_rate,
    bool pcg_rand = true, std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool verbose = false) {

  auto coords = r_to_coords(head_embedding);

  UmapFactory umap_factory(
      true, pcg_rand, coords.get_head_embedding(), coords.get_tail_embedding(),
      positive_head, positive_tail, n_epochs, n_vertices, epochs_per_sample,
      initial_alpha, negative_sample_rate, n_threads, grain_size, verbose);
  const uwot::largevis_gradient gradient(gamma);
  umap_factory.create(gradient);

  return NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                       coords.get_head_embedding().begin());
}
