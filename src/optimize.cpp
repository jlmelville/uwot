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
#include <memory>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <utility>

#include "gradient.h"
#include "sampler.h"
#include "tauprng.h"

// Must come after any include that transitively include dqrng
#include "RcppPerpendicular.h"

// Function to decide whether to move both vertices in an edge
// Default empty version does nothing: used in umap_transform when
// some of the vertices should be held fixed
template <bool DoMoveVertex = false>
void move_other_vertex(std::vector<float> &embedding, const float grad_d,
                       const std::size_t i, const std::size_t nrj) {}

// Specialization to move the vertex: used in umap when both
// vertices in an edge should be moved
template <>
void move_other_vertex<true>(std::vector<float> &embedding, const float grad_d,
                             const std::size_t i, const std::size_t nrj) {
  embedding[nrj + i] -= grad_d;
}

auto clamp(const float v, const float lo, const float hi) -> const float {
  const float t = v < lo ? lo : v;
  return t > hi ? hi : t;
}

// Gradient: the type of gradient used in the optimization
// DoMoveVertex: true if both ends of a positive edge should be updated
template <typename Gradient, bool DoMoveVertex = true,
          typename RngFactory = pcg_factory>
struct SgdWorker {
  int n; // epoch counter
  float alpha;
  const Gradient gradient;
  const std::vector<unsigned int> positive_head;
  const std::vector<unsigned int> positive_tail;
  Sampler sampler;
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  const std::size_t ndim;
  const std::size_t head_nvert;
  const std::size_t tail_nvert;
  const float dist_eps;
  RngFactory rng_factory;

  SgdWorker(const Gradient &gradient, std::vector<unsigned int> positive_head,
            std::vector<unsigned int> positive_tail, Sampler &sampler,
            std::vector<float> &head_embedding,
            std::vector<float> &tail_embedding, const std::size_t ndim)
      :

        n(0), alpha(0.0), gradient(gradient),
        positive_head(std::move(positive_head)),
        positive_tail(std::move(positive_tail)),

        sampler(sampler),

        head_embedding(head_embedding), tail_embedding(tail_embedding),
        ndim(ndim), head_nvert(head_embedding.size() / ndim),
        tail_nvert(tail_embedding.size() / ndim),
        dist_eps(std::numeric_limits<float>::epsilon()),

        rng_factory() {}

  void operator()(std::size_t begin, std::size_t end) {
    // Each window gets its own PRNG state, to prevent locking inside the loop.
    auto prng = rng_factory.create(end);

    std::vector<float> dys(ndim);
    for (std::size_t i = begin; i < end; i++) {
      if (!sampler.is_sample_edge(i, n)) {
        continue;
      }
      const std::size_t dj = ndim * positive_head[i];
      const std::size_t dk = ndim * positive_tail[i];

      float dist_squared = 0.0;
      for (std::size_t d = 0; d < ndim; d++) {
        const float diff = head_embedding[dj + d] - tail_embedding[dk + d];
        dys[d] = diff;
        dist_squared += diff * diff;
      }
      dist_squared = (std::max)(dist_eps, dist_squared);
      const float grad_coeff = gradient.grad_attr(dist_squared);

      for (std::size_t d = 0; d < ndim; d++) {
        const float grad_d =
            alpha *
            clamp(grad_coeff * dys[d], Gradient::clamp_lo, Gradient::clamp_hi);
        head_embedding[dj + d] += grad_d;
        move_other_vertex<DoMoveVertex>(tail_embedding, grad_d, d, dk);
      }

      const std::size_t n_neg_samples = sampler.get_num_neg_samples(i, n);
      for (std::size_t p = 0; p < n_neg_samples; p++) {
        const std::size_t dkn = prng(tail_nvert) * ndim;
        if (dj == dkn) {
          continue;
        }
        float dist_squared = 0.0;
        for (std::size_t d = 0; d < ndim; d++) {
          const float diff = head_embedding[dj + d] - tail_embedding[dkn + d];
          dys[d] = diff;
          dist_squared += diff * diff;
        }
        dist_squared = (std::max)(dist_eps, dist_squared);
        const float grad_coeff = gradient.grad_rep(dist_squared);

        for (std::size_t d = 0; d < ndim; d++) {
          const float grad_d =
              alpha * clamp(grad_coeff * dys[d], Gradient::clamp_lo,
                            Gradient::clamp_hi);
          head_embedding[dj + d] += grad_d;
        }
      }
      sampler.next_sample(i, n_neg_samples);
    }
  }

  void set_n(int n) { this->n = n; }

  void set_alpha(float alpha) { this->alpha = alpha; }

  void reseed() { this->rng_factory.reseed(); }
};

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
  Sampler sampler(epochs_per_sample, negative_sample_rate);

  SgdWorker<T, DoMove, RandFactory> worker(
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
Rcpp::NumericMatrix optimize_layout_umap(
    Rcpp::NumericMatrix head_embedding,
    Rcpp::Nullable<Rcpp::NumericMatrix> tail_embedding,
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
  auto head_vec = Rcpp::as<std::vector<float>>(head_embedding);
  std::vector<float> *tail_vec_ptr = nullptr;
  bool delete_tail_ptr = false;
  if (tail_embedding.isNull()) {
    tail_vec_ptr = &head_vec;
  } else {
    tail_vec_ptr =
        new std::vector<float>(Rcpp::as<std::vector<float>>(tail_embedding));
    delete_tail_ptr = true;
  }

  std::vector<float> result;
  if (approx_pow) {
    const apumap_gradient gradient(a, b, gamma);
    if (move_other) {
      if (pcg_rand) {
        result = optimize_layout<apumap_gradient, true, pcg_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      } else {
        result = optimize_layout<apumap_gradient, true, tau_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      }
    } else {
      if (pcg_rand) {
        result = optimize_layout<apumap_gradient, false, pcg_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      } else {
        result = optimize_layout<apumap_gradient, false, tau_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      }
    }
  } else {
    const umap_gradient gradient(a, b, gamma);
    if (move_other) {
      if (pcg_rand) {
        result = optimize_layout<umap_gradient, true, pcg_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      } else {
        result = optimize_layout<umap_gradient, true, tau_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      }
    } else {
      if (pcg_rand) {
        result = optimize_layout<umap_gradient, false, pcg_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      } else {
        result = optimize_layout<umap_gradient, false, tau_factory>(
            gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
            n_epochs, n_vertices, epochs_per_sample, initial_alpha,
            negative_sample_rate, n_threads, grain_size, verbose);
      }
    }
  }

  if (delete_tail_ptr) {
    delete (tail_vec_ptr);
  }

  return Rcpp::NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                             result.begin());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix optimize_layout_tumap(
    Rcpp::NumericMatrix head_embedding,
    Rcpp::Nullable<Rcpp::NumericMatrix> tail_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float initial_alpha, float negative_sample_rate, bool pcg_rand = true,
    std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool move_other = true, bool verbose = false) {
  const tumap_gradient gradient;
  auto head_vec = Rcpp::as<std::vector<float>>(head_embedding);
  std::vector<float> *tail_vec_ptr = nullptr;
  bool delete_tail_ptr = false;
  if (tail_embedding.isNull()) {
    tail_vec_ptr = &head_vec;
  } else {
    tail_vec_ptr =
        new std::vector<float>(Rcpp::as<std::vector<float>>(tail_embedding));
    delete_tail_ptr = true;
  }

  std::vector<float> result;

  if (move_other) {
    if (pcg_rand) {
      result = optimize_layout<tumap_gradient, true, pcg_factory>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size, verbose);
    } else {
      result = optimize_layout<tumap_gradient, true, tau_factory>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size, verbose);
    }
  } else {
    if (pcg_rand) {
      result = optimize_layout<tumap_gradient, false, pcg_factory>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size, verbose);
    } else {
      result = optimize_layout<tumap_gradient, false, tau_factory>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, n_threads, grain_size, verbose);
    }
  }

  if (delete_tail_ptr) {
    delete (tail_vec_ptr);
  }

  return Rcpp::NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                             result.begin());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix optimize_layout_largevis(
    Rcpp::NumericMatrix head_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail, unsigned int n_epochs,
    unsigned int n_vertices, const std::vector<float> epochs_per_sample,
    float gamma, float initial_alpha, float negative_sample_rate,
    bool pcg_rand = true, std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool verbose = false) {
  // We don't support adding extra points for LargeVis, so this is much simpler
  // than the UMAP case
  const largevis_gradient gradient(gamma);
  auto head_vec = Rcpp::as<std::vector<float>>(head_embedding);

  std::vector<float> result;

  if (pcg_rand) {
    result = optimize_layout<largevis_gradient, true, pcg_factory>(
        gradient, head_vec, head_vec, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        n_threads, grain_size, verbose);
  } else {
    result = optimize_layout<largevis_gradient, true, tau_factory>(
        gradient, head_vec, head_vec, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        n_threads, grain_size, verbose);
  }

  return Rcpp::NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                             result.begin());
}
