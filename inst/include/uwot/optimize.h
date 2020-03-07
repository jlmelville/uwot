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

#ifndef UWOT_OPTIMIZE_H
#define UWOT_OPTIMIZE_H

#include <limits>
#include <utility>

#include "sampler.h"

namespace uwot {

// Function to decide whether to move both vertices in an edge
// Default empty version does nothing: used in umap_transform when
// some of the vertices should be held fixed
template <bool DoMoveVertex = false>
void move_other_vertex(std::vector<float> &, float, std::size_t, std::size_t) {}

// Specialization to move the vertex: used in umap when both
// vertices in an edge should be moved
template <>
void move_other_vertex<true>(std::vector<float> &embedding, float grad_d,
                             std::size_t i, std::size_t nrj) {
  embedding[nrj + i] -= grad_d;
}

inline auto clamp(float v, float lo, float hi) -> float {
  float t = v < lo ? lo : v;
  return t > hi ? hi : t;
}

// Gradient: the type of gradient used in the optimization
// DoMoveVertex: true if both ends of a positive edge should be updated
template <typename Gradient, bool DoMoveVertex, typename RngFactory>
struct SgdWorker {
  int n; // epoch counter
  float alpha;
  const Gradient gradient;
  const std::vector<unsigned int> positive_head;
  const std::vector<unsigned int> positive_tail;
  uwot::Sampler sampler;
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  std::size_t ndim;
  std::size_t head_nvert;
  std::size_t tail_nvert;
  float dist_eps;
  RngFactory rng_factory;

  SgdWorker(const Gradient &gradient, std::vector<unsigned int> positive_head,
            std::vector<unsigned int> positive_tail, uwot::Sampler &sampler,
            std::vector<float> &head_embedding,
            std::vector<float> &tail_embedding, std::size_t ndim)
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
    for (auto i = begin; i < end; i++) {
      if (!sampler.is_sample_edge(i, n)) {
        continue;
      }
      std::size_t dj = ndim * positive_head[i];
      std::size_t dk = ndim * positive_tail[i];

      float dist_squared = 0.0;
      for (std::size_t d = 0; d < ndim; d++) {
        float diff = head_embedding[dj + d] - tail_embedding[dk + d];
        dys[d] = diff;
        dist_squared += diff * diff;
      }
      dist_squared = (std::max)(dist_eps, dist_squared);
      float grad_coeff = gradient.grad_attr(dist_squared);

      for (std::size_t d = 0; d < ndim; d++) {
        float grad_d = alpha * clamp(grad_coeff * dys[d], Gradient::clamp_lo,
                                     Gradient::clamp_hi);
        head_embedding[dj + d] += grad_d;
        move_other_vertex<DoMoveVertex>(tail_embedding, grad_d, d, dk);
      }

      std::size_t n_neg_samples = sampler.get_num_neg_samples(i, n);
      for (std::size_t p = 0; p < n_neg_samples; p++) {
        std::size_t dkn = prng(tail_nvert) * ndim;
        if (dj == dkn) {
          continue;
        }
        float dist_squared = 0.0;
        for (std::size_t d = 0; d < ndim; d++) {
          float diff = head_embedding[dj + d] - tail_embedding[dkn + d];
          dys[d] = diff;
          dist_squared += diff * diff;
        }
        dist_squared = (std::max)(dist_eps, dist_squared);
        float grad_coeff = gradient.grad_rep(dist_squared);

        for (std::size_t d = 0; d < ndim; d++) {
          float grad_d = alpha * clamp(grad_coeff * dys[d], Gradient::clamp_lo,
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
} // namespace uwot

#endif // UWOT_OPTIMIZE_H
