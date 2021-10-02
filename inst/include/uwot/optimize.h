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
#include <memory>
#include <utility>

#include "sampler.h"

namespace uwot {

// For normal UMAP, tail_embedding is NULL and we want to pass
// a shallow copy of head_embedding as tail_embedding.
// When updating new values, tail_embedding is the new coordinate to optimize
// and gets passed as normal.
struct Coords {
  std::vector<float> head_embedding;
  std::unique_ptr<std::vector<float>> tail_vec_ptr;

  Coords(std::vector<float> &head_embedding)
      : head_embedding(head_embedding), tail_vec_ptr(nullptr) {}

  Coords(std::vector<float> &head_embedding, std::vector<float> &tail_embedding)
      : head_embedding(head_embedding),
        tail_vec_ptr(new std::vector<float>(tail_embedding)) {}

  auto get_tail_embedding() -> std::vector<float> & {
    if (tail_vec_ptr) {
      return *tail_vec_ptr;
    } else {
      return head_embedding;
    }
  }

  auto get_head_embedding() -> std::vector<float> & { return head_embedding; }
};

// Function to decide whether to move both vertices in an edge
// Default empty version does nothing: used in umap_transform when
// some of the vertices should be held fixed
template <bool DoMoveVertex = false>
void move_other_vertex(std::vector<float> &, float, std::size_t, std::size_t) {}

// Specialization to move the vertex: used in umap when both
// vertices in an edge should be moved
template <>
void move_other_vertex<true>(std::vector<float> &embedding, float update_d,
                             std::size_t i, std::size_t nrj) {
  embedding[nrj + i] -= update_d;
}

// DoMoveVertex: true if both ends of a positive edge should be updated
template <bool DoMoveVertex> struct InPlaceUpdate {
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;

  InPlaceUpdate(std::vector<float> &head_embedding,
                std::vector<float> &tail_embedding)
      : head_embedding(head_embedding), tail_embedding(tail_embedding) {}
  void attract(std::size_t dj, std::size_t dk, std::size_t d, float update_d) {
    head_embedding[dj + d] += update_d;
    // we don't only always want points in the tail to move
    // e.g. if adding new points to an existing embedding
    move_other_vertex<DoMoveVertex>(tail_embedding, update_d, d, dk);
  }
  void repel(std::size_t dj, std::size_t dk, std::size_t d, float update_d) {
    head_embedding[dj + d] += update_d;
  }
};

// Gradient: the type of gradient used in the optimization
// Update: type of update to the embedding coordinates
template <typename Gradient, typename Update, typename RngFactory>
struct SgdWorker {
  int n; // epoch counter
  float alpha;
  const Gradient gradient;
  Update update;
  const std::vector<unsigned int> positive_head;
  const std::vector<unsigned int> positive_tail;
  uwot::Sampler sampler;
  const std::vector<float> &head_embedding;
  const std::vector<float> &tail_embedding;
  std::size_t ndim;
  std::size_t tail_nvert;
  RngFactory rng_factory;

  SgdWorker(const Gradient &gradient, Update &update,
            std::vector<unsigned int> positive_head,
            std::vector<unsigned int> positive_tail, uwot::Sampler &sampler,
            const std::vector<float> &head_embedding,
            const std::vector<float> &tail_embedding, std::size_t ndim)
      :

        n(0), alpha(0.0), gradient(gradient), update(update),
        positive_head(std::move(positive_head)),
        positive_tail(std::move(positive_tail)),

        sampler(sampler),

        head_embedding(head_embedding), tail_embedding(tail_embedding),
        ndim(ndim), tail_nvert(tail_embedding.size() / ndim),

        rng_factory() {}

  void operator()(std::size_t begin, std::size_t end) {
    // Each window gets its own PRNG state, to prevent locking inside the loop.
    auto prng = rng_factory.create(end);

    // the displacement between two points
    std::vector<float> disp(ndim);
    for (auto i = begin; i < end; i++) {
      if (!sampler.is_sample_edge(i, n)) {
        continue;
      }
      const std::size_t dj = ndim * positive_head[i];

      // j and k are joined by an edge: push them together
      const std::size_t dk = ndim * positive_tail[i];
      float grad_coeff = grad_attr(gradient, head_embedding, dj, tail_embedding,
                                   dk, ndim, disp);
      for (std::size_t d = 0; d < ndim; d++) {
        float update_d = alpha * grad_d<Gradient>(disp, d, grad_coeff);
        update.attract(dj, dk, d, update_d);
      }

      // Negative sampling step: assume any other point is a negative example
      // Push them apart
      std::size_t n_neg_samples = sampler.get_num_neg_samples(i, n);
      for (std::size_t p = 0; p < n_neg_samples; p++) {
        const std::size_t dkn = prng(tail_nvert) * ndim;
        if (dj == dkn) {
          continue;
        }
        float grad_coeff = grad_rep(gradient, head_embedding, dj,
                                    tail_embedding, dkn, ndim, disp);
        for (std::size_t d = 0; d < ndim; d++) {
          float update_d = alpha * grad_d<Gradient>(disp, d, grad_coeff);
          update.repel(dj, dkn, d, update_d);
        }
      }
      sampler.next_sample(i, n_neg_samples);
    }
  }

  void reseed() { this->rng_factory.reseed(); }
};
} // namespace uwot

#endif // UWOT_OPTIMIZE_H
