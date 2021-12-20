// BSD 2-Clause License
//
// Copyright 2021 James Melville
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

#ifndef UWOT_UPDATE_H
#define UWOT_UPDATE_H

#include "gradient.h"
#include "optimize.h"
#include "sampler.h"

namespace uwot {

// (User-supplied) callback to be invoked at the end of each epoch
// e.g. to plot or save coordinates for diagnostic purposes
struct EpochCallback {
  virtual void operator()(std::size_t epoch, std::size_t n_epochs,
                          const std::vector<float> &head_embedding,
                          const std::vector<float> &tail_embedding) = 0;
  virtual ~EpochCallback() {}
};

struct DoNothingCallback : EpochCallback {
  void operator()(std::size_t, std::size_t, const std::vector<float> &,
                  const std::vector<float> &) override {}
};

template <typename Update, typename Gradient, typename Prng>
inline void
process_edge(Update &update, Gradient &gradient, Sampler &sampler, Prng &prng,
             const std::vector<unsigned int> &positive_head,
             const std::vector<unsigned int> &positive_tail, std::size_t ndim,
             std::size_t n_tail_vertices, std::size_t edge,
             std::size_t thread_id, std::vector<float> &disp) {
  if (!sampler.is_sample_edge(edge)) {
    return;
  }
  const std::size_t dj = ndim * positive_head[edge];
  const std::size_t dk = ndim * positive_tail[edge];
  // j and k are joined by an edge: push them together
  update_attract(update, gradient, dj, dk, ndim, disp, thread_id);

  // Negative sampling step: assume any other point (dkn) is a -ve example
  std::size_t n_neg_samples = sampler.get_num_neg_samples(edge);
  for (std::size_t p = 0; p < n_neg_samples; p++) {
    const std::size_t dkn = prng(n_tail_vertices) * ndim;
    if (dj == dkn) {
      continue;
    }
    // push them apart
    update_repel(update, gradient, dj, dkn, ndim, disp, thread_id);
  }
  sampler.next_sample(edge, n_neg_samples);
}

template <typename Update, typename Gradient>
inline void update_attract(Update &update, const Gradient &gradient,
                           std::size_t dj, std::size_t dk, std::size_t ndim,
                           std::vector<float> &disp, std::size_t key) {
  float grad_coeff = grad_attr(gradient, update.head_embedding, dj,
                               update.tail_embedding, dk, ndim, disp);
  for (std::size_t d = 0; d < ndim; d++) {
    update.attract(dj, dk, d, grad_d(gradient, disp, d, grad_coeff), key);
  }
}

template <typename Update, typename Gradient>
inline void update_repel(Update &update, const Gradient &gradient,
                         std::size_t dj, std::size_t dk, std::size_t ndim,
                         std::vector<float> &disp, std::size_t key) {
  float grad_coeff = grad_rep(gradient, update.head_embedding, dj,
                              update.tail_embedding, dk, ndim, disp);
  for (std::size_t d = 0; d < ndim; d++) {
    update.repel(dj, dk, d, grad_d(gradient, disp, d, grad_coeff), key);
  }
}

// Function to decide whether to move both vertices in an edge
// Default empty version does nothing: used in umap_transform when
// some of the vertices should be held fixed
template <bool DoMoveVertex = false>
inline void update_vec(std::vector<float> &, float, std::size_t, std::size_t) {}

// Specialization to move vertex/update gradient: used in umap when both
// vertices in an edge should be moved
template <>
inline void update_vec<true>(std::vector<float> &vec, float val, std::size_t i,
                             std::size_t j) {
  vec[i + j] += val;
}

// If DoMoveTailVertex = true, graph is symmetric and head and tail point to the
// same data. So we can just update the head coord i with double the gradient
// now and not worry about updating it when it shows up in the edge list as tail
// point j
template <bool DoMoveTailVertex = true>
inline void update_head_grad_vec(std::vector<float> &vec, std::size_t i,
                                 float val) {
  vec[i] += 2.0 * val;
}

// Specialization for DoMoveTailVertex = true. In this case the edges are not
// symmetric and the tail embedding should be held fixed, so the head node only
// get one lot of gradient updating
template <>
inline void update_head_grad_vec<false>(std::vector<float> &vec, std::size_t i,
                                        float val) {
  vec[i] += val;
}

// DoMoveVertex: true if both ends of a positive edge should be updated
template <bool DoMoveVertex> struct InPlaceUpdate {
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;

  Sgd opt;
  std::unique_ptr<EpochCallback> epoch_callback;

  InPlaceUpdate(std::vector<float> &head_embedding,
                std::vector<float> &tail_embedding, float alpha,
                EpochCallback *epoch_callback)
      : head_embedding(head_embedding), tail_embedding(tail_embedding),
        opt(alpha), epoch_callback(std::move(epoch_callback)) {}

  inline void attract(std::size_t dj, std::size_t dk, std::size_t d,
                      float grad_d, std::size_t) {
    float update_d = opt.alpha * grad_d;
    head_embedding[dj + d] += update_d;
    // we don't only always want points in the tail to move
    // e.g. if adding new points to an existing embedding
    update_vec<DoMoveVertex>(tail_embedding, -update_d, d, dk);
  }
  inline void repel(std::size_t dj, std::size_t dk, std::size_t d, float grad_d,
                    std::size_t) {
    head_embedding[dj + d] += opt.alpha * grad_d;
    // Python implementation doesn't move the negative sample but as Damrich
    // and Hamprecht (2021) note, it ought to. However they also note it has
    // no qualitative effect on the results. This is presumably because
    // including this repulsion is the same as doubling the negative sample rate
    // which doesn't have a huge effect on going from the default of 5 to 10
    // update_vec<DoMoveVertex>(tail_embedding, opt.alpha * -grad_d, d, dk);
  }

  inline void epoch_begin(std::size_t, std::size_t) {}
  template <typename Parallel>
  void epoch_end(std::size_t epoch, std::size_t n_epochs, Parallel &) {
    opt.epoch_end(epoch, n_epochs);
    (*epoch_callback)(epoch, n_epochs, head_embedding, tail_embedding);
  }
};

// 1. When DoMoveVertex is true then we want to update the head and tail nodes
// of an edge. In this case the head and tail coordinates point to the same data
// so it doesn't matter whether we calculate the gradient for or update the
// coordinates in head or tail.
// 2. When DoMoveVertex is false then the head and tail coordinates point to
// different data. The tail coordinates are fixed in this case, so again they
// do not move. Hence both so in both cases we only ever need to update the head
// coordinates.
template <bool DoMoveVertex, typename Opt> struct BatchUpdate {
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  Opt &opt;
  std::vector<float> gradient;
  std::unique_ptr<EpochCallback> epoch_callback;

  BatchUpdate(std::vector<float> &head_embedding,
              std::vector<float> &tail_embedding, Opt &opt,
              EpochCallback *epoch_callback)
      : head_embedding(head_embedding), tail_embedding(tail_embedding),
        opt(opt), gradient(head_embedding.size()),
        epoch_callback(std::move(epoch_callback)) {}

  BatchUpdate(BatchUpdate &&other)
      : head_embedding(other.head_embedding),
        tail_embedding(other.tail_embedding), opt(other.opt),
        gradient(other.head_embedding.size()),
        epoch_callback(std::move(other.epoch_callback)) {}

  inline void attract(std::size_t dj, std::size_t dk, std::size_t d,
                      float grad_d, std::size_t) {
    update_head_grad_vec<DoMoveVertex>(gradient, dj + d, grad_d);
  }

  inline void repel(std::size_t dj, std::size_t dk, std::size_t d, float grad_d,
                    std::size_t) {
    gradient[dj + d] += grad_d;
  }

  void epoch_begin(std::size_t, std::size_t) {
    std::fill(gradient.begin(), gradient.end(), 0.0f);
  }

  template <typename Parallel>
  void epoch_end(std::size_t epoch, std::size_t n_epochs, Parallel &parallel) {
    auto worker = [&](std::size_t begin, std::size_t end, std::size_t) {
      for (std::size_t i = begin; i < end; i++) {
        // head_embedding[i] += opt.alpha * gradient[i];
        opt.update(head_embedding, gradient, i);
      }
    };
    parallel.pfor(head_embedding.size(), worker);

    opt.epoch_end(epoch, n_epochs);
    (*epoch_callback)(epoch, n_epochs, head_embedding, tail_embedding);
  }
};

} // namespace uwot

#endif // UWOT_UPDATE_H