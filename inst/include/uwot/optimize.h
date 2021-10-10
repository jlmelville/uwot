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

struct Sgd {
  float initial_alpha;
  float alpha;

  Sgd(float alpha) : initial_alpha(alpha), alpha(alpha){};

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    alpha = initial_alpha * (1.0 - (float(epoch) / float(n_epochs)));
  }
};

// Function to decide whether to move both vertices in an edge
// Default empty version does nothing: used in umap_transform when
// some of the vertices should be held fixed
template <bool DoMoveVertex = false>
void update_vec(std::vector<float> &, float, std::size_t, std::size_t) {}

// Specialization to move vertex/update gradient: used in umap when both
// vertices in an edge should be moved
template <>
void update_vec<true>(std::vector<float> &vec, float val, std::size_t i,
                      std::size_t j) {
  vec[i + j] += val;
}

// If DoMoveTailVertex = true, graph is symmetric and head and tail point to the
// same data. So we can just update the head coord i with double the gradient
// now and not worry about updating it when it shows up in the edge list as tail
// point j
template <bool DoMoveTailVertex = true>
void update_head_grad_vec(std::vector<float> &head_grad_vec, float val,
                          std::size_t i) {
  head_grad_vec[i] += 2.0 * val;
}

// Specialization for DoMoveTailVertex = true. In this case the edges are not
// symmetric and the tail embedding should be held fixed, so the head node only
// get one lot of gradient updating
template <>
void update_head_grad_vec<false>(std::vector<float> &head_grad_vec, float val,
                                 std::size_t i) {
  head_grad_vec[i] += val;
}

// DoMoveVertex: true if both ends of a positive edge should be updated
template <bool DoMoveVertex> struct InPlaceUpdate {
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  float alpha; // learning rate

  Sgd opt;

  InPlaceUpdate(std::vector<float> &head_embedding,
                std::vector<float> &tail_embedding, float alpha)
      : head_embedding(head_embedding), tail_embedding(tail_embedding),
        opt(alpha) {}
  void attract(std::size_t dj, std::size_t dk, std::size_t d, float grad_d,
               std::size_t) {
    float update_d = opt.alpha * grad_d;
    head_embedding[dj + d] += update_d;
    // we don't only always want points in the tail to move
    // e.g. if adding new points to an existing embedding
    update_vec<DoMoveVertex>(tail_embedding, -update_d, d, dk);
  }
  void repel(std::size_t dj, std::size_t dk, std::size_t d, float grad_d,
             std::size_t) {
    head_embedding[dj + d] += opt.alpha * grad_d;
    // Python implementation doesn't move the negative sample but as Damrich
    // and Hamprecht (2021) note, it ought to. However they also note it has
    // no qualitative effect on the results. This is presumably because
    // including this repulsion is the same as doubling the negative sample rate
    // which doesn't have a huge effect on going from the default of 5 to 10
    // update_vec<DoMoveVertex>(tail_embedding, opt.alpha * -grad_d, d, dk);
  }

  void epoch_begin(std::size_t, std::size_t) {}
  template <typename Parallel>
  void epoch_end(std::size_t epoch, std::size_t n_epochs, Parallel &) {
    opt.epoch_end(epoch, n_epochs);
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
template <bool DoMoveVertex> struct BatchUpdate {
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;

  Sgd opt;
  const std::size_t head_1d_length;
  std::vector<float> head_gupd;

  BatchUpdate(std::vector<float> &head_embedding,
              std::vector<float> &tail_embedding, float alpha)
      : head_embedding(head_embedding), tail_embedding(tail_embedding),
        opt(alpha), head_1d_length(head_embedding.size()),
        head_gupd(head_1d_length) {}

  void attract(std::size_t dj, std::size_t dk, std::size_t d, float grad_d,
               std::size_t) {
    update_head_grad_vec<DoMoveVertex>(head_gupd, grad_d, dj + d);
  }
  void repel(std::size_t dj, std::size_t dk, std::size_t d, float grad_d,
             std::size_t key) {
    head_gupd[dj + d] += grad_d;
  }

  void epoch_begin(std::size_t, std::size_t) {
    std::fill(head_gupd.begin(), head_gupd.end(), 0.0);
  }

  template <typename Parallel>
  void epoch_end(std::size_t epoch, std::size_t n_epochs, Parallel &parallel) {
    auto worker = [&](std::size_t begin, std::size_t end, std::size_t) {
      for (std::size_t i = begin; i < end; i++) {
        head_embedding[i] += opt.alpha * head_gupd[i];
      }
    };
    parallel.pfor(head_1d_length, worker);

    opt.epoch_end(epoch, n_epochs);
  }
};

template <typename Update, typename Gradient>
void update_attract(Update &update, const Gradient &gradient, std::size_t dj,
                    std::size_t dk, std::size_t ndim, std::vector<float> &disp,
                    std::size_t key) {
  float grad_coeff = grad_attr(gradient, update.head_embedding, dj,
                               update.tail_embedding, dk, ndim, disp);
  for (std::size_t d = 0; d < ndim; d++) {
    update.attract(dj, dk, d, grad_d<Gradient>(disp, d, grad_coeff), key);
  }
}

template <typename Update, typename Gradient>
void update_repel(Update &update, const Gradient &gradient, std::size_t dj,
                  std::size_t dk, std::size_t ndim, std::vector<float> &disp,
                  std::size_t key) {
  float grad_coeff = grad_rep(gradient, update.head_embedding, dj,
                              update.tail_embedding, dk, ndim, disp);
  for (std::size_t d = 0; d < ndim; d++) {
    update.repel(dj, dk, d, grad_d<Gradient>(disp, d, grad_coeff), key);
  }
}

template <typename Worker, typename Progress, typename Parallel>
void optimize_layout(Worker &worker, Progress &progress, unsigned int n_epochs,
                     Parallel &parallel) {
  for (auto n = 0U; n < n_epochs; n++) {
    epoch(worker, n, n_epochs, parallel);

    if (progress.is_aborted()) {
      break;
    }
    progress.report();
  }
}

template <typename Worker, typename Parallel>
void epoch(Worker &worker, std::size_t n, std::size_t n_epochs,
           Parallel &parallel) {
  worker.epoch_begin(n, n_epochs);

  parallel.pfor(worker.n_items, worker);

  worker.epoch_end(n, n_epochs, parallel);
}

template <typename Update, typename Gradient, typename Prng>
void process_edge(Update &update, Gradient &gradient, Sampler &sampler,
                  Prng &prng, const std::vector<unsigned int> &positive_head,
                  const std::vector<unsigned int> &positive_tail,
                  std::size_t ndim, std::size_t tail_nvert, std::size_t epoch,
                  std::size_t edge, std::size_t thread_id,
                  std::vector<float> &disp) {
  if (!sampler.is_sample_edge(edge, epoch)) {
    return;
  }
  const std::size_t dj = ndim * positive_head[edge];
  const std::size_t dk = ndim * positive_tail[edge];
  // j and k are joined by an edge: push them together
  update_attract(update, gradient, dj, dk, ndim, disp, thread_id);

  // Negative sampling step: assume any other point (dkn) is a -ve example
  std::size_t n_neg_samples = sampler.get_num_neg_samples(edge, epoch);
  for (std::size_t p = 0; p < n_neg_samples; p++) {
    const std::size_t dkn = prng(tail_nvert) * ndim;
    if (dj == dkn) {
      continue;
    }
    // push them apart
    update_repel(update, gradient, dj, dkn, ndim, disp, thread_id);
  }
  sampler.next_sample(edge, n_neg_samples);
}

// Gradient: the type of gradient used in the optimization
// Update: type of update to the embedding coordinates
template <typename Gradient, typename Update, typename RngFactory>
struct EdgeWorker {
  int n; // epoch counter
  const Gradient gradient;
  Update &update;
  const std::vector<unsigned int> &positive_head;
  const std::vector<unsigned int> &positive_tail;
  uwot::Sampler sampler;
  std::size_t ndim;
  std::size_t tail_nvert;
  std::size_t n_items;
  std::size_t n_threads;
  RngFactory rng_factory;

  EdgeWorker(const Gradient &gradient, Update &update,
             const std::vector<unsigned int> &positive_head,
             const std::vector<unsigned int> &positive_tail,
             uwot::Sampler &sampler, std::size_t ndim, std::size_t tail_nvert,
             std::size_t n_threads)
      : n(0), gradient(gradient), update(update), positive_head(positive_head),
        positive_tail(positive_tail), sampler(sampler), ndim(ndim),
        tail_nvert(tail_nvert), n_items(positive_head.size()),
        n_threads(std::max(n_threads, std::size_t{1})),
        rng_factory(this->n_threads) {}

  void epoch_begin(std::size_t epoch, std::size_t n_epochs) {
    n = epoch;
    rng_factory.reseed();

    update.epoch_begin(epoch, n_epochs);
  }

  template <typename Parallel>
  void epoch_end(std::size_t epoch, std::size_t n_epochs, Parallel &parallel) {
    update.epoch_end(epoch, n_epochs, parallel);
  }

  void operator()(std::size_t begin, std::size_t end, std::size_t thread_id) {
    // Each window gets its own PRNG state, to prevent locking inside the loop.
    auto prng = rng_factory.create(end);

    // displacement between two points, cost of reallocating inside the loop
    // is noticeable, also cheaper to calculate it once in the d2 calc
    std::vector<float> disp(ndim);
    for (auto i = begin; i < end; i++) {
      process_edge(update, gradient, sampler, prng, positive_head,
                   positive_tail, ndim, tail_nvert, n, i, thread_id, disp);
    }
  }
};

template <typename Gradient, typename Update, typename RngFactory>
struct NodeWorker {
  int n; // epoch counter
  const Gradient gradient;
  Update &update;
  const std::vector<unsigned int> &positive_head;
  const std::vector<unsigned int> &positive_tail;
  const std::vector<unsigned int> &positive_ptr;
  uwot::Sampler sampler;
  std::size_t ndim;
  std::size_t tail_nvert;
  std::size_t n_items;
  RngFactory rng_factory;

  NodeWorker(const Gradient &gradient, Update &update,
             const std::vector<unsigned int> &positive_head,
             const std::vector<unsigned int> &positive_tail,
             const std::vector<unsigned int> &positive_ptr,
             uwot::Sampler &sampler, std::size_t ndim, std::size_t tail_nvert)
      : n(0), gradient(gradient), update(update), positive_head(positive_head),
        positive_tail(positive_tail), positive_ptr(positive_ptr),
        sampler(sampler), ndim(ndim), tail_nvert(tail_nvert),
        n_items(positive_ptr.size() - 1), rng_factory(n_items) {}

  void epoch_begin(std::size_t epoch, std::size_t n_epochs) {
    n = epoch;
    rng_factory.reseed();

    update.epoch_begin(epoch, n_epochs);
  }

  template <typename Parallel>
  void epoch_end(std::size_t epoch, std::size_t n_epochs, Parallel &parallel) {
    update.epoch_end(epoch, n_epochs, parallel);
  }

  void operator()(std::size_t begin, std::size_t end, std::size_t thread_id) {
    std::vector<float> disp(ndim);
    for (auto p = begin; p < end; p++) {
      auto prng = rng_factory.create(p);
      for (auto i = positive_ptr[p]; i < positive_ptr[p + 1]; i++) {
        process_edge(update, gradient, sampler, prng, positive_head,
                     positive_tail, ndim, tail_nvert, n, i, thread_id, disp);
      }
    }
  }
};

} // namespace uwot

#endif // UWOT_OPTIMIZE_H
