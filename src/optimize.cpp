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
#include <random>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "gradient.h"
#include "tauprng.h"

// Function to decide whether to move both vertices in an edge
// Default empty version does nothing: used in umap_transform when
// some of the vertices should be held fixed
template <bool DoMoveVertex = false>
void move_other_vertex(arma::mat& embedding,
                       const double& grad_d,
                       const arma::uword& k,
                       const arma::uword& d) {
}

// Specialization to move the vertex: used in umap when both
// vertices in an edge should be moved
template <>
void move_other_vertex<true>(arma::mat& embedding,
                             const double& grad_d,
                             const arma::uword& k,
                             const arma::uword& d) {
  embedding.at(k, d) -= grad_d;
}

// Gradient: the type of gradient used in the optimization
// DoMoveVertex: true if both ends of a positive edge should be updated
template <typename Gradient,
          bool DoMoveVertex = true>
struct SgdWorker : public RcppParallel::Worker {

  int n; // epoch counter
  double alpha;
  const Gradient gradient;
  const arma::uvec positive_head;
  const arma::uvec positive_tail;
  const arma::vec epochs_per_sample;

  arma::vec epoch_of_next_sample;
  const arma::vec epochs_per_negative_sample;
  arma::vec epoch_of_next_negative_sample;
  arma::mat& head_embedding;
  arma::mat& tail_embedding;
  unsigned int n_vertices;
  const arma::uword nrow;
  const arma::uword ncol;
  tthread::mutex mutex;
  std::mt19937 rng;
  std::uniform_int_distribution<long> gen;
  const double dist_eps;

  SgdWorker(
    const Gradient& gradient,
    const arma::uvec& positive_head,
    const arma::uvec& positive_tail,
    const arma::vec& epochs_per_sample,

    arma::vec& epoch_of_next_sample,
    const arma::vec& epochs_per_negative_sample,
    arma::vec& epoch_of_next_negative_sample,

    arma::mat& head_embedding,
    arma::mat& tail_embedding,
    unsigned int n_vertices,
    const arma::uword nrow,
    const arma::uword ncol,
    unsigned int seed) :

    n(0), alpha(0.0), gradient(gradient),
    positive_head(positive_head), positive_tail(positive_tail),
    epochs_per_sample(epochs_per_sample),
    epoch_of_next_sample(epoch_of_next_sample),
    epochs_per_negative_sample(epochs_per_negative_sample),
    epoch_of_next_negative_sample(epoch_of_next_negative_sample),
    head_embedding(head_embedding),
    tail_embedding(tail_embedding),
    n_vertices(n_vertices), nrow(nrow), ncol(ncol),
    rng(seed), gen(-2147483647, 2147483646),
    dist_eps(std::numeric_limits<double>::epsilon())
  {  }

  void operator()(std::size_t begin, std::size_t end) {
    // Each window gets its own fast PRNG state, so no locking needed inside the loop.
    // Want separate seeds though, so seed the fast PRNG with three random numbers
    // taken from the mt19937 generator, which is shared across windows, so is locked.
    // Probably this is a bit of a waste of time:
    // Could use the mt19937 seed, the begin and the epoch number as seeds?
    // Doesn't waste much time, though.
    long s1, s2, s3;
    {
      tthread::lock_guard<tthread::mutex> guard(mutex);
      s1 = gen(rng);
      s2 = gen(rng); // technically this needs to always be > 7
      s3 = gen(rng); // should be > 15
    }
    tau_prng prng(s1, s2, s3);

    for (std::size_t i = begin; i < end; i++) {
      if (epoch_of_next_sample[i] <= n) {
        arma::uword j = positive_head[i];
        arma::uword k = positive_tail[i];

        const double dist_squared = std::max(rdist(
          head_embedding, tail_embedding, j, k, ncol), dist_eps);
        const double grad_coeff = gradient.grad_attr(dist_squared);

        for (arma::uword d = 0; d < ncol; d++) {
          double grad_d = alpha * clip(grad_coeff * 
            (head_embedding.at(j, d) - tail_embedding.at(k, d)), 
            gradient.clip_max);
          head_embedding.at(j, d) += grad_d;
          move_other_vertex<DoMoveVertex>(tail_embedding, grad_d, k, d);
        }

        epoch_of_next_sample[i] += epochs_per_sample[i];

        unsigned int n_neg_samples = static_cast<unsigned int>(
          (n - epoch_of_next_negative_sample[i]) / 
            epochs_per_negative_sample[i]);

        for (unsigned int p = 0; p < n_neg_samples; p++) {
          arma::uword k = prng() % n_vertices;

          if (j == k) {
            continue;
          }

          const double dist_squared = std::max(
            rdist(head_embedding, tail_embedding, j, k, ncol), dist_eps);
          const double grad_coeff = gradient.grad_rep(dist_squared);

          for (arma::uword d = 0; d < ncol; d++) {
            head_embedding.at(j, d) += alpha * clip(grad_coeff * 
              (head_embedding.at(j, d) - tail_embedding.at(k, d)), 
              gradient.clip_max);
          }
        }
        epoch_of_next_negative_sample[i] += 
          n_neg_samples * epochs_per_negative_sample[i];
      }
    }
  }

  void set_n(int n) {
    this->n = n;
  }

  void set_alpha(double alpha) {
    this->alpha = alpha;
  }

  double clip(double val, double clip_max) {
    return std::max(std::min(val, clip_max), -clip_max);
  }

  double rdist(const arma::mat& m,
               const arma::mat& n,
               const arma::uword a,
               const arma::uword b,
               const arma::uword ndim) {
    double sum = 0.0;
    for (arma::uword i = 0; i < ndim; i++) {
      sum += (m.at(a, i) - n.at(b, i)) * (m.at(a, i) - n.at(b, i));
    }

    return sum;
  }

};


template<typename T, bool DoMove = true>
arma::mat optimize_layout(const T& gradient,
                          arma::mat& head_embedding,
                          arma::mat& tail_embedding,
                          const arma::uvec& positive_head,
                          const arma::uvec& positive_tail,
                          unsigned int n_epochs, unsigned int n_vertices,
                          const arma::vec& epochs_per_sample,
                          double initial_alpha,
                          double negative_sample_rate,
                          unsigned int seed,
                          bool parallelize = true,
                          std::size_t grain_size = 1,
                          bool verbose = false) {
Progress progress(n_epochs, verbose);

  const auto n_epochs_per_sample = epochs_per_sample.size();
  double alpha = initial_alpha;

  arma::vec epochs_per_negative_sample(epochs_per_sample / negative_sample_rate);
  arma::vec epoch_of_next_negative_sample(epochs_per_negative_sample);
  arma::vec epoch_of_next_sample(epochs_per_sample);

  SgdWorker<T, DoMove> worker(gradient, positive_head, positive_tail, 
                              epochs_per_sample, epoch_of_next_sample, 
                              epochs_per_negative_sample, 
                              epoch_of_next_negative_sample,
                              head_embedding, tail_embedding, n_vertices, 
                              head_embedding.n_rows, head_embedding.n_cols,
                              seed);

  for (auto n = 0U; n < n_epochs; n++) {
    worker.set_alpha(alpha);
    worker.set_n(n);

    if (parallelize) {
      RcppParallel::parallelFor(0, n_epochs_per_sample, worker, grain_size);
    }
    else {
      worker(0, n_epochs_per_sample);
    }
    alpha = initial_alpha * (1.0 - (double(n) / double(n_epochs)));

    if (Progress::check_abort()) {
      return head_embedding;
    }
    if (verbose) {
      progress.increment();
    }
  }
  return head_embedding;
}

// Reasoning that may come back to haunt me:
// positive_head, positive_tail and epochs_per_sample are read from multiple
// threads: naively, I am hoping that passing by copy as arma vecs should
// prevent R garbage collection from moving this data or causing other issues
// [[Rcpp::export]]
arma::mat optimize_layout_umap(arma::mat head_embedding,
                               Rcpp::Nullable<Rcpp::NumericMatrix> tail_embedding,
                               const arma::uvec positive_head,
                               const arma::uvec positive_tail,
                               unsigned int n_epochs, unsigned int n_vertices,
                               const arma::vec epochs_per_sample,
                               double a, double b,
                               double gamma, double initial_alpha,
                               double negative_sample_rate,
                               unsigned int seed,
                               bool approx_pow,
                               bool parallelize = true,
                               std::size_t grain_size = 1,
                               bool move_other = true,
                               bool verbose = false) {
  
  // For normal UMAP, tail_embedding is NULL and we want to pass 
  // a shallow copy of head_embedding as tail_embedding.
  // When updating new values, tail_embedding is the new coordinate to optimize
  // and gets passed as normal.
  arma::mat tail = tail_embedding.isNull() ? 
  arma::mat(head_embedding.memptr(), head_embedding.n_rows, 
            head_embedding.n_cols, false) :
  Rcpp::as<arma::mat>(tail_embedding);  

  if (approx_pow) {
    const apumap_gradient gradient(a, b, gamma);
    if (move_other) {
      return optimize_layout<apumap_gradient, true>(
          gradient, head_embedding, tail, positive_head, positive_tail, 
          n_epochs, n_vertices, epochs_per_sample, initial_alpha, 
          negative_sample_rate, seed, parallelize, grain_size, verbose);
    }
    else {
      return optimize_layout<apumap_gradient, false>(
          gradient, head_embedding, tail, positive_head, positive_tail, 
          n_epochs, n_vertices, epochs_per_sample, initial_alpha, 
          negative_sample_rate, seed, parallelize, grain_size, verbose);
    }
  }
  else {
    const umap_gradient gradient(a, b, gamma);
    if (move_other) {
      return optimize_layout<umap_gradient, true>(
          gradient, head_embedding, tail, positive_head, positive_tail, 
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, seed, parallelize, grain_size, verbose);
    }
    else {
      return optimize_layout<umap_gradient, false>(
          gradient, head_embedding, tail, positive_head, positive_tail, 
          n_epochs,n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, seed, parallelize, grain_size, verbose);
    }
  }
}

// [[Rcpp::export]]
arma::mat optimize_layout_tumap(arma::mat& head_embedding,
                                Rcpp::Nullable<Rcpp::NumericMatrix> tail_embedding,
                                const arma::uvec positive_head,
                                const arma::uvec positive_tail,
                                unsigned int n_epochs, unsigned int n_vertices,
                                const arma::vec epochs_per_sample,
                                double initial_alpha,
                                double negative_sample_rate,
                                unsigned int seed,
                                bool parallelize = true,
                                std::size_t grain_size = 1,
                                bool move_other = true,
                                bool verbose = false) {
  const tumap_gradient gradient;
  arma::mat tail = tail_embedding.isNull() ? 
  arma::mat(head_embedding.memptr(), head_embedding.n_rows, 
            head_embedding.n_cols, false) :
    Rcpp::as<arma::mat>(tail_embedding);
  if (move_other) {
    return optimize_layout<tumap_gradient, true>(
        gradient, head_embedding, tail, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate, 
        seed, parallelize, grain_size, verbose);
  }
  else {
    return optimize_layout<tumap_gradient, false>(
        gradient, head_embedding, tail, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate, 
        seed, parallelize, grain_size, verbose);
  }
}

// [[Rcpp::export]]
arma::mat optimize_layout_largevis(arma::mat& head_embedding,
                                   Rcpp::Nullable<Rcpp::NumericMatrix> tail_embedding,
                                   const arma::uvec positive_head,
                                   const arma::uvec positive_tail,
                                   unsigned int n_epochs, unsigned int n_vertices,
                                   const arma::vec epochs_per_sample,
                                   double gamma, double initial_alpha,
                                   double negative_sample_rate,
                                   unsigned int seed,
                                   bool parallelize = true,
                                   std::size_t grain_size = 1,
                                   bool move_other = true,
                                   bool verbose = false) {
  const largevis_gradient gradient(gamma);
  arma::mat tail = tail_embedding.isNull() ? 
  arma::mat(head_embedding.memptr(), head_embedding.n_rows, 
            head_embedding.n_cols, false) :
    Rcpp::as<arma::mat>(tail_embedding);
  if (move_other) {
    return optimize_layout<largevis_gradient, true>(
        gradient, head_embedding, tail, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate, 
        seed, parallelize, grain_size, verbose);
  }
  else {
    return optimize_layout<largevis_gradient, false>(
        gradient, head_embedding, tail, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate, 
        seed, parallelize, grain_size, verbose);
  }
}
