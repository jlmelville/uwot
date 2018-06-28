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
#include "tthread/fast_mutex.h"

template <typename T>
struct SgdWorker : public RcppParallel::Worker {

  int n; // epoch counter
  double alpha;
  const T gradient;
  const arma::uvec positive_head;
  const arma::uvec positive_tail;
  const arma::vec epochs_per_sample;

  arma::vec epoch_of_next_sample;
  const arma::vec epochs_per_negative_sample;
  arma::vec epoch_of_next_negative_sample;
  arma::mat& embedding;
  unsigned int n_vertices;
  const arma::uword nrow;
  const arma::uword ncol;
  tthread::fast_mutex mutex;
  std::mt19937 rng;
  std::uniform_int_distribution<long> gen;
  const double dist_eps;

  SgdWorker(
    const T& gradient,
    const arma::uvec& positive_head,
    const arma::uvec& positive_tail,
    const arma::vec& epochs_per_sample,

    arma::vec& epoch_of_next_sample,
    const arma::vec& epochs_per_negative_sample,
    arma::vec& epoch_of_next_negative_sample,

    arma::mat& embedding,
    unsigned int n_vertices,
    const arma::uword nrow,
    const arma::uword ncol,
    unsigned int seed) :

    n(0), alpha(0.0), gradient(gradient), positive_head(positive_head), positive_tail(positive_tail),
    epochs_per_sample(epochs_per_sample),
    epoch_of_next_sample(epoch_of_next_sample),
    epochs_per_negative_sample(epochs_per_negative_sample),
    epoch_of_next_negative_sample(epoch_of_next_negative_sample),
    embedding(embedding), n_vertices(n_vertices), nrow(nrow), ncol(ncol), rng(seed),
    gen(-2147483647, 2147483646),
    dist_eps(std::numeric_limits<double>::epsilon()) {  }

  void operator()(std::size_t begin, std::size_t end) {
    // Each window gets its own fast PRNG state, so no locking needed inside the loop.
    // Want separate seeds though, so seed the fast PRNG with three random numbers
    // taken from the mt19937 generator, which is shared across windows, so is locked.
    // Probably this is a bit of a waste of time:
    // Could use the mt19937 seed, the begin and the epoch number as seeds?
    // Doesn't waste much time, though.
    long s1, s2, s3;
    {
      tthread::lock_guard<tthread::fast_mutex> guard(mutex);
      s1 = gen(rng);
      s2 = gen(rng); // technically this needs to always be > 7
      s3 = gen(rng); // should be > 15
    }
    tau_prng prng(s1, s2, s3);

    for (std::size_t i = begin; i < end; i++) {
      if (epoch_of_next_sample[i] <= n) {
        arma::uword j = positive_head[i];
        arma::uword k = positive_tail[i];

        const double dist_squared = std::max(rdist(embedding, j, k, ncol), dist_eps);
        const double grad_coeff = gradient.grad_attr(dist_squared);

        for (arma::uword d = 0; d < ncol; d++) {
          double grad_d = clip(grad_coeff * (embedding.at(j, d) - embedding.at(k, d))) * alpha;
          embedding.at(j, d) += grad_d;
          embedding.at(k, d) -= grad_d;
        }

        epoch_of_next_sample[i] += epochs_per_sample[i];

        unsigned int n_neg_samples = static_cast<unsigned int>((n - epoch_of_next_negative_sample[i]) /
                                                               epochs_per_negative_sample[i]);

        for (unsigned int p = 0; p < n_neg_samples; p++) {
          arma::uword k = prng() % n_vertices;

          if (j == k) {
            continue;
          }

          const double dist_squared = std::max(rdist(embedding, j, k, ncol), dist_eps);
          const double grad_coeff = gradient.grad_rep(dist_squared);

          for (arma::uword d = 0; d < ncol; d++) {
            embedding.at(j, d) +=
              clip(grad_coeff * (embedding.at(j, d) - embedding.at(k, d))) * alpha;
          }
        }
        epoch_of_next_negative_sample[i] += n_neg_samples * epochs_per_negative_sample[i];
      }
    }
  }

  void set_n(int n) {
    this->n = n;
  }

  void set_alpha(double alpha) {
    this->alpha = alpha;
  }

  double clip(double val) {
    return std::max(std::min(val, 4.0), -4.0);
  }

  double rdist(const arma::mat& mat,
               const arma::uword a,
               const arma::uword b,
               const arma::uword ndim) {
    double sum = 0.0;
    for (arma::uword i = 0; i < ndim; i++) {
      sum += (mat.at(a, i) - mat.at(b, i)) * (mat.at(a, i) - mat.at(b, i));
    }

    return sum;
  }

};

template<typename T>
void optimize_layout_parallel(const T& gradient,
                              arma::mat& embedding,
                              const arma::uvec& positive_head,
                              const arma::uvec& positive_tail,
                              int n_epochs, unsigned int n_vertices,
                              const arma::vec& epochs_per_sample,
                              double initial_alpha,
                              double negative_sample_rate,
                              unsigned int seed,
                              std::size_t grain_size = 1000,
                              bool verbose = false) {
  Progress progress(n_epochs, verbose);

  const auto n_epochs_per_sample = epochs_per_sample.size();
  double alpha = initial_alpha;

  arma::vec epochs_per_negative_sample(epochs_per_sample / negative_sample_rate);
  arma::vec epoch_of_next_negative_sample(epochs_per_negative_sample);
  arma::vec epoch_of_next_sample(epochs_per_sample);

  SgdWorker<T> worker(gradient, positive_head, positive_tail, epochs_per_sample,
                      epoch_of_next_sample, epochs_per_negative_sample,
                      epoch_of_next_negative_sample, embedding,
                      n_vertices, embedding.n_rows, embedding.n_cols, seed);
  for (auto n = 0; n < n_epochs; n++) {
    worker.set_alpha(alpha);
    worker.set_n(n);

    RcppParallel::parallelFor(0, n_epochs_per_sample, worker, grain_size);

    alpha = initial_alpha * (1.0 - (double(n) / double(n_epochs)));

    if (Progress::check_abort()) {
      return;
    }
    if (verbose) {
      progress.increment();
    }
  }
}

// [[Rcpp::export]]
void optimize_layout_umap_parallel(arma::mat& embedding,
                          const arma::uvec& positive_head,
                          const arma::uvec& positive_tail,
                          int n_epochs, unsigned int n_vertices,
                          const arma::vec& epochs_per_sample,
                          double a, double b,
                          double gamma, double initial_alpha,
                          double negative_sample_rate,
                          unsigned int seed,
                          bool approx_pow,
                          std::size_t grain_size = 1000,
                          bool verbose = false) {
  if (approx_pow) {
    const apumap_gradient gradient(a, b, gamma);
    optimize_layout_parallel(gradient, embedding, positive_head, positive_tail, n_epochs,
                             n_vertices, epochs_per_sample, initial_alpha,
                             negative_sample_rate, seed, grain_size, verbose);
  }
  else {
    const umap_gradient gradient(a, b, gamma);
    optimize_layout_parallel(gradient, embedding, positive_head, positive_tail, n_epochs,
                             n_vertices, epochs_per_sample, initial_alpha,
                             negative_sample_rate, seed, grain_size, verbose);
  }
}

// [[Rcpp::export]]
void optimize_layout_tumap_parallel(arma::mat& embedding,
                           const arma::uvec& positive_head,
                           const arma::uvec& positive_tail,
                           int n_epochs, unsigned int n_vertices,
                           const arma::vec& epochs_per_sample,
                           double initial_alpha,
                           double negative_sample_rate,
                           unsigned int seed,
                           std::size_t grain_size = 1000,
                           bool verbose = false) {
  const tumap_gradient gradient;
  optimize_layout_parallel(gradient, embedding, positive_head, positive_tail, n_epochs,
                  n_vertices, epochs_per_sample, initial_alpha,
                  negative_sample_rate, seed, grain_size, verbose);
}

// [[Rcpp::export]]
void optimize_layout_largevis_parallel(arma::mat& embedding,
                              const arma::uvec& positive_head,
                              const arma::uvec& positive_tail,
                              int n_epochs, unsigned int n_vertices,
                              const arma::vec& epochs_per_sample,
                              double gamma, double initial_alpha,
                              double negative_sample_rate,
                              unsigned int seed,
                              std::size_t grain_size = 1000,
                              bool verbose = false) {
  const largevis_gradient gradient(gamma);
  optimize_layout_parallel(gradient, embedding, positive_head, positive_tail, n_epochs,
                  n_vertices, epochs_per_sample, initial_alpha,
                  negative_sample_rate, seed, grain_size, verbose);
}
