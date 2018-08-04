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
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "gradient.h"
#include "tauprng.h"

// Clip a numeric value to within [-clip_max, clip_max]
double clip(const double val, const double clip_max) {
  return std::max(std::min(val, clip_max), -clip_max);
}

// The squared Euclidean distance between rows at indexes a and b
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



template<typename T>
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
                     bool move_other,
                     bool verbose) {
  Progress progress(n_epochs, verbose);

  const auto ndim = head_embedding.n_cols;
  const auto n_epochs_per_sample = epochs_per_sample.size();

  const double dist_eps = std::numeric_limits<double>::epsilon();
  double alpha = initial_alpha;

  arma::vec epochs_per_negative_sample(epochs_per_sample / negative_sample_rate);
  arma::vec epoch_of_next_negative_sample(epochs_per_negative_sample);
  arma::vec epoch_of_next_sample(epochs_per_sample);

  // Reproducibly initialize the state of the Tausworthe PRNG.
  // seed is an unsigned integer sampled from R, and hence is reproducible.
  // (It is dependent on, but NOT the same as, any value that might be passed
  // to set.seed)
  // Use that seed to initialize the standard library Mersenne Twister PRNG.
  // Generate three random integers. Use those to initialize the Tausworthe
  // PRNG.
  std::mt19937 rng(seed);
  std::uniform_int_distribution<long> gen(-2147483647, 2147483646);
  long s1 = gen(rng);
  long s2 = gen(rng);
  long s3 = gen(rng);
  tau_prng prng(s1, s2, s3);

  for (auto n = 0U; n < n_epochs; n++) {
    for (arma::uword i = 0; i < n_epochs_per_sample; i++) {
      if (epoch_of_next_sample[i] <= n) {
        arma::uword j = positive_head[i];
        arma::uword k = positive_tail[i];

        const double dist_squared = std::max(rdist(head_embedding, tail_embedding, j, k, ndim), dist_eps);
        const double grad_coeff = gradient.grad_attr(dist_squared);

        for (arma::uword d = 0; d < ndim; d++) {
          double grad_d = clip(grad_coeff * (head_embedding.at(j, d) - tail_embedding.at(k, d)), gradient.clip_max) * alpha;
          head_embedding.at(j, d) += grad_d;
          if (move_other) {
            tail_embedding.at(k, d) -= grad_d;
          }
        }

        epoch_of_next_sample[i] += epochs_per_sample[i];

        unsigned int n_neg_samples = static_cast<unsigned int>((n - epoch_of_next_negative_sample[i]) /
                                epochs_per_negative_sample[i]);

        for (unsigned int p = 0; p < n_neg_samples; p++) {
          arma::uword k = prng() % n_vertices;
          if (j == k) {
            continue;
          }

          const double dist_squared = std::max(rdist(head_embedding, tail_embedding, j, k, ndim), dist_eps);
          const double grad_coeff = gradient.grad_rep(dist_squared);

          // This is in the original code, but I strongly suspect this can never happen
          // if (!arma::is_finite(grad_coeff)) {
          //   grad_coeff = 4.0;
          // }

          for (arma::uword d = 0; d < ndim; d++) {
            head_embedding.at(j, d) +=
              clip(grad_coeff * (head_embedding.at(j, d) - tail_embedding.at(k, d)), gradient.clip_max) * alpha;
          }
        }
        epoch_of_next_negative_sample[i] += n_neg_samples * epochs_per_negative_sample[i];
      }
    }
    if (Progress::check_abort()) {
      return head_embedding;
    }

    alpha = initial_alpha * (1.0 - (double(n) / double(n_epochs)));
    if (verbose) {
      progress.increment();
    }
  } // next epoch
  return head_embedding;
}

// [[Rcpp::export]]
arma::mat optimize_layout_umap(arma::mat& head_embedding,
                               arma::mat& tail_embedding,
                          const arma::uvec& positive_head,
                          const arma::uvec& positive_tail,
                          unsigned int n_epochs, unsigned int n_vertices,
                          const arma::vec& epochs_per_sample,
                          double a, double b,
                          double gamma, double initial_alpha,
                          double negative_sample_rate,
                          unsigned int seed,
                          bool approx_pow,
                          bool move_other,
                          bool verbose) {
  if (approx_pow) {
    const apumap_gradient gradient(a, b, gamma);
    return optimize_layout(gradient, head_embedding, tail_embedding, positive_head, positive_tail, n_epochs,
                    n_vertices, epochs_per_sample, initial_alpha,
                    negative_sample_rate, seed, move_other, verbose);
  }
  else {
    const umap_gradient gradient(a, b, gamma);
    return optimize_layout(gradient, head_embedding, tail_embedding, positive_head, positive_tail, n_epochs,
                    n_vertices, epochs_per_sample, initial_alpha,
                    negative_sample_rate, seed, move_other, verbose);
  }

}

// [[Rcpp::export]]
arma::mat optimize_layout_tumap(arma::mat& head_embedding,
                                arma::mat& tail_embedding,
                           const arma::uvec& positive_head,
                           const arma::uvec& positive_tail,
                           unsigned int n_epochs, unsigned int n_vertices,
                           const arma::vec& epochs_per_sample,
                           double initial_alpha,
                           double negative_sample_rate,
                           unsigned int seed,
                           bool move_other,
                           bool verbose) {
  const tumap_gradient gradient;
  return optimize_layout(gradient, head_embedding, tail_embedding, positive_head, positive_tail, n_epochs,
                  n_vertices, epochs_per_sample, initial_alpha,
                  negative_sample_rate, seed, move_other, verbose);
}


// [[Rcpp::export]]
arma::mat optimize_layout_largevis(arma::mat& head_embedding,
                                   arma::mat& tail_embedding,
                          const arma::uvec& positive_head,
                          const arma::uvec& positive_tail,
                          unsigned int n_epochs, unsigned int n_vertices,
                          const arma::vec& epochs_per_sample,
                          double gamma, double initial_alpha,
                          double negative_sample_rate,
                          unsigned int seed,
                          bool move_other,
                          bool verbose) {
  const largevis_gradient gradient(gamma);

  return optimize_layout(gradient, head_embedding, tail_embedding, positive_head, positive_tail, n_epochs,
                  n_vertices, epochs_per_sample, initial_alpha,
                  negative_sample_rate, seed, move_other, verbose);
}
