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

// Clip a numeric value to within [-4, 4]
double clip(double val) {
  return std::max(std::min(val, 4.0), -4.0);
}

// The squared Euclidean distance between rows at indexes a and b
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



template<typename T>
void optimize_layout(const T& gradient,
                     arma::mat& embedding,
                     const arma::uvec& positive_head,
                     const arma::uvec& positive_tail,
                     int n_epochs, unsigned int n_vertices,
                     const arma::vec& epochs_per_sample,
                     double initial_alpha,
                     double negative_sample_rate,
                     unsigned int seed,
                     bool verbose) {
  Progress progress(n_epochs, verbose);

  const auto ndim = embedding.n_cols;
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

  for (auto n = 0; n < n_epochs; n++) {
    for (arma::uword i = 0; i < n_epochs_per_sample; i++) {
      if (epoch_of_next_sample[i] <= n) {
        arma::uword j = positive_head[i];
        arma::uword k = positive_tail[i];

        const double dist_squared = std::max(rdist(embedding, j, k, ndim), dist_eps);
        const double grad_coeff = gradient.grad_attr(dist_squared);

        for (arma::uword d = 0; d < ndim; d++) {
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

          const double dist_squared = std::max(rdist(embedding, j, k, ndim), dist_eps);
          const double grad_coeff = gradient.grad_rep(dist_squared);

          // This is in the original code, but I strongly suspect this can never happen
          // if (!arma::is_finite(grad_coeff)) {
          //   grad_coeff = 4.0;
          // }

          for (arma::uword d = 0; d < ndim; d++) {
            embedding.at(j, d) +=
              clip(grad_coeff * (embedding.at(j, d) - embedding.at(k, d))) * alpha;
          }
        }
        epoch_of_next_negative_sample[i] += n_neg_samples * epochs_per_negative_sample[i];
      }
    }
    if (Progress::check_abort()) {
      return;
    }

    alpha = initial_alpha * (1.0 - (double(n) / double(n_epochs)));
    if (verbose) {
      progress.increment();
    }
  } // next epoch
}

// [[Rcpp::export]]
void optimize_layout_umap(arma::mat& embedding,
                          const arma::uvec& positive_head,
                          const arma::uvec& positive_tail,
                          int n_epochs, unsigned int n_vertices,
                          const arma::vec& epochs_per_sample,
                          double a, double b,
                          double gamma, double initial_alpha,
                          double negative_sample_rate,
                          unsigned int seed,
                          bool approx_pow,
                          bool verbose) {
  if (approx_pow) {
    const apumap_gradient gradient(a, b, gamma);
    optimize_layout(gradient, embedding, positive_head, positive_tail, n_epochs,
                    n_vertices, epochs_per_sample, initial_alpha,
                    negative_sample_rate, seed, verbose);
  }
  else {
    const umap_gradient gradient(a, b, gamma);
    optimize_layout(gradient, embedding, positive_head, positive_tail, n_epochs,
                    n_vertices, epochs_per_sample, initial_alpha,
                    negative_sample_rate, seed, verbose);
  }

}

// [[Rcpp::export]]
void optimize_layout_tumap(arma::mat& embedding,
                           const arma::uvec& positive_head,
                           const arma::uvec& positive_tail,
                           int n_epochs, unsigned int n_vertices,
                           const arma::vec& epochs_per_sample,
                           double initial_alpha,
                           double negative_sample_rate,
                           unsigned int seed,
                           bool verbose) {
  const tumap_gradient gradient;
  optimize_layout(gradient, embedding, positive_head, positive_tail, n_epochs,
                  n_vertices, epochs_per_sample, initial_alpha,
                  negative_sample_rate, seed, verbose);
}


// [[Rcpp::export]]
void optimize_layout_largevis(arma::mat& embedding,
                          const arma::uvec& positive_head,
                          const arma::uvec& positive_tail,
                          int n_epochs, unsigned int n_vertices,
                          const arma::vec& epochs_per_sample,
                          double gamma, double initial_alpha,
                          double negative_sample_rate,
                          unsigned int seed,
                          bool verbose) {
  const largevis_gradient gradient(gamma);
  optimize_layout(gradient, embedding, positive_head, positive_tail, n_epochs,
                  n_vertices, epochs_per_sample, initial_alpha,
                  negative_sample_rate, seed, verbose);
}
