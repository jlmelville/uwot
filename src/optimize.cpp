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
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "gradient.h"

// Clip a numeric value to within [-4, 4]
double clip(double val) {
  return std::max(std::min(val, 4.0), -4.0);
}

// The squared Euclidean distance between vectors a and b
double rdist(const arma::rowvec& a,
             const arma::rowvec& b,
             const arma::uword ndim) {
  double sum = 0.0;
  for (arma::uword i = 0; i < ndim; i++) {
    sum += (a[i] - b[i]) * (a[i] - b[i]);
  }

  return sum;
}


template<typename T>
void optimize_layout(const T& gradient,
                     arma::mat& embedding,
                     const arma::uvec& positive_head,
                     const arma::uvec& positive_tail,
                     int n_epochs, int n_vertices,
                     const arma::vec& epochs_per_sample,
                     double initial_alpha,
                     double negative_sample_rate, bool verbose) {
  Progress progress(n_epochs, verbose);

  const arma::uword ndim = embedding.n_cols;
  const arma::uword n_epochs_per_sample = epochs_per_sample.size();

  const double dist_eps = std::numeric_limits<double>::epsilon();
  double alpha = initial_alpha;

  arma::vec epochs_per_negative_sample(epochs_per_sample / negative_sample_rate);
  arma::vec epoch_of_next_negative_sample(epochs_per_negative_sample);
  arma::vec epoch_of_next_sample(epochs_per_sample);

  for (int n = 0; n < n_epochs; n++) {
    for (arma::uword i = 0; i < n_epochs_per_sample; i++) {
      if (epoch_of_next_sample[i] <= n) {
        arma::uword j = positive_head[i];
        arma::uword k = positive_tail[i];

        arma::subview_row<double> current = embedding.row(j);
        arma::subview_row<double> other = embedding.row(k);
        const double dist_squared = std::max(rdist(current, other, ndim), dist_eps);

        const double grad_coeff = gradient.grad_attr(dist_squared);

        for (arma::uword d = 0; d < ndim; d++) {
          double grad_d = clip(grad_coeff * (current[d] - other[d])) * alpha;
          current[d] += grad_d;
          other[d] -= grad_d;
        }

        epoch_of_next_sample[i] += epochs_per_sample[i];

        int n_neg_samples = int((n - epoch_of_next_negative_sample[i]) /
                                epochs_per_negative_sample[i]);

        Rcpp::IntegerVector ks = Rcpp::sample(n_vertices, n_neg_samples, true) - 1;
        for (int p = 0; p < ks.size(); p++) {
          arma::uword k = ks[p];

          if (j == k) {
            continue;
          }

          arma::subview_row<double> other_neg = embedding.row(k);

          const double dist_squared = std::max(rdist(current, other_neg, ndim), dist_eps);
          const double grad_coeff = gradient.grad_rep(dist_squared);

          // This is in the original code, but I strongly suspect this can never happen
          // if (!arma::is_finite(grad_coeff)) {
          //   grad_coeff = 4.0;
          // }

          for (arma::uword d = 0; d < ndim; d++) {
            current[d] += clip(grad_coeff * (current[d] - other_neg[d])) * alpha;
          }
        }

        epoch_of_next_negative_sample[i] +=
          n_neg_samples * epochs_per_negative_sample[i];
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
                          int n_epochs, int n_vertices,
                          const arma::vec& epochs_per_sample,
                          double a, double b,
                          double gamma, double initial_alpha,
                          double negative_sample_rate, bool verbose) {
  const umap_gradient gradient(a, b, gamma);
  optimize_layout(gradient, embedding, positive_head, positive_tail, n_epochs,
                  n_vertices, epochs_per_sample, initial_alpha,
                  negative_sample_rate, verbose);
}

// [[Rcpp::export]]
void optimize_layout_tumap(arma::mat& embedding,
                           const arma::uvec& positive_head,
                           const arma::uvec& positive_tail,
                           int n_epochs, int n_vertices,
                           const arma::vec& epochs_per_sample,
                           double initial_alpha,
                           double negative_sample_rate, bool verbose) {
  const tumap_gradient gradient;
  optimize_layout(gradient, embedding, positive_head, positive_tail, n_epochs,
                  n_vertices, epochs_per_sample, initial_alpha,
                  negative_sample_rate, verbose);
}
