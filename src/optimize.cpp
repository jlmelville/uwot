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

// Clip a numeric vector so that all values are within [-4, 4]
Rcpp::NumericVector clip(const Rcpp::NumericVector& vec) {
  return Rcpp::pmax(Rcpp::pmin(vec, 4.0), -4.0);
}

// The squared Euclidean distance between vectors a and b
double rdist(const Rcpp::NumericVector& a, const Rcpp::NumericVector& b) {
  return Rcpp::sum((a - b) * (a - b));
}

// [[Rcpp::export]]
void optimize_layout_cpp(Rcpp::NumericMatrix& embedding,
                         const Rcpp::IntegerVector& positive_head,
                         const Rcpp::IntegerVector& positive_tail,
                         int n_epochs, int n_vertices,
                         const Rcpp::NumericVector& epochs_per_sample,
                         double a, double b, double gamma, double initial_alpha,
                         double negative_sample_rate, bool verbose) {
  double alpha = initial_alpha;

  Rcpp::NumericVector epochs_per_negative_sample =
    Rcpp::clone(epochs_per_sample);
  epochs_per_negative_sample =
    epochs_per_negative_sample / negative_sample_rate;

  Rcpp::NumericVector epoch_of_next_negative_sample =
    Rcpp::clone(epochs_per_negative_sample);
  Rcpp::NumericVector epoch_of_next_sample = Rcpp::clone(epochs_per_sample);

  const double dist_eps = std::numeric_limits<double>::epsilon();
  Progress progress(n_epochs, verbose);

  for (int n = 0; n < n_epochs; n++) {
    for (int i = 0; i < epochs_per_sample.size(); i++) {
      if (epoch_of_next_sample[i] <= n) {

        int j = positive_head[i];
        int k = positive_tail[i];
        Rcpp::NumericMatrix::Row current = embedding.row(j);
        Rcpp::NumericMatrix::Row other = embedding.row(k);

        double dist_squared = std::max(rdist(current, other), dist_eps);

        double grad_coeff = -2.0 * a * b * pow(dist_squared, b - 1.0);
        grad_coeff /= a * pow(dist_squared, b) + 1.0;

        Rcpp::NumericVector grad_d = clip(grad_coeff * (current - other));
        current = current + grad_d * alpha;
        other = other - grad_d * alpha;

        epoch_of_next_sample[i] += epochs_per_sample[i];

        int n_neg_samples = int((n - epoch_of_next_negative_sample[i]) /
          epochs_per_negative_sample[i]);

        Rcpp::IntegerVector ks =
          Rcpp::sample(n_vertices, n_neg_samples, true) - 1;
        for (int p = 0; p < ks.size(); p++) {
          int k = ks[p];

          Rcpp::NumericMatrix::Row other_neg = embedding.row(k);

          dist_squared = std::max(rdist(current, other_neg), dist_eps);
          double grad_coeff = (2.0 * gamma * b);
          grad_coeff /= (0.001 + dist_squared) * (a * pow(dist_squared, b) + 1);

          if (!arma::is_finite(grad_coeff)) {
            grad_coeff = 4.0;
          }

          grad_d = clip(grad_coeff * (current - other_neg));
          current = current + grad_d * alpha;
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

