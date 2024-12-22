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

#include <cmath>
#include <vector>

#include "uwot/smooth_knn.h"

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List smooth_knn_distances_parallel(
    NumericVector nn_dist, IntegerVector nn_ptr, bool skip_first,
    NumericVector target, std::size_t n_iter = 64,
    double local_connectivity = 1.0, double tol = 1e-5,
    double min_k_dist_scale = 1e-3, bool ret_sigma = false,
    std::size_t n_threads = 0, std::size_t grain_size = 1) {

  std::size_t n_neighbors = 0;
  std::size_t n_vertices = 0;
  std::vector<std::size_t> nn_ptrv(0);
  if (nn_ptr.size() == 0) {
    stop("nn_ptr cannot be empty");
  }
  if (nn_ptr.size() == 1) {
    // Size optimization for the typical kNN graph case:
    // all points have the same number of neighbors so just store that number
    // as the single entry in the nn_ptr vector
    n_neighbors = nn_ptr[0];
    if (nn_dist.size() % n_neighbors != 0) {
      stop("Invalid n_neighbors for nn_dist size");
    }
    nn_ptrv = std::vector<std::size_t>{n_neighbors};
    n_vertices = nn_dist.size() / n_neighbors;
  } else {
    nn_ptrv = as<std::vector<std::size_t>>(nn_ptr);
    n_vertices = nn_ptrv.size() - 1;
  }

  auto targetv = as<std::vector<double>>(target);

  auto nn_distv = as<std::vector<double>>(nn_dist);
  double mean_distances = uwot::mean_average(nn_distv);

  std::atomic_size_t n_search_fails{0};
  std::vector<double> nn_weights(nn_dist.size());
  std::vector<double> sigmas(ret_sigma ? n_vertices : 0);
  std::vector<double> rhos(ret_sigma ? n_vertices : 0);

  auto worker = [&](std::size_t begin, std::size_t end) {
    uwot::smooth_knn(begin, end, nn_distv, nn_ptrv, skip_first, targetv,
                     local_connectivity, tol, n_iter, min_k_dist_scale,
                     mean_distances, ret_sigma, nn_weights, sigmas, rhos,
                     n_search_fails);
  };

  RcppPerpendicular::parallel_for(n_vertices, worker, n_threads, grain_size);

  auto res = List::create(
      _("matrix") = NumericVector(nn_weights.begin(), nn_weights.end()),
      _("n_failures") = static_cast<std::size_t>(n_search_fails));
  if (ret_sigma) {
    res["sigma"] = sigmas;
    res["rho"] = rhos;
  }
  return res;
}

// [[Rcpp::export]]
List reset_local_metrics_parallel(IntegerVector indptr,
                                  NumericVector probabilities,
                                  std::size_t n_iter = 32, double tol = 1e-5,
                                  double num_local_metric_neighbors = 15.0,
                                  std::size_t n_threads = 0) {

  auto n_vertices = indptr.size() - 1;
  double target = std::log2(num_local_metric_neighbors);
  std::atomic_size_t n_search_fails{0};
  auto prob_ptrv = as<std::vector<std::size_t>>(indptr);
  auto probabilitiesv = as<std::vector<double>>(probabilities);
  auto worker = [&](std::size_t begin, std::size_t end) {
    uwot::reset_local_metric(begin, end, probabilitiesv, prob_ptrv, target, tol,
                             n_iter, n_search_fails);
  };
  RcppPerpendicular::parallel_for(n_vertices, worker, n_threads);

  auto res = List::create(
      _("values") = NumericVector(probabilitiesv.begin(), probabilitiesv.end()),
      _("n_failures") = static_cast<std::size_t>(n_search_fails));

  return res;
}
