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

#include <vector>

#include "uwot/smooth_knn.h"

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List smooth_knn_distances_parallel(
    NumericVector nn_dist, std::size_t n_vertices, double target = -1.0,
    std::size_t n_iter = 64, double local_connectivity = 1.0,
    double bandwidth = 1.0, double tol = 1e-5, double min_k_dist_scale = 1e-3,
    bool ret_sigma = false, std::size_t n_threads = 0,
    std::size_t grain_size = 1) {

  std::size_t n_neighbors = nn_dist.size() / n_vertices;

  if (target < 0.0) {
    target = std::log2(n_neighbors);
  }

  auto nn_distv = as<std::vector<double>>(nn_dist);
  double mean_distances = uwot::mean_average(nn_distv);

  std::atomic_size_t n_search_fails{0};
  std::vector<double> nn_weights(n_vertices * n_neighbors);
  std::vector<double> sigmas(ret_sigma ? n_vertices : 0);
  std::vector<double> rhos(ret_sigma ? n_vertices : 0);

  auto worker = [&](std::size_t begin, std::size_t end) {
    uwot::smooth_knn(begin, end, nn_distv, n_neighbors, target,
                     local_connectivity, tol, n_iter, bandwidth,
                     min_k_dist_scale, mean_distances, ret_sigma, nn_weights,
                     sigmas, rhos, n_search_fails);
  };

  RcppPerpendicular::parallel_for(0, n_vertices, worker, n_threads, grain_size);

  auto res = List::create(
      _("matrix") = NumericMatrix(n_neighbors, n_vertices, nn_weights.begin()),
      _("n_failures") = static_cast<std::size_t>(n_search_fails));
  if (ret_sigma) {
    res["sigma"] = sigmas;
    res["rho"] = rhos;
  }
  return res;
}
