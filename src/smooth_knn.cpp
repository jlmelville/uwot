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

#include "uwot/matrix.h"
#include "uwot/smooth_knn.h"

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List smooth_knn_distances_parallel(NumericMatrix nn_dist,
                                   std::size_t n_iter = 64,
                                   double local_connectivity = 1.0,
                                   double bandwidth = 1.0, double tol = 1e-5,
                                   double min_k_dist_scale = 1e-3,
                                   std::size_t n_threads = 0,
                                   std::size_t grain_size = 1) {

  std::size_t n_vertices = nn_dist.nrow();
  std::size_t n_neighbors = nn_dist.ncol();

  auto nn_distv = as<std::vector<double>>(nn_dist);
  uwot::SmoothKnnWorker worker(nn_distv, n_vertices, n_iter, local_connectivity,
                               bandwidth, tol, min_k_dist_scale);

  if (n_threads > 0) {
    RcppPerpendicular::parallel_for(0, n_vertices, worker, n_threads,
                                    grain_size);
  } else {
    worker(0, n_vertices);
  }

  return List::create(
      _("matrix") =
          NumericMatrix(n_vertices, n_neighbors, worker.nn_weights.begin()),
      _("n_failures") = static_cast<std::size_t>(worker.n_search_fails));
}
