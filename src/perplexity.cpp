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

#include <Rcpp.h>

#include "RcppPerpendicular.h"
#include "uwot/perplexity.h"

using namespace Rcpp;

// [[Rcpp::export]]
List calc_row_probabilities_parallel(NumericVector nn_dist,
                                     std::size_t n_vertices, double perplexity,
                                     std::size_t n_iter = 200,
                                     double tol = 1e-5, bool ret_sigma = false,
                                     std::size_t n_threads = 0,
                                     std::size_t grain_size = 1) {
  std::size_t n_neighbors = nn_dist.size() / n_vertices;
  auto nn_distv = as<std::vector<double>>(nn_dist);

  double target = std::log(perplexity);
  std::atomic_size_t n_search_fails{0};
  std::vector<double> nn_weights(n_vertices * n_neighbors);
  std::vector<double> sigmas(ret_sigma ? n_vertices : 0);

  auto worker = [&](std::size_t begin, std::size_t end) {
    uwot::perplexity_search(begin, end, nn_distv, n_neighbors, target, tol,
                            n_iter, nn_weights, ret_sigma, sigmas,
                            n_search_fails);
  };

  RcppPerpendicular::parallel_for(0, n_vertices, worker, n_threads, grain_size);

  auto res = List::create(
      _("matrix") = NumericMatrix(n_neighbors, n_vertices, nn_weights.begin()),
      _("n_failures") = static_cast<std::size_t>(n_search_fails));
  if (ret_sigma) {
    res["sigma"] = sigmas;
  }
  return res;
}
