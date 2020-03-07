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
List calc_row_probabilities_parallel(NumericMatrix nn_dist,
                                     IntegerMatrix nn_idx, double perplexity,
                                     std::size_t n_iter = 200,
                                     double tol = 1e-5,
                                     std::size_t n_threads = 0,
                                     std::size_t grain_size = 1) {

  std::size_t n_vertices = nn_dist.nrow();
  std::size_t n_neighbors = nn_dist.ncol();

  auto nn_distv = as<std::vector<double>>(nn_dist);
  auto nn_idxv = as<std::vector<int>>(nn_idx);

  NumericMatrix res(n_vertices, n_neighbors);
  uwot::PerplexityWorker worker(nn_distv, nn_idxv, n_vertices, perplexity,
                                n_iter, tol);

  if (n_threads > 0) {
    RcppPerpendicular::parallel_for(0, n_vertices, worker, n_threads,
                                    grain_size);
  } else {
    worker(0, n_vertices);
  }

  return List::create(
      _("matrix") = NumericMatrix(n_vertices, n_neighbors, worker.res.begin()),
      _("n_failures") = static_cast<std::size_t>(worker.n_search_fails));
}
