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

#include <Rcpp.h>

#include "uwot/connected_components.h"

using namespace Rcpp;

// [[Rcpp::export]]
List connected_components_undirected(std::size_t N,
                                     const std::vector<int> &indices1,
                                     const std::vector<int> &indptr1,
                                     const std::vector<int> &indices2,
                                     const std::vector<int> &indptr2) {

  std::pair<unsigned int, std::vector<int>> result =
      uwot::connected_components_undirected(N, indices1, indptr1, indices2,
                                            indptr2);

  return List::create(_["n_components"] = result.first,
                      _["labels"] = result.second);
}
