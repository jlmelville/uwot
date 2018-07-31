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
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fast_intersection_cpp(const IntegerVector rows,
                                    const IntegerVector cols,
                                    NumericVector values,
                                    const IntegerVector target,
                                    double unknown_dist = 1.0,
                                    double far_dist = 5.0) {
    double ex_unknown = std::exp(-unknown_dist);
    double ex_far = std::exp(-far_dist);

    auto len = values.length();

    for (auto nz = 0; nz < len; ++nz) {
      auto i = rows[nz];
      auto j = cols[nz];
      if (IntegerVector::is_na(target[i]) || IntegerVector::is_na(target[j])) {
        values[nz] = values[nz] * ex_unknown;
      }
      else if (target[i] != target[j]) {
        values[nz] = values[nz] * ex_far;
      }
    }

    return values;
}
