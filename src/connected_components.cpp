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

// Translated from the Python source code of:
//   scipy.sparse.csgraph.connected_components
//   Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
//   License: BSD, (C) 2012
// You may also use the (non-Rcpp) parts of the C++ algorithm under the same 
// license.

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List connected_components_undirected(
    const unsigned long N,
    const IntegerVector indices1,
    const IntegerVector indptr1,
    const IntegerVector indices2,
    const IntegerVector indptr2
)
{
  const int VOID = -1;
  const int END = -2;
  std::vector<int> labels(N, VOID);
  std::vector<int> SS(labels);
  int label = 0;
  int SS_head = END;
  for (unsigned int v = 0; v < N; ++v) {
    unsigned int vv = v;
    if (labels[vv] == VOID) {
      SS_head = vv;
      SS[vv] = END;
      while (SS_head != END) {
        vv = SS_head;
        SS_head = SS[vv];
        labels[vv] = label;
        for (int jj = indptr1[vv]; jj < indptr1[vv + 1]; ++jj) {
          int ww = indices1[jj];
          if (SS[ww] == VOID) {
            SS[ww] = SS_head;
            SS_head = ww;
          }
        }
        for (int jj = indptr2[vv]; jj < indptr2[vv + 1]; ++jj) {
          int ww = indices2[jj];
          if (SS[ww] == VOID) {
            SS[ww] = SS_head;
            SS_head = ww;
          }
        }
      }
      ++label;
    }
  }
  return List::create(
    _["n_components"] = label,
    _["labels"] = labels
  );
}