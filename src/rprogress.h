//  UWOT -- An R package for dimensionality reduction using UMAP
//
//  Copyright (C) 2021 James Melville
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
//  along with UWOT.  If not, see <http://www.gnu.org/licenses

#ifndef UWOT_RPROGRESS_H
#define UWOT_RPROGRESS_H

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

struct RProgress {

  Progress progress;
  bool verbose;

  RProgress(std::size_t n_epochs, bool verbose)
      : progress(n_epochs, verbose), verbose(verbose) {}

  auto is_aborted() -> bool {
    bool aborted = Progress::check_abort();
    if (aborted) {
      progress.cleanup();
    }
    return aborted;
  }

  void report() {
    if (verbose) {
      progress.increment();
    }
  }
};

#endif // UWOT_RPROGRESS_H