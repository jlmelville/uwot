//  UWOT -- An R package for dimensionality reduction using UMAP
//
//  Copyright (C) 2020 James Melville
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

#include "RcppPerpendicular.h"
#include "uwot/matrix.h"

#if defined(__MINGW32__)
#undef Realloc
#undef Free
#endif

#define __ERROR_PRINTER_OVERRIDE__ REprintf

#include <annoylib.h>
#include <kissrandom.h>

struct UwotAnnoyEuclidean {
  using Distance = Euclidean;
  using S = int32_t;
  using T = float;
};

struct UwotAnnoyCosine {
  using Distance = Angular;
  using S = int32_t;
  using T = float;
};

struct UwotAnnoyManhattan {
  using Distance = Manhattan;
  using S = int32_t;
  using T = float;
};

struct UwotAnnoyHamming {
  using Distance = Hamming;
  using S = int32_t;
  using T = uint64_t;
};

template <typename UwotAnnoyDistance> struct NNWorker {
  const std::string &index_name;
  const std::vector<double> &mat;
  std::size_t nrow;
  std::size_t ncol;
  std::size_t n_neighbors;
  std::size_t search_k;
  std::vector<int> idx;
  std::vector<typename UwotAnnoyDistance::T> dists;

  NNWorker(const std::string &index_name, const std::vector<double> &mat,
           std::size_t ncol, std::size_t n_neighbors, std::size_t search_k)
      : index_name(index_name), mat(mat), nrow(mat.size() / ncol), ncol(ncol),
        n_neighbors(n_neighbors), search_k(search_k),
        idx(nrow * n_neighbors, -1), dists(nrow * n_neighbors) {}

  void operator()(std::size_t begin, std::size_t end) {
    AnnoyIndex<typename UwotAnnoyDistance::S, typename UwotAnnoyDistance::T,
               typename UwotAnnoyDistance::Distance, Kiss64Random>
        index(ncol);
    index.load(index_name.c_str());

    for (auto i = begin; i < end; i++) {
      std::vector<typename UwotAnnoyDistance::T> fv(ncol);
      uwot::get_row(mat, nrow, ncol, i, fv);

      std::vector<typename UwotAnnoyDistance::S> result;
      std::vector<typename UwotAnnoyDistance::T> distances;

      index.get_nns_by_vector(fv.data(), n_neighbors, search_k, &result,
                              &distances);
      if (result.size() != n_neighbors || distances.size() != n_neighbors) {
        break;
      }

      for (std::size_t j = 0; j < n_neighbors; j++) {
        uwot::set_row(dists, nrow, n_neighbors, i, distances);
        uwot::set_row(idx, nrow, n_neighbors, i, result);
      }
    }
  }
};
