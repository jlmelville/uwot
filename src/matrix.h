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

#ifndef UWOT_MATRIX_H
#define UWOT_MATRIX_H

#include <vector>

template <typename T>
void row(const std::vector<double> &m, std::size_t nrow, std::size_t ncol,
         std::size_t r, std::vector<T> &out) {
  for (std::size_t j = 0; j < ncol; j++) {
    out[j] = m[r + j * nrow];
  }
}

template <typename T>
void set_row(std::vector<T> &m, std::size_t nrow, std::size_t ncol,
             std::size_t r, const std::vector<T> &out) {
  for (std::size_t j = 0; j < ncol; j++) {
    m[r + j * nrow] = out[j];
  }
}

#endif // UWOT_MATRIX_H
