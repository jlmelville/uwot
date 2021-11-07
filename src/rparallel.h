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

#ifndef UWOT_RPARALLEL_H
#define UWOT_RPARALLEL_H

#include "RcppPerpendicular.h"

struct RParallel {
  std::size_t n_threads;
  std::size_t grain_size;

  RParallel(std::size_t n_threads, std::size_t grain_size)
      : n_threads(n_threads), grain_size(grain_size) {}

  template <typename Worker> void pfor(std::size_t n_items, Worker &worker) {
    RcppPerpendicular::pfor(n_items, worker, n_threads, grain_size);
  }

  template <typename Worker>
  void pfor(std::size_t begin, std::size_t end, Worker &worker) {
    RcppPerpendicular::pfor(begin, end, worker, n_threads, grain_size);
  }
};

struct RSerial {
  template <typename Worker> void pfor(std::size_t n_items, Worker &worker) {
    pfor(0, n_items, worker);
  }

  template <typename Worker>
  void pfor(std::size_t begin, std::size_t end, Worker &worker) {
    worker(begin, end, 0);
  }
};

#endif // UWOT_RPARALLEL_H
