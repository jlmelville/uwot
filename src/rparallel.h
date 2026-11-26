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

#include "pforr.h"
#include "uwot/epoch.h"
#include <Rcpp.h>

struct RParallel {
  std::size_t n_threads;
  std::size_t grain_size;

  RParallel(std::size_t n_threads, std::size_t grain_size)
      : n_threads(n_threads), grain_size(grain_size) {}

  template <typename Worker> void pfor(std::size_t n_items, Worker &worker) {
    pforr::parallel_for_indexed(n_items, worker, n_threads, grain_size);
  }

  template <typename G, typename U, typename F>
  void pfor(std::size_t n, uwot::EdgeWorker<G, U, F> &worker) {
    visit_windows(n, worker);
  }

  template <typename G, typename U, typename F>
  void pfor(std::size_t n, uwot::NodeWorker<G, U, F> &worker) {
    visit_windows(n, worker);
  }

  template <typename Worker>
  void pfor(std::size_t begin, std::size_t end, Worker &worker) {
    pforr::parallel_for_indexed(begin, end, worker, n_threads, grain_size);
  }

private:
  template <typename Worker> void visit_windows(std::size_t n, Worker &worker) {
    constexpr std::size_t edges_per_poll = 65536;
    const auto n_edges = worker.positive_head.size();
    if (n_edges <= edges_per_poll || n == 0) {
      pforr::parallel_for_indexed(n, worker, n_threads, grain_size);
      return;
    }
    // Preserve each worker's original partition and random stream. Only the
    // main thread polls R, after all workers have finished the current slice.
    const auto ranges = pforr::split_input_range({0, n}, n_threads, grain_size);
    using Window = decltype(worker.make_window(0, 0));
    std::vector<Window> windows;
    windows.reserve(ranges.size());
    for (std::size_t i = 0; i < ranges.size(); ++i) {
      windows.emplace_back(worker.make_window(ranges[i].second, i));
    }
    const auto degree = (std::max)(std::size_t{1}, n_edges / n);
    const auto stride = (std::max)(std::size_t{1}, edges_per_poll / degree);
    for (std::size_t offset = 0; offset < ranges.front().second;
         offset += stride) {
      auto visit = [&](std::size_t begin, std::size_t end, std::size_t) {
        for (auto i = begin; i < end; ++i) {
          const auto first = ranges[i].first + offset;
          const auto last = ranges[i].second;
          if (first < last)
            windows[i](first, first + (std::min)(stride, last - first));
        }
      };
      pforr::parallel_for_indexed(ranges.size(), visit, n_threads);
      Rcpp::checkUserInterrupt();
    }
  }
};

#endif // UWOT_RPARALLEL_H
