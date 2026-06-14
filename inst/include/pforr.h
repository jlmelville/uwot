// Taken from RcppParallel.h and then modified slightly to rename header guards
// and namespaces to avoid any potential clashes. RcppParallel is licensed under
// GPLv2 or later:

// pfor.h a version of parallel for based on RcppParallel
// Copyright (C) 2020 James Melville
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
// USA.

#ifndef PFORR
#define PFORR

#include <algorithm>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace pforr {

using IndexRange = std::pair<std::size_t, std::size_t>;

class ThreadJoiner {
public:
  explicit ThreadJoiner(std::vector<std::thread> &threads) : threads(threads) {}

  ~ThreadJoiner() {
    for (auto &thread : threads) {
      if (thread.joinable()) {
        thread.join();
      }
    }
  }

  ThreadJoiner(const ThreadJoiner &) = delete;
  auto operator=(const ThreadJoiner &) -> ThreadJoiner & = delete;

private:
  std::vector<std::thread> &threads;
};

template <typename Worker>
auto worker_thread(Worker &worker, const IndexRange &range,
                   std::exception_ptr &worker_exception,
                   std::mutex &worker_exception_mutex) -> void {
  try {
    worker(range.first, range.second);
  } catch (...) {
    std::lock_guard<std::mutex> guard(worker_exception_mutex);
    if (worker_exception == nullptr) {
      worker_exception = std::current_exception();
    }
  }
}

template <typename Worker>
auto worker_thread_indexed(Worker &worker, const IndexRange &range,
                           std::size_t chunk_id,
                           std::exception_ptr &worker_exception,
                           std::mutex &worker_exception_mutex) -> void {
  try {
    worker(range.first, range.second, chunk_id);
  } catch (...) {
    std::lock_guard<std::mutex> guard(worker_exception_mutex);
    if (worker_exception == nullptr) {
      worker_exception = std::current_exception();
    }
  }
}

// Function to calculate the ranges for a given input
inline auto effective_grain_size(const IndexRange &range, std::size_t n_threads,
                                 std::size_t grain_size) -> std::size_t {
  if (range.first >= range.second) {
    return 0;
  }
  const auto length = range.second - range.first;

  grain_size = (std::max)(grain_size, static_cast<std::size_t>(1));
  if (n_threads <= 1) {
    return length;
  }

  const auto max_chunks = (std::min)(length, n_threads);
  const auto even_chunk = static_cast<std::size_t>(1) +
    ((length - 1) / max_chunks);
  return (std::max)(grain_size, even_chunk);
}

inline auto split_input_range(const IndexRange &range, std::size_t n_threads,
                              std::size_t grain_size)
    -> std::vector<IndexRange> {
  if (range.first >= range.second) {
    return {};
  }
  grain_size = effective_grain_size(range, n_threads, grain_size);

  // allocate ranges
  std::vector<IndexRange> ranges;
  ranges.reserve(1 + ((range.second - range.first - 1) / grain_size));
  std::size_t begin = range.first;
  while (begin < range.second) {
    const auto remaining = range.second - begin;
    const auto end =
      remaining <= grain_size ? range.second : begin + grain_size;
    ranges.emplace_back(begin, end);
    begin = end;
  }

  return ranges;
}

// Execute the Worker over the IndexRange in parallel.
template <typename Worker>
inline void parallel_for(std::size_t begin, std::size_t end, Worker &worker,
                         std::size_t n_threads, std::size_t grain_size = 1) {
  if (begin >= end) {
    return;
  }
  if (n_threads <= 1) {
    worker(begin, end);
    return;
  }
  // split the work
  IndexRange input_range(begin, end);
  std::vector<IndexRange> ranges =
    split_input_range(input_range, n_threads, grain_size);
  if (ranges.size() <= 1) {
    worker(begin, end);
    return;
  }

  std::exception_ptr worker_exception = nullptr;
  std::mutex worker_exception_mutex;
  std::vector<std::thread> threads;
  threads.reserve(ranges.size());
  ThreadJoiner thread_joiner(threads);
  for (auto &range : ranges) {
    threads.emplace_back(&worker_thread<Worker>, std::ref(worker), range,
                         std::ref(worker_exception),
                         std::ref(worker_exception_mutex));
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  if (worker_exception != nullptr) {
    std::rethrow_exception(worker_exception);
  }

  return;
}

template <typename Worker>
inline void parallel_for(std::size_t end, Worker &worker, std::size_t n_threads,
                         std::size_t grain_size = 1) {
  parallel_for(0, end, worker, n_threads, grain_size);
}

// Execute the Worker over the IndexRange in parallel, passing a zero-based
// chunk index as the third worker argument.
template <typename Worker>
inline void parallel_for_indexed(std::size_t begin, std::size_t end,
                                 Worker &worker, std::size_t n_threads,
                                 std::size_t grain_size = 1) {
  if (begin >= end) {
    return;
  }
  if (n_threads <= 1) {
    worker(begin, end, 0);
    return;
  }
  IndexRange input_range(begin, end);
  std::vector<IndexRange> ranges =
    split_input_range(input_range, n_threads, grain_size);
  if (ranges.size() <= 1) {
    worker(begin, end, 0);
    return;
  }

  std::exception_ptr worker_exception = nullptr;
  std::mutex worker_exception_mutex;
  std::vector<std::thread> threads;
  threads.reserve(ranges.size());
  ThreadJoiner thread_joiner(threads);
  for (std::size_t chunk_id = 0; chunk_id < ranges.size(); ++chunk_id) {
    threads.emplace_back(&worker_thread_indexed<Worker>, std::ref(worker),
                         ranges[chunk_id], chunk_id,
                         std::ref(worker_exception),
                         std::ref(worker_exception_mutex));
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  if (worker_exception != nullptr) {
    std::rethrow_exception(worker_exception);
  }

  return;
}

template <typename Worker>
inline void parallel_for_indexed(std::size_t end, Worker &worker,
                                 std::size_t n_threads,
                                 std::size_t grain_size = 1) {
  parallel_for_indexed(0, end, worker, n_threads, grain_size);
}

} // namespace pforr

#endif // PFORR
