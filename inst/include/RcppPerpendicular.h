// Taken from RcppParallel.h and then modified slightly to rename header guards
// and namespaces to avoid any potential clashes. RcppParallel is licensed under
// GPLv2

#ifndef RCPP_PERPENDICULAR
#define RCPP_PERPENDICULAR

#include <thread>
#include <vector>

namespace RcppPerpendicular {

// Class which represents a range of indexes to perform work on
// (worker functions are passed this range so they know which
// elements are safe to read/write to)
class IndexRange {
public:
  // Initialize with a begin and (exclusive) end index
  IndexRange(std::size_t begin, std::size_t end) : begin_(begin), end_(end) {}

  // Access begin() and end()
  auto begin() const -> std::size_t { return begin_; }
  auto end() const -> std::size_t { return end_; }
  auto size() const -> std::size_t { return end_ - begin_; }

private:
  std::size_t begin_;
  std::size_t end_;
};

template <typename Worker>
auto worker_thread(Worker &worker, IndexRange range) -> void {
  try {
    worker(range.begin(), range.end());
  } catch (...) {
  }
}

// Function to calculate the ranges for a given input
inline auto split_input_range(const IndexRange &range, std::size_t grain_size)
    -> std::vector<IndexRange> {

  // determine max number of threads
  std::size_t threads = std::thread::hardware_concurrency();
  char *numThreads = ::getenv("RCPP_PERPENDICULAR_NUM_THREADS");
  if (numThreads != nullptr) {
    int parsedThreads = ::atoi(numThreads);
    if (parsedThreads > 0)
      threads = parsedThreads;
  }

  // compute grain_size (including enforcing requested minimum)
  std::size_t length = range.end() - range.begin();
  if (threads == 1)
    grain_size = length;
  else if ((length % threads) == 0) // perfect division
    grain_size = (std::max)(length / threads, grain_size);
  else // imperfect division, divide by threads - 1
    grain_size = (std::max)(length / (threads - 1), grain_size);

  // allocate ranges
  std::vector<IndexRange> ranges;
  std::size_t begin = range.begin();
  while (begin < range.end()) {
    std::size_t end = (std::min)(begin + grain_size, range.end());
    ranges.emplace_back(IndexRange(begin, end));
    begin = end;
  }

  return ranges;
}

// Execute the Worker over the IndexRange in parallel
template <typename Worker>
inline void parallel_for(std::size_t begin, std::size_t end, Worker &worker,
                        std::size_t grain_size = 1) {
  // split the work
  IndexRange input_range(begin, end);
  std::vector<IndexRange> ranges = split_input_range(input_range, grain_size);

  std::vector<std::thread> threads;
  for (auto &range : ranges) {
    threads.push_back(std::thread(&worker_thread<Worker>, std::ref(worker), range));
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

} // namespace RcppPerpendicular

#endif // RCPP_PERPENDICULAR
