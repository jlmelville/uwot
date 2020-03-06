// Taken from RcppParallel.h and then modified slightly to rename header guards
// and namespaces to avoid any potential clashes. RcppParallel is licensed under
// GPLv2

#ifndef RCPP_PERPENDICULAR
#define RCPP_PERPENDICULAR

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

namespace RcppPerpendicular {

// Work executed within a background thread. We implement dynamic
// dispatch using vtables so we can have a stable type to cast
// to from the void* passed to the worker thread (required because
// the tinythreads interface allows to pass only a void* to the
// thread main rather than a generic type / template)

struct Worker {
  // construct and destruct (delete virtually)
  Worker() = default;
  virtual ~Worker() = default;

  // dispatch work over a range of values
  virtual void operator()(std::size_t begin, std::size_t end) = 0;

  // disable copying and assignment
  Worker(const Worker &) = delete;
  void operator=(const Worker &) = delete;
};

namespace {

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

// Because tinythread allows us to pass only a plain C function
// we need to pass our worker and range within a struct that we
// can cast to/from void*
struct Work {
  Work(IndexRange range, Worker &worker) : range(range), worker(worker) {}
  IndexRange range;
  Worker &worker;
};

// Thread which performs work (then deletes the work object
// when it's done)
extern "C" inline void workerThread(void *data) {
  try {
    Work *pWork = static_cast<Work *>(data);
    pWork->worker(pWork->range.begin(), pWork->range.end());
    delete pWork;
  } catch (...) {
  }
}

// Function to calculate the ranges for a given input
auto splitInputRange(const IndexRange &range, std::size_t grainSize)
    -> std::vector<IndexRange> {

  // determine max number of threads
  std::size_t threads = std::thread::hardware_concurrency();
  char *numThreads = ::getenv("RCPP_PERPENDICULAR_NUM_THREADS");
  if (numThreads != nullptr) {
    int parsedThreads = ::atoi(numThreads);
    if (parsedThreads > 0)
      threads = parsedThreads;
  }

  // compute grainSize (including enforcing requested minimum)
  std::size_t length = range.end() - range.begin();
  if (threads == 1)
    grainSize = length;
  else if ((length % threads) == 0) // perfect division
    grainSize = (std::max)(length / threads, grainSize);
  else // imperfect division, divide by threads - 1
    grainSize = (std::max)(length / (threads - 1), grainSize);

  // allocate ranges
  std::vector<IndexRange> ranges;
  std::size_t begin = range.begin();
  while (begin < range.end()) {
    std::size_t end = (std::min)(begin + grainSize, range.end());
    ranges.emplace_back(IndexRange(begin, end));
    begin = end;
  }

  return ranges;
}

} // anonymous namespace

// Execute the Worker over the IndexRange in parallel
inline void parallelFor(std::size_t begin, std::size_t end, Worker &worker,
                        std::size_t grainSize = 1) {
  // split the work
  IndexRange inputRange(begin, end);
  std::vector<IndexRange> ranges = splitInputRange(inputRange, grainSize);

  // create threads
  std::vector<std::thread *> threads;
  for (auto &range : ranges) {
    threads.push_back(new std::thread(workerThread, new Work(range, worker)));
  }

  // join and delete them
  for (auto &thread : threads) {
    thread->join();
    delete thread;
  }
}

} // namespace RcppPerpendicular

#endif // RCPP_PERPENDICULAR
