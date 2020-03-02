// Taken from RcppParallel's TinyThread.h and then modified slightly to work
// with std::thread and to rename header guards and namespaces to avoid any
// potential clashes. RcppParallel is licensed under GPLv2
#ifndef RCPP_PERPENDICULAR_STDTHREAD
#define RCPP_PERPENDICULAR_STDTHREAD

#include <cstdlib>
#include <stdio.h>
#include <thread>
#include <vector>

#include "Common.h"

namespace RcppPerpendicular {

namespace {

// Class which represents a range of indexes to perform work on
// (worker functions are passed this range so they know which
// elements are safe to read/write to)
class IndexRange {
public:

   // Initizlize with a begin and (exclusive) end index
   IndexRange(std::size_t begin, std::size_t end)
    : begin_(begin), end_(end)
   {
   }

   // Access begin() and end()
   std::size_t begin() const { return begin_; }
   std::size_t end() const { return end_; }
   std::size_t size() const { return end_ - begin_ ; }

private:
   std::size_t begin_;
   std::size_t end_;
};


// Because tinythread allows us to pass only a plain C function
// we need to pass our worker and range within a struct that we
// can cast to/from void*
struct Work {
   Work(IndexRange range, Worker& worker)
    :  range(range), worker(worker)
   {
   }
   IndexRange range;
   Worker& worker;
};

// Thread which performs work (then deletes the work object
// when it's done)
extern "C" inline void workerThread(void* data) {
   try
   {
      Work* pWork = static_cast<Work*>(data);
      pWork->worker(pWork->range.begin(), pWork->range.end());
      delete pWork;
   }
   catch(...)
   {
   }
}

// Function to calculate the ranges for a given input
std::vector<IndexRange> splitInputRange(const IndexRange& range,
                                        std::size_t grainSize) {

   // determine max number of threads
   std::size_t threads = std::thread::hardware_concurrency();
   char* numThreads = ::getenv("RCPP_PERPENDICULAR_NUM_THREADS");
   if (numThreads != NULL) {
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
      grainSize = (std::max)(length / (threads-1), grainSize);

   // allocate ranges
   std::vector<IndexRange> ranges;
   std::size_t begin = range.begin();
   while (begin < range.end()) {
      std::size_t end = (std::min)(begin + grainSize, range.end());
      ranges.push_back(IndexRange(begin, end));
      begin = end;
   }

   return ranges;
}

} // anonymous namespace

// Execute the Worker over the IndexRange in parallel
inline void stParallelFor(std::size_t begin, std::size_t end,
                          Worker& worker, std::size_t grainSize = 1) {

   // split the work
   IndexRange inputRange(begin, end);
   std::vector<IndexRange> ranges = splitInputRange(inputRange, grainSize);

   // create threads
   std::vector<std::thread*> threads;
   for (std::size_t i = 0; i<ranges.size(); ++i) {
      threads.push_back(new std::thread(workerThread, new Work(ranges[i], worker)));
   }

   // join and delete them
   for (std::size_t i = 0; i<threads.size(); ++i) {
      threads[i]->join();
      delete threads[i];
   }
}

} // namespace RcppPerpendicular

#endif // RCPP_PERPENDICULAR_STDTHREAD
