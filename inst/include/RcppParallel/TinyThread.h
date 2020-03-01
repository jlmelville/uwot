#ifndef __RCPP_PARALLEL_TINYTHREAD__
#define __RCPP_PARALLEL_TINYTHREAD__

#include <cstdlib>

#include "Common.h"

#include <tthread/tinythread.h>

#include <stdio.h>

#include <vector>

namespace RcppParallel {

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
   std::size_t threads = tthread::thread::hardware_concurrency();
   char* numThreads = ::getenv("RCPP_PARALLEL_NUM_THREADS");
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
      grainSize = std::max(length / threads, grainSize);
   else // imperfect division, divide by threads - 1
      grainSize = std::max(length / (threads-1), grainSize);
  
   // allocate ranges
   std::vector<IndexRange> ranges;
   std::size_t begin = range.begin();
   while (begin < range.end()) {
      std::size_t end = std::min(begin + grainSize, range.end());
      ranges.push_back(IndexRange(begin, end));      
      begin = end;
   }
   
   // return ranges  
   return ranges;
}

} // anonymous namespace

// Execute the Worker over the IndexRange in parallel
inline void ttParallelFor(std::size_t begin, std::size_t end, 
                          Worker& worker, std::size_t grainSize = 1) {
  
   // split the work
   IndexRange inputRange(begin, end);
   std::vector<IndexRange> ranges = splitInputRange(inputRange, grainSize);
   
   // create threads
   std::vector<tthread::thread*> threads;
   for (std::size_t i = 0; i<ranges.size(); ++i) {
      threads.push_back(new tthread::thread(workerThread, new Work(ranges[i], worker)));
   }
   
   // join and delete them
   for (std::size_t i = 0; i<threads.size(); ++i) {
      threads[i]->join();
      delete threads[i];
   }
}

// Execute the IWorker over the range in parallel then join results
template <typename Reducer>
inline void ttParallelReduce(std::size_t begin, std::size_t end, 
                             Reducer& reducer, std::size_t grainSize = 1) {
  
   // split the work
   IndexRange inputRange(begin, end);
   std::vector<IndexRange> ranges = splitInputRange(inputRange, grainSize);
   
   // create threads (split for each thread and track the allocated workers)
   std::vector<tthread::thread*> threads;
   std::vector<Worker*> workers;
   for (std::size_t i = 0; i<ranges.size(); ++i) {
      Reducer* pReducer = new Reducer(reducer, RcppParallel::Split());
      workers.push_back(pReducer);
      threads.push_back(new tthread::thread(workerThread, new Work(ranges[i], *pReducer)));
   }
   
   // wait for each thread, join it's results, then delete the worker & thread
   for (std::size_t i = 0; i<threads.size(); ++i) {
      // wait for thread
      threads[i]->join();
      
      // join the results
      reducer.join(static_cast<Reducer&>(*workers[i]));
      
      // delete the worker (which we split above) and the thread
      delete workers[i];
      delete threads[i];
   }
}

} // namespace RcppParallel

#endif // __RCPP_PARALLEL_TINYTHREAD__
