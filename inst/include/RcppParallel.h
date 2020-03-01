
#ifndef __RCPP_PARALLEL__
#define __RCPP_PARALLEL__

// TinyThread implementation
#include "RcppParallel/TinyThread.h"

#include "RcppParallel/RVector.h"
#include "RcppParallel/RMatrix.h"

namespace RcppParallel {

inline void parallelFor(std::size_t begin, std::size_t end, 
                        Worker& worker, std::size_t grainSize = 1) {
   ttParallelFor(begin, end, worker, grainSize);
}

template <typename Reducer>
inline void parallelReduce(std::size_t begin, std::size_t end, 
                           Reducer& reducer, std::size_t grainSize = 1) {
   ttParallelReduce(begin, end, reducer, grainSize);
}

} // namespace RcppParallel

#endif // __RCPP_PARALLEL__
