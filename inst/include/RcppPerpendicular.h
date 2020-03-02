// Taken from RcppParallel.h and then modified slightly to rename header guards
// and namespaces to avoid any potential clashes. RcppParallel is licensed under
// GPLv2


#ifndef __RCPP_PERPENDICULAR__
#define __RCPP_PERPENDICULAR__

// TinyThread implementation
#include "RcppPerpendicular/StdThread.h"

#include "RcppPerpendicular/RMatrix.h"

namespace RcppPerpendicular {

inline void parallelFor(std::size_t begin, std::size_t end, 
                        Worker& worker, std::size_t grainSize = 1) {
   stParallelFor(begin, end, worker, grainSize);
}

} // namespace RcppPerpendicular

#endif // __RCPP_PERPENDICULAR__
