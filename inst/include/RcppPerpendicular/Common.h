// Taken from RcppParallel's Common.h and then modified slightly to rename
// header guards and namespaces to avoid any potential clashes. RcppParallel is
// licensed under GPLv2

#ifndef RCPP_PERPENDICULAR_COMMON
#define RCPP_PERPENDICULAR_COMMON

#include <cstddef>

namespace RcppPerpendicular {

// Work executed within a background thread. We implement dynamic
// dispatch using vtables so we can have a stable type to cast
// to from the void* passed to the worker thread (required because
// the tinythreads interface allows to pass only a void* to the
// thread main rather than a generic type / template)

struct Worker
{
   // construct and destruct (delete virtually)
   Worker() {}
   virtual ~Worker() {}

   // dispatch work over a range of values
   virtual void operator()(std::size_t begin, std::size_t end) = 0;

   // disable copying and assignment
private:
   Worker(const Worker&);
   void operator=(const Worker&);
};

} // namespace RcppPerpendicular


#endif // RCPP_PERPENDICULAR_COMMON
