#include <thread>

#include <Rcpp.h>

// [[Rcpp::export]]
unsigned int hardware_concurrency() {
  return std::thread::hardware_concurrency();
}
