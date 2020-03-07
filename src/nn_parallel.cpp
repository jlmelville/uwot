#include <vector>

#include "RcppPerpendicular.h"
#include <Rcpp.h>

#include "nn_parallel.h"

template <typename UwotAnnoyDistance>
auto annoy_nns_impl(const std::string &index_name,
                    const Rcpp::NumericMatrix &mat, std::size_t n_neighbors,
                    std::size_t search_k, std::size_t n_threads = 0,
                    std::size_t grain_size = 1, bool verbose = false)
    -> Rcpp::List {

  std::size_t nrow = mat.rows();
  std::size_t ncol = mat.cols();

  std::vector<double> vmat = Rcpp::as<std::vector<double>>(mat);

  NNWorker<UwotAnnoyDistance> worker(index_name, vmat, ncol, n_neighbors,
                                     search_k);
  RcppPerpendicular::parallel_for(0, nrow, worker, n_threads, grain_size);

  return Rcpp::List::create(Rcpp::_("item") = Rcpp::IntegerMatrix(
                                nrow, n_neighbors, worker.idx.begin()),
                            Rcpp::_("distance") = Rcpp::NumericMatrix(
                                nrow, n_neighbors, worker.dists.begin()));
}

// [[Rcpp::export]]
Rcpp::List
annoy_search_parallel_cpp(const std::string &index_name,
                          const Rcpp::NumericMatrix &mat,
                          std::size_t n_neighbors, std::size_t search_k,
                          const std::string &metric, std::size_t n_threads = 0,
                          std::size_t grain_size = 1, bool verbose = false) {
  if (metric == "euclidean") {
    return annoy_nns_impl<UwotAnnoyEuclidean>(index_name, mat, n_neighbors,
                                              search_k, grain_size, verbose);
  } else if (metric == "cosine") {
    return annoy_nns_impl<UwotAnnoyCosine>(index_name, mat, n_neighbors,
                                           search_k, grain_size, verbose);
  } else if (metric == "manhattan") {
    return annoy_nns_impl<UwotAnnoyManhattan>(index_name, mat, n_neighbors,
                                              search_k, grain_size, verbose);
  } else if (metric == "hamming") {
    return annoy_nns_impl<UwotAnnoyHamming>(index_name, mat, n_neighbors,
                                            search_k, grain_size, verbose);
  } else {
    Rcpp::stop("Unknown metric '", metric, "'");
  }
}
