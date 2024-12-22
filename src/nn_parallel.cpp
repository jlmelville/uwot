#include <vector>

#include <Rcpp.h>

#include "RcppPerpendicular.h"
#include "nn_parallel.h"

using namespace Rcpp;

template <typename UwotAnnoyDistance>
auto annoy_nns_impl(const std::string &index_name, NumericMatrix mat,
                    std::size_t n_neighbors, std::size_t search_k,
                    std::size_t n_threads = 0,
                    std::size_t grain_size = 1) -> List {

  std::size_t nrow = mat.rows();
  std::size_t ncol = mat.cols();

  std::vector<double> vmat = as<std::vector<double>>(mat);

  NNWorker<UwotAnnoyDistance> worker(index_name, vmat, ncol, n_neighbors,
                                     search_k);
  RcppPerpendicular::parallel_for(0, nrow, worker, n_threads, grain_size);

  return List::create(
      _("item") = IntegerMatrix(nrow, n_neighbors, worker.idx.begin()),
      _("distance") = NumericMatrix(nrow, n_neighbors, worker.dists.begin()));
}

// [[Rcpp::export]]
List annoy_search_parallel_cpp(const std::string &index_name, NumericMatrix mat,
                               std::size_t n_neighbors, std::size_t search_k,
                               const std::string &metric,
                               std::size_t n_threads = 0,
                               std::size_t grain_size = 1) {
  if (metric == "euclidean") {
    return annoy_nns_impl<UwotAnnoyEuclidean>(index_name, mat, n_neighbors,
                                              search_k, n_threads, grain_size);
  } else if (metric == "cosine") {
    return annoy_nns_impl<UwotAnnoyCosine>(index_name, mat, n_neighbors,
                                           search_k, n_threads, grain_size);
  } else if (metric == "manhattan") {
    return annoy_nns_impl<UwotAnnoyManhattan>(index_name, mat, n_neighbors,
                                              search_k, n_threads, grain_size);
  } else if (metric == "hamming") {
    return annoy_nns_impl<UwotAnnoyHamming>(index_name, mat, n_neighbors,
                                            search_k, n_threads, grain_size);
  } else {
    stop("Unknown metric '", metric, "'");
  }
}
