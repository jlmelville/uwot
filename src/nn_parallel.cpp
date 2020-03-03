#include <utility>
#include <vector>

#include <Rcpp.h>
#include "RcppPerpendicular.h"


#if defined(__MINGW32__)
#undef Realloc
#undef Free
#endif

#define __ERROR_PRINTER_OVERRIDE__ REprintf

#include <annoylib.h>
#include <kissrandom.h>

#include "matrix.h"

template <typename S, typename T, typename Distance, typename Random>
struct NNWorker : public RcppPerpendicular::Worker {
  std::string index_name;
  const std::vector<double> &mat;
  std::vector<T> &dists;
  std::vector<int> &idx;
  const std::size_t nrow;
  std::size_t ncol;
  std::size_t n_neighbors;
  std::size_t search_k;

  NNWorker(const std::string &index_name, const std::vector<double> &mat,
           std::vector<T> &dists, std::vector<int> &idx,
           std::size_t ncol, std::size_t n_neighbors, std::size_t search_k)
      : index_name(index_name), mat(mat), dists(dists), idx(idx), 
        nrow(mat.size() / ncol), ncol(ncol), n_neighbors(n_neighbors), 
        search_k(search_k)
 {}

  void operator()(std::size_t begin, std::size_t end) {
    AnnoyIndex<S, T, Distance, Random> index(ncol);
    index.load(index_name.c_str());

    for (std::size_t i = begin; i < end; i++) {
      std::vector<T> fv(ncol);
      row(mat, nrow, ncol, i, fv);

      std::vector<S> result;
      std::vector<T> distances;

      index.get_nns_by_vector(fv.data(), n_neighbors, search_k, &result,
                              &distances);
      if (result.size() != n_neighbors || distances.size() != n_neighbors) {
        break;
      }

      for (std::size_t j = 0; j < n_neighbors; j++) {
        set_row(dists, nrow, n_neighbors, i, distances);
        set_row(idx, nrow, n_neighbors, i, result);
      }
    }
  }
};

// [[Rcpp::export]]
Rcpp::List annoy_euclidean_nns(const std::string &index_name,
                               const Rcpp::NumericMatrix &mat,
                               std::size_t n_neighbors, std::size_t search_k,
                               std::size_t grain_size = 1,
                               bool verbose = false) {
  std::size_t nrow = mat.rows();
  std::size_t ncol = mat.cols();

  std::vector<float> dist(nrow * n_neighbors);
  std::vector<int> idx(nrow * n_neighbors, -1);
  
  std::vector<double> vmat = Rcpp::as<std::vector<double>>(mat);

  NNWorker<int32_t, float, Euclidean, Kiss64Random> worker(
      index_name, vmat, dist, idx, ncol, n_neighbors, search_k);
  RcppPerpendicular::parallelFor(0, nrow, worker, grain_size);

  return Rcpp::List::create(Rcpp::Named("item") = std::move(Rcpp::IntegerMatrix(nrow, n_neighbors, idx.begin())),
                            Rcpp::Named("distance") = std::move(Rcpp::NumericMatrix(nrow, n_neighbors, dist.begin())));
}

// [[Rcpp::export]]
Rcpp::List annoy_cosine_nns(const std::string &index_name,
                            const Rcpp::NumericMatrix &mat,
                            std::size_t n_neighbors, std::size_t search_k,
                            std::size_t grain_size = 1, bool verbose = false) {
  std::size_t nrow = mat.rows();
  std::size_t ncol = mat.cols();

  std::vector<float> dist(nrow * n_neighbors);
  std::vector<int> idx(nrow * n_neighbors, -1);
  
  std::vector<double> vmat = Rcpp::as<std::vector<double>>(mat);
  
  NNWorker<int32_t, float, Angular, Kiss64Random> worker(
      index_name, vmat, dist, idx, ncol, n_neighbors, search_k);
  RcppPerpendicular::parallelFor(0, nrow, worker, grain_size);

  return Rcpp::List::create(Rcpp::Named("item") = std::move(Rcpp::IntegerMatrix(nrow, n_neighbors, idx.begin())),
                            Rcpp::Named("distance") = std::move(Rcpp::NumericMatrix(nrow, n_neighbors, dist.begin())));
}

// [[Rcpp::export]]
Rcpp::List annoy_manhattan_nns(const std::string &index_name,
                               const Rcpp::NumericMatrix &mat,
                               std::size_t n_neighbors, std::size_t search_k,
                               std::size_t grain_size = 1,
                               bool verbose = false) {
  std::size_t nrow = mat.rows();
  std::size_t ncol = mat.cols();

  std::vector<float> dist(nrow * n_neighbors);
  std::vector<int> idx(nrow * n_neighbors, -1);
  
  std::vector<double> vmat = Rcpp::as<std::vector<double>>(mat);

  NNWorker<int32_t, float, Manhattan, Kiss64Random> worker(
      index_name, vmat, dist, idx, ncol, n_neighbors, search_k);

  RcppPerpendicular::parallelFor(0, nrow, worker, grain_size);

  return Rcpp::List::create(Rcpp::Named("item") = std::move(Rcpp::IntegerMatrix(nrow, n_neighbors, idx.begin())),
                            Rcpp::Named("distance") = std::move(Rcpp::NumericMatrix(nrow, n_neighbors, dist.begin())));
}

// [[Rcpp::export]]
Rcpp::List annoy_hamming_nns(const std::string &index_name,
                             const Rcpp::NumericMatrix &mat,
                             std::size_t n_neighbors, std::size_t search_k,
                             std::size_t grain_size = 1, bool verbose = false) {
  std::size_t nrow = mat.rows();
  std::size_t ncol = mat.cols();

  std::vector<uint64_t> dist(nrow * n_neighbors);
  std::vector<int> idx(nrow * n_neighbors, -1);
  
  std::vector<double> vmat = Rcpp::as<std::vector<double>>(mat);
  
  NNWorker<int32_t, uint64_t, Hamming, Kiss64Random> worker(
      index_name, vmat, dist, idx, ncol, n_neighbors, search_k);

  RcppPerpendicular::parallelFor(0, nrow, worker, grain_size);

  return Rcpp::List::create(Rcpp::Named("item") = std::move(Rcpp::IntegerMatrix(nrow, n_neighbors, idx.begin())),
                            Rcpp::Named("distance") = std::move(Rcpp::NumericMatrix(nrow, n_neighbors, dist.begin())));
}
