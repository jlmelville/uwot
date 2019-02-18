//  UWOT -- An R package for dimensionality reduction using UMAP
//
//  Copyright (C) 2019 James Melville
//
//  This file is part of UWOT
//
//  UWOT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  UWOT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with UWOT.  If not, see <http://www.gnu.org/licenses/>.

#include <bitset>
#include <Rcpp.h>

struct Heap
{
  std::vector<std::vector<int>> idx;
  std::vector<std::vector<double>> dist;
  std::vector<std::vector<bool>> flags; // vector of bool, yes ugh
  Heap(const std::size_t n_points, const std::size_t size) {
    for (std::size_t i = 0; i < n_points; i++) {
      idx.push_back(std::vector<int>(size, -1));
      dist.push_back(std::vector<double>(size, std::numeric_limits<double>::max()));
      flags.push_back(std::vector<bool>(size, false));
    }
  }
  
  unsigned int push(std::size_t row, double weight, std::size_t index, bool flag) {
    std::vector<int>& indices = idx[row];
    std::vector<double>& weights = dist[row];
    std::vector<bool>& is_new = flags[row];
    
    if (weight >= weights[0]) {
      return 0;
    }
    
    // break if we already have this element
    int iindex = static_cast<int>(index);
    const std::size_t n_nbrs = indices.size();
    for (std::size_t i = 0; i < n_nbrs; i++) {
      if (iindex == indices[i]) {
        return 0;
      }
    }
    
    // insert val at position zero
    weights[0] = weight;
    indices[0] = iindex;
    is_new[0] = flag;
    
    // descend the heap, swapping values until the max heap criterion is met
    std::size_t i = 0;
    std::size_t i_swap = 0;
    while (true) {
      std::size_t ic1 = 2 * i + 1;
      std::size_t ic2 = ic1 + 1;
      
      if (ic1 >= n_nbrs) {
        break;
      }
      else if (ic2 >= n_nbrs) {
        if (weights[ic1] >= weight) {
          i_swap = ic1;
        }
        else {
          break;
        }
      }
      else if (weights[ic1] >= weights[ic2]) {
        if (weight < weights[ic1]) {
          i_swap = ic1;
        }
        else {
          break;
        }
      }
      else {
        if (weight < weights[ic2]) {
          i_swap = ic2;
        }
        else {
          break;
        }
      }
      
      weights[i] = weights[i_swap];
      indices[i] = indices[i_swap];
      is_new[i] = is_new[i_swap];
      
      i = i_swap;
    }
    
    weights[i] = weight;
    indices[i] = iindex;
    is_new[i] = flag;
    
    return 1;
  }
  
  void deheap_sort() {
    const std::size_t npoints = idx.size();
    
    for (std::size_t i = 0; i < npoints; i++) {
      std::vector<int>& ind_heap = idx[i];
      std::vector<double>& dist_heap = dist[i];
      
      const std::size_t nnbrs = ind_heap.size();
      for (std::size_t j = 0; j < nnbrs - 1; j++) {
        std::swap(ind_heap[0], ind_heap[nnbrs - j - 1]);
        std::swap(dist_heap[0], dist_heap[nnbrs - j - 1]);
        siftdown(dist_heap, ind_heap, nnbrs - j - 1, 0);
      }
    }
  }
  
  void siftdown(std::vector<double>& dist_heap,
                std::vector<int>& ind_heap,
                const std::size_t len,
                std::size_t elt) {
    
    while (elt * 2 + 1 < len) {
      std::size_t left_child = elt * 2 + 1;
      std::size_t right_child = left_child + 1;
      std::size_t swap = elt;
      
      if (dist_heap[swap] < dist_heap[left_child]) {
        swap = left_child;
      }
      
      if (right_child < len && dist_heap[swap] < dist_heap[right_child]) {
        swap = right_child;
      }
      
      if (swap == elt) {
        break;
      }
      else {
        std::swap(dist_heap[elt], dist_heap[swap]);
        std::swap(ind_heap[elt], ind_heap[swap]);
        elt = swap;
      }
    }
  }
};

template <typename In, typename Out>
struct Euclidean
{
  Euclidean(const std::vector<In>& data, std::size_t ndim)
    : data(data), ndim(ndim) { }
  
  Out operator()(std::size_t i, std::size_t j) {
    Out sum = 0.0;
    const std::size_t di = ndim * i;
    const std::size_t dj = ndim * j;
    
    for (std::size_t d = 0; d < ndim; d++) {
      const Out diff = data[di + d] - data[dj + d];
      sum += diff * diff;
    }
    
    return std::sqrt(sum);
  }
  
  std::vector<In> data;
  std::size_t ndim;
  
  typedef In in_type;
};

template <typename In, typename Out>
struct L2
{
  L2(const std::vector<In>& data, std::size_t ndim)
    : data(data), ndim(ndim) { }
  
  Out operator()(std::size_t i, std::size_t j) {
    Out sum = 0.0;
    const std::size_t di = ndim * i;
    const std::size_t dj = ndim * j;
    
    for (std::size_t d = 0; d < ndim; d++) {
      const Out diff = data[di + d] - data[dj + d];
      sum += diff * diff;
    }
    
    return sum;
  }
  
  std::vector<In> data;
  std::size_t ndim;
  
  typedef In in_type;
};

template <typename In, typename Out>
struct Cosine
{
  Cosine(const std::vector<In>& _data, std::size_t ndim)
    : data(_data), ndim(ndim)
  
  {
    // normalize data on input
    const std::size_t npoints = data.size() / ndim;
    for (std::size_t i = 0; i < npoints; i++) {
      const std::size_t di = ndim * i;
      In norm = 0.0;
      
      for (std::size_t d = 0; d < ndim; d++) {
        const auto val = data[di + d];
        norm += val * val;
      }
      norm = 1.0 / (std::sqrt(norm) + 1e-30);
      for (std::size_t d = 0; d < ndim; d++) {
        data[di + d] *= norm;
      }
    }
  }
  
  Out operator()(std::size_t i, std::size_t j) {
    const std::size_t di = ndim * i;
    const std::size_t dj = ndim * j;
    
    Out sum = 0.0;
    for (std::size_t d = 0; d < ndim; d++) {
      sum += data[di + d] * data[dj + d];
    }
    
    return 1.0 - sum;
  }
  
  std::vector<In> data;
  const std::size_t ndim;
  
  typedef In in_type;
};


template <typename In, typename Out>
struct Manhattan
{
  
  Manhattan(const std::vector<In>& data, std::size_t ndim)
    : data(data), ndim(ndim) { }
  
  Out operator()(std::size_t i, std::size_t j) {
    Out sum = 0.0;
    const std::size_t di = ndim * i;
    const std::size_t dj = ndim * j;
    
    for (std::size_t d = 0; d < ndim; d++) {
      sum += std::abs(data[di + d] - data[dj + d]);
    }
    
    return sum;
  }
  
  std::vector<In> data;
  std::size_t ndim;
  
  typedef In in_type;
};

template<typename In, typename Out>
struct Hamming
{
  Hamming(const std::vector<In>& vdata, std::size_t vndim) {
    // Instead of storing each bit as an element, we will pack them
    // into a series of 64-bit bitsets. Possibly compilers are smart enough
    // to use built in integer popcount routines for the bitset count()
    // method.
    std::bitset<64> bits;
    std::size_t bit_count = 0;
    std::size_t vd_count = 0;
    
    for (std::size_t i = 0; i < vdata.size(); i++) {
      if (bit_count == 64 || vd_count == vndim) {
        // filled up current bitset
        data.push_back(bits);
        bit_count = 0;
        bits.reset();
        
        if (vd_count == vndim) {
          // end of item
          vd_count = 0;
        }
      }
      bits[bit_count] = vdata[i];
      
      ++vd_count;
      ++bit_count;
    }
    if (bit_count > 0) {
      data.push_back(bits);
    }
    
    ndim = std::ceil(vndim / 64.0);
  }
  
  Out operator()(std::size_t i, std::size_t j) {
    Out sum = 0;
    const std::size_t di = ndim * i;
    const std::size_t dj = ndim * j;
    
    for (std::size_t d = 0; d < ndim; d++) {
      sum += (data[di + d] ^ data[dj + d]).count();
    }
    
    return sum;
  }
  
  std::vector<std::bitset<64>> data;
  std::size_t ndim;
  
  typedef In in_type;
};


struct RRand {
  // a random uniform value between 0 and 1
  double unif() {
    return Rcpp::runif(1, 0.0, 1.0)[0];
  }
};

struct RProgress {
  void iter(std::size_t n, std::size_t n_iters, const Heap& heap) {
    double sum = 0.0;
    for (std::size_t i = 0; i < heap.dist.size(); i++) {
      for (std::size_t j = 0; j < heap.dist[i].size(); j++) {
        sum += heap.dist[i][j];
      }
    }
    Rcpp::Rcout << (n + 1) << " / " << n_iters << " " << sum << std::endl;
  }
  void converged(const std::size_t c, const double tol) {
    Rcpp::Rcout << "c = " << c << " tol = " << tol << std::endl;
  }
  void check_interrupt() {
    Rcpp::checkUserInterrupt();
  }
};

template <typename Rand>
void build_candidates(Heap& current_graph, Heap& candidate_neighbors,
                      const std::size_t npoints, const std::size_t nnbrs) {
  Rand rand;
  for (std::size_t i = 0; i < npoints; i++) {
    for (std::size_t j = 0; j < nnbrs; j++) {
      if (current_graph.idx[i][j] < 0) {
        continue;
      }
      int idx = current_graph.idx[i][j];
      bool isn = current_graph.flags[i][j];
      double d = rand.unif();
      
      candidate_neighbors.push(i, d, idx, isn);
      candidate_neighbors.push(static_cast<std::size_t>(idx), d, i, isn);
      
      current_graph.flags[i][j] = false;
    }
  }
}

Heap r_to_heap(
    Rcpp::IntegerMatrix idx,
    Rcpp::NumericMatrix dist
) {
  const std::size_t npoints = idx.nrow();
  const std::size_t nnbrs = idx.ncol();
  
  Heap heap(npoints, nnbrs);
  const int max_idx = npoints - 1; // internally we need to be 0-indexed
  for (std::size_t i = 0; i < npoints; i++) {
    for (std::size_t j = 0; j < nnbrs; j++) {
      const int k = idx(i, j);
      if (k < 0 || k > max_idx) {
        Rcpp::stop("Bad indexes in input");
      }
      const double d = dist(i, j);
      heap.push(i, d, k, true);
      heap.push(k, d, i, true);
    }
  }
  
  return heap;
}

// transfer data into R Matrices
Rcpp::List heap_to_r(const Heap& heap) {
  const std::size_t npoints = heap.idx.size();
  const std::size_t nnbrs = heap.idx[0].size();
  
  Rcpp::IntegerMatrix idxres(npoints, nnbrs);
  Rcpp::NumericMatrix distres(npoints, nnbrs);
  for (std::size_t i = 0; i < npoints; i++) {
    for (std::size_t j = 0; j < nnbrs; j++) {
      idxres(i, j) = heap.idx[i][j];
      distres(i, j) = heap.dist[i][j];
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("idx") = idxres,
    Rcpp::Named("dist") = distres
  );
}

template <typename Distance,
          typename Rand,
          typename Progress>
Rcpp::List nn_descent_impl(
    Rcpp::NumericMatrix data,
    Rcpp::IntegerMatrix idx,
    Rcpp::NumericMatrix dist,
    const std::size_t max_candidates = 50,
    const std::size_t n_iters = 10,
    const double delta = 0.001,
    const double rho = 0.5,
    bool verbose = false) {
  const std::size_t npoints = idx.nrow();
  const std::size_t nnbrs = idx.ncol();
  
  Heap heap = r_to_heap(idx, dist);
  
  const std::size_t ndim = data.ncol();
  data = Rcpp::transpose(data);
  auto data_vec = Rcpp::as<std::vector<typename Distance::in_type>>(data);
  
  Progress progress;
  Rand rand;
  Distance distance(data_vec, ndim);
  const double tol = delta * nnbrs * npoints;
  
  for (std::size_t n = 0; n < n_iters; n++) {
    if (verbose) {
      progress.iter(n, n_iters, heap);
    }
    
    Heap candidate_neighbors(npoints, max_candidates);
    
    build_candidates<Rand>(heap, candidate_neighbors, npoints, nnbrs);
    
    std::size_t c = 0;
    for (std::size_t i = 0; i < npoints; i++) {
      for (std::size_t j = 0; j < max_candidates; j++) {
        int p = candidate_neighbors.idx[i][j];
        if (p < 0 || rand.unif() < rho) {
          continue;
        }
        
        for (std::size_t k = 0; k < max_candidates; k++) {
          int q = candidate_neighbors.idx[i][k];
          if (q < 0 || (!candidate_neighbors.flags[i][j] &&
              !candidate_neighbors.flags[i][k])) {
              continue;
          }
          double d = distance(p, q);
          c += heap.push(p, d, q, true);
          c += heap.push(q, d, p, true);
        }
      }
      progress.check_interrupt();
    }
    if (static_cast<double>(c) <= tol) {
      if (verbose) {
        progress.converged(c, tol);
      }
      break;
    }
  }
  
  // sort data
  heap.deheap_sort();
  
  return heap_to_r(heap);
}

// [[Rcpp::export]]
Rcpp::List nn_descent(
    Rcpp::NumericMatrix data,
    Rcpp::IntegerMatrix idx,
    Rcpp::NumericMatrix dist,
    const std::string metric = "euclidean",
    const std::size_t max_candidates = 50,
    const std::size_t n_iters = 10,
    const double delta = 0.001,
    const double rho = 0.5,
    bool verbose = false) {
  
  if (metric == "euclidean") {
    return nn_descent_impl<Euclidean<float, float>,
                           RRand,
                           RProgress>
    (data, idx, dist, max_candidates, n_iters, delta, rho, verbose);
  }
  else if (metric == "l2") {
    return nn_descent_impl<L2<float, float>,
                           RRand,
                           RProgress>
    (data, idx, dist, max_candidates, n_iters, delta, rho, verbose);
  }
  else if (metric == "cosine") {
    return nn_descent_impl<Cosine<float, float>,
                           RRand,
                           RProgress>
    (data, idx, dist, max_candidates, n_iters, delta, rho, verbose);
  }
  else if (metric == "manhattan") {
    return nn_descent_impl<Manhattan<float, float>,
                           RRand,
                           RProgress>
    (data, idx, dist, max_candidates, n_iters, delta, rho, verbose);
  }
  else if (metric == "hamming") {
    return nn_descent_impl<Hamming<uint8_t, std::size_t>,
                           RRand,
                           RProgress>
    (data, idx, dist, max_candidates, n_iters, delta, rho, verbose);
  }
  else {
    Rcpp::stop("Bad metric");
  }
}
