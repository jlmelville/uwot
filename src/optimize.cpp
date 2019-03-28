//  UWOT -- An R package for dimensionality reduction using UMAP
//
//  Copyright (C) 2018 James Melville
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

#include <limits>
#include <memory>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "gradient.h"
#include "sampler.h"
#include "tauprng.h"

// Function to decide whether to move both vertices in an edge
// Default empty version does nothing: used in umap_transform when
// some of the vertices should be held fixed
template <bool DoMoveVertex = false>
void move_other_vertex(
    std::vector<double>& embedding,
    const double grad_d,
    const std::size_t i,
    const std::size_t nrj) {
}

// Specialization to move the vertex: used in umap when both
// vertices in an edge should be moved
template <>
void move_other_vertex<true>(
    std::vector<double>& embedding,
    const double grad_d,
    const std::size_t i,
    const std::size_t nrj) {
  embedding[nrj + i] -= grad_d;
}

const double clamp(const double v, const double lo, const double hi) {
  const double t = v < lo ? lo : v;
  return t > hi ? hi : t;
}

// Gradient: the type of gradient used in the optimization
// DoMoveVertex: true if both ends of a positive edge should be updated
template <typename Gradient,
          bool DoMoveVertex = true,
          typename Rand = pcg_prng>
struct SgdWorker : public RcppParallel::Worker {
  int n; // epoch counter
  double alpha;
  const Gradient gradient;
  const std::vector<unsigned int> positive_head;
  const std::vector<unsigned int> positive_tail;
  Sampler sampler;
  std::vector<double>& head_embedding;
  std::vector<double>& tail_embedding;
  const std::size_t ndim;
  const std::size_t head_nvert;
  const std::size_t tail_nvert;
  tthread::mutex mutex;
  const double dist_eps;
  
  SgdWorker(
    const Gradient& gradient,
    const std::vector<unsigned int>& positive_head,
    const std::vector<unsigned int>& positive_tail,
    Sampler& sampler,
    std::vector<double>& head_embedding,
    std::vector<double>& tail_embedding,
    const std::size_t ndim) :
    
    n(0), alpha(0.0), gradient(gradient),
    positive_head(positive_head), positive_tail(positive_tail),
    
    sampler(sampler),
    
    head_embedding(head_embedding),
    tail_embedding(tail_embedding),
    ndim(ndim), 
    head_nvert(head_embedding.size() / ndim), 
    tail_nvert(tail_embedding.size() / ndim),
    dist_eps(std::numeric_limits<double>::epsilon())
  { }
  
  void  operator()(std::size_t begin, std::size_t end) {
    std::unique_ptr<Rand> prng(nullptr);
   
    // Each window gets its own PRNG state, to prevent locking inside the loop.
    {
      tthread::lock_guard<tthread::mutex> guard(mutex);
      prng.reset(new Rand());
      
    }

    std::vector<double> dys(ndim);
    for (std::size_t i = begin; i < end; i++) {
      if (!sampler.is_sample_edge(i, n)) {
        continue;
      }
      
      const std::size_t dj = ndim * positive_head[i];
      const std::size_t dk = ndim * positive_tail[i];

      double dist_squared = 0.0;
      for (std::size_t d = 0; d < ndim; d++) {
        const double diff = head_embedding[dj + d] - tail_embedding[dk + d];
        dys[d] = diff;
        dist_squared += diff * diff;
      }
      dist_squared = std::max(dist_eps, dist_squared);
      
      const double grad_coeff = gradient.grad_attr(dist_squared);
      for (std::size_t d = 0; d < ndim; d++) {
        const double grad_d = alpha * clamp(grad_coeff * dys[d], 
                                            Gradient::clamp_lo,
                                            Gradient::clamp_hi);
        head_embedding[dj + d] += grad_d;
        move_other_vertex<DoMoveVertex>(tail_embedding, grad_d, d, dk);
      }

      const std::size_t n_neg_samples = sampler.get_num_neg_samples(i, n);
      for (std::size_t p = 0; p < n_neg_samples; p++) {
        const std::size_t dkn = (*prng)(tail_nvert) * ndim;
        if (dj == dkn) {
          continue;
        }

        double dist_squared = 0.0;
        for (std::size_t d = 0; d < ndim; d++) {
          const double diff = head_embedding[dj + d] - tail_embedding[dkn + d];
          dys[d] = diff;
          dist_squared += diff * diff;
        }
        dist_squared = std::max(dist_eps, dist_squared);
        
        const double grad_coeff = gradient.grad_rep(dist_squared);
        for (std::size_t d = 0; d < ndim; d++) {
          const double grad_d = alpha * clamp(grad_coeff * dys[d], 
                                              Gradient::clamp_lo,
                                              Gradient::clamp_hi);
          head_embedding[dj + d] += grad_d;
        }
      }
      
      sampler.next_sample(i, n_neg_samples);
    }
  }
  
  void set_n(int n) {
    this->n = n;
  }
  
  void set_alpha(double alpha) {
    this->alpha = alpha;
  }
};


template<typename T, bool DoMove = true, typename Rand = pcg_prng>
std::vector<double> optimize_layout(
    const T& gradient,
    std::vector<double>& head_embedding,
    std::vector<double>& tail_embedding,
    const std::vector<unsigned int>& positive_head,
    const std::vector<unsigned int>& positive_tail,
    unsigned int n_epochs, 
    unsigned int n_vertices,
    const std::vector<double>& epochs_per_sample,
    double initial_alpha,
    double negative_sample_rate,
    bool parallelize = true,
    std::size_t grain_size = 1,
    bool verbose = false) 
{
  Sampler sampler(epochs_per_sample, negative_sample_rate);
  
  SgdWorker<T, DoMove, Rand> worker(gradient, 
                              positive_head, positive_tail,
                              sampler,
                              head_embedding, tail_embedding,
                              head_embedding.size() / n_vertices);
  
  
  Progress progress(n_epochs, verbose);
  const auto n_epochs_per_sample = epochs_per_sample.size();
  double alpha = initial_alpha;
  
  for (auto n = 0U; n < n_epochs; n++) {
    worker.set_alpha(alpha);
    worker.set_n(n);
    
    if (parallelize) {
      RcppParallel::parallelFor(0, n_epochs_per_sample, worker, grain_size);
    }
    else {
      worker(0, n_epochs_per_sample);
    }
    alpha = initial_alpha * (1.0 - (double(n) / double(n_epochs)));
    
    if (Progress::check_abort()) {
      return head_embedding;
    }
    if (verbose) {
      progress.increment();
    }
  }
  return head_embedding;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix optimize_layout_umap(
    Rcpp::NumericMatrix head_embedding,
    Rcpp::Nullable<Rcpp::NumericMatrix> tail_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail,
    unsigned int n_epochs, 
    unsigned int n_vertices,
    const std::vector<double> epochs_per_sample,
    double a, 
    double b,
    double gamma, 
    double initial_alpha,
    double negative_sample_rate,
    bool approx_pow,
    bool pcg_rand = true,
    bool parallelize = true,
    std::size_t grain_size = 1,
    bool move_other = true,
    bool verbose = false) 
{
  // For normal UMAP, tail_embedding is NULL and we want to pass 
  // a shallow copy of head_embedding as tail_embedding.
  // When updating new values, tail_embedding is the new coordinate to optimize
  // and gets passed as normal.
  auto head_vec = Rcpp::as<std::vector<double>>(head_embedding);
  std::vector<double>* tail_vec_ptr = nullptr;
  bool delete_tail_ptr = false;
  if (tail_embedding.isNull()) {
    tail_vec_ptr = &head_vec;
  }
  else {
    tail_vec_ptr = new std::vector<double>(
      Rcpp::as<std::vector<double>>(tail_embedding));
    delete_tail_ptr = true;
  }
  
  std::vector<double> result;
  if (approx_pow) {
    const apumap_gradient gradient(a, b, gamma);
    if (move_other) {
      if (pcg_rand) {
        result = optimize_layout<apumap_gradient, true, pcg_prng>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, parallelize, grain_size, verbose);
      }
      else {
        result = optimize_layout<apumap_gradient, true, tau_prng>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, parallelize, grain_size, verbose);
      }
    }
    else {
      if (pcg_rand) {
        result = optimize_layout<apumap_gradient, false, pcg_prng>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, parallelize, grain_size, verbose);
      }
      else {
        result = optimize_layout<apumap_gradient, false, tau_prng>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, parallelize, grain_size, verbose);
      }
    }
  }
  else {
    const umap_gradient gradient(a, b, gamma);
    if (move_other) {
      if (pcg_rand) {
        result = optimize_layout<umap_gradient, true, pcg_prng>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, parallelize, grain_size, verbose);
      }
      else {
        result = optimize_layout<umap_gradient, true, tau_prng>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs, n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, parallelize, grain_size, verbose);
      }
    }
    else {
      if (pcg_rand) {
        result = optimize_layout<umap_gradient, false, pcg_prng>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs,n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, parallelize, grain_size, verbose);
      }
      else {
        result = optimize_layout<umap_gradient, false, tau_prng>(
          gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail,
          n_epochs,n_vertices, epochs_per_sample, initial_alpha,
          negative_sample_rate, parallelize, grain_size, verbose);
      }
    }
  }
  
  if (delete_tail_ptr) {
    delete(tail_vec_ptr);
  }
  
  return Rcpp::NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                             result.begin());
}



// [[Rcpp::export]]
Rcpp::NumericMatrix optimize_layout_tumap(
    Rcpp::NumericMatrix head_embedding,
    Rcpp::Nullable<Rcpp::NumericMatrix> tail_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail,
    unsigned int n_epochs, 
    unsigned int n_vertices,
    const std::vector<double> epochs_per_sample,
    double initial_alpha,
    double negative_sample_rate,
    bool pcg_rand = true,
    bool parallelize = true,
    std::size_t grain_size = 1,
    bool move_other = true,
    bool verbose = false) 
{
  const tumap_gradient gradient;
  auto head_vec = Rcpp::as<std::vector<double>>(head_embedding);
  std::vector<double>* tail_vec_ptr = nullptr;
  bool delete_tail_ptr = false;
  if (tail_embedding.isNull()) {
    tail_vec_ptr = &head_vec;
  }
  else {
    tail_vec_ptr = new std::vector<double>(
      Rcpp::as<std::vector<double>>(tail_embedding));
    delete_tail_ptr = true;
  }
  
  std::vector<double> result;
  
  if (move_other) {
    if (pcg_rand) {
      result = optimize_layout<tumap_gradient, true, pcg_prng>(
        gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        parallelize, grain_size, verbose);
    }
    else {
      result = optimize_layout<tumap_gradient, true, tau_prng>(
        gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        parallelize, grain_size, verbose);      
    }
  }
  else {
    if (pcg_rand) {
      result = optimize_layout<tumap_gradient, false, pcg_prng>(
        gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        parallelize, grain_size, verbose);
    }
    else {
      result = optimize_layout<tumap_gradient, false, tau_prng>(
        gradient, head_vec, *tail_vec_ptr, positive_head, positive_tail, n_epochs,
        n_vertices, epochs_per_sample, initial_alpha, negative_sample_rate,
        parallelize, grain_size, verbose);      
    }
  }
  
  if (delete_tail_ptr) {
    delete(tail_vec_ptr);
  }
  
  return Rcpp::NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                             result.begin());
}

// [[Rcpp::export]]
Rcpp::NumericMatrix optimize_layout_largevis(
    Rcpp::NumericMatrix head_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail,
    unsigned int n_epochs, 
    unsigned int n_vertices,
    const std::vector<double> epochs_per_sample,
    double gamma, double initial_alpha,
    double negative_sample_rate,
    bool pcg_rand = true,
    bool parallelize = true,
    std::size_t grain_size = 1,
    bool verbose = false) 
{
  // We don't support adding extra points for LargeVis, so this is much simpler
  // than the UMAP case
  const largevis_gradient gradient(gamma);
  auto head_vec = Rcpp::as<std::vector<double>>(head_embedding);
  
  std::vector<double> result;
  
  if (pcg_rand) {
    result = optimize_layout<largevis_gradient, true, pcg_prng>(
      gradient, head_vec, head_vec, positive_head, positive_tail,
      n_epochs, n_vertices, epochs_per_sample, initial_alpha,
      negative_sample_rate, parallelize, grain_size, verbose);
  }
  else {
    result = optimize_layout<largevis_gradient, true, tau_prng>(
      gradient, head_vec, head_vec, positive_head, positive_tail,
      n_epochs, n_vertices, epochs_per_sample, initial_alpha,
      negative_sample_rate, parallelize, grain_size, verbose);
  }
  
  return Rcpp::NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                             result.begin());
}
