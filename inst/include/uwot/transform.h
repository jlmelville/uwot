// BSD 2-Clause License
//
// Copyright 2020 James Melville
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// OF SUCH DAMAGE.

#ifndef UWOT_TRANSFORM_H
#define UWOT_TRANSFORM_H

#include <algorithm>
#include <cstddef>
#include <vector>

namespace uwot {

template <typename Weight>
inline double accumulate_neighbors(std::size_t first, std::size_t last,
                                   const std::vector<float> &training,
                                   const std::vector<int> &indices,
                                   std::size_t ndim, const Weight *weights,
                                   std::vector<double> &sum,
                                   double scale = 1.0) {
  double total = 0;
  for (auto j = first; j < last; ++j) {
    const double weight = weights ? weights[j] / scale : 1.0;
    total += weight;
    for (std::size_t k = 0; k < ndim; ++k) {
      sum[k] += training[ndim * indices[j] + k] * weight;
    }
  }
  return total;
}

inline double sum_neighbors(std::size_t first, std::size_t last,
                            const std::vector<float> &training,
                            const std::vector<int> &indices, std::size_t ndim,
                            const float *weights, std::vector<double> &sum) {
  return accumulate_neighbors(first, last, training, indices, ndim, weights,
                              sum);
}

inline void sum_nbrs(std::size_t i, const std::vector<float> &training,
                     const std::vector<int> &indices, std::size_t k,
                     std::size_t ndim, std::vector<double> &sum) {
  sum_neighbors(i * k, (i + 1) * k, training, indices, ndim, nullptr, sum);
}

inline void sum_nbrs_weighted(std::size_t i, const std::vector<float> &training,
                              const std::vector<int> &indices, std::size_t k,
                              std::size_t ndim,
                              const std::vector<float> &weights,
                              std::vector<double> &sum, double &total) {
  total += sum_neighbors(i * k, (i + 1) * k, training, indices, ndim,
                         weights.data(), sum);
}

template <typename Weight = float>
inline void init_by_mean(std::size_t begin, std::size_t end, std::size_t ndim,
                         const std::vector<std::size_t> &ptr,
                         const std::vector<int> &nn_index,
                         const std::vector<Weight> &nn_weights,
                         const std::vector<float> &train_embedding,
                         std::vector<float> &embedding) {
  std::vector<double> sum(ndim);
  for (auto i = begin; i < end; ++i) {
    const auto first = ptr.size() == 1 ? i * ptr[0] : ptr[i];
    const auto last = ptr.size() == 1 ? first + ptr[0] : ptr[i + 1];
    std::fill(sum.begin(), sum.end(), 0.0);
    // Scale within each query so even subnormal weights retain their ratios.
    // If every weight underflowed, use the ordinary neighbor mean.
    const double scale = nn_weights.empty()
                             ? 0.0
                             : *std::max_element(nn_weights.begin() + first,
                                                 nn_weights.begin() + last);
    const double total = accumulate_neighbors(
        first, last, train_embedding, nn_index, ndim,
        scale > 0 ? nn_weights.data() : nullptr, sum, scale);
    for (std::size_t k = 0; k < ndim; ++k) {
      embedding[ndim * i + k] = sum[k] / total;
    }
  }
}

inline void init_by_mean(std::size_t begin, std::size_t end, std::size_t ndim,
                         std::size_t n_neighbors,
                         const std::vector<int> &nn_index,
                         const std::vector<float> &nn_weights,
                         std::size_t /* n_test_vertices */,
                         const std::vector<float> &train_embedding,
                         std::size_t /* n_train_vertices */,
                         std::vector<float> &embedding) {
  init_by_mean(begin, end, ndim, std::vector<std::size_t>{n_neighbors},
               nn_index, nn_weights, train_embedding, embedding);
}

} // namespace uwot

#endif // UWOT_TRANSFORM_H
