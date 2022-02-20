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

namespace uwot {

inline auto sum_nbrs(std::size_t i, const std::vector<float> &train_embedding,
                     const std::vector<int> &nn_index, std::size_t n_neighbors,
                     std::size_t ndim, std::vector<double> &sumc) -> void {
  for (std::size_t j = 0; j < n_neighbors; j++) {
    auto nbr = nn_index[i * n_neighbors + j];
    for (std::size_t k = 0; k < ndim; k++) {
      sumc[k] += train_embedding[ndim * nbr + k];
    }
  }
}

inline auto sum_nbrs_weighted(std::size_t i,
                              const std::vector<float> &train_embedding,
                              const std::vector<int> &nn_index,
                              std::size_t n_neighbors, std::size_t ndim,
                              const std::vector<float> &nn_weights,
                              std::vector<double> &sumc, double &sumw) -> void {
  std::size_t i_nbrs = i * n_neighbors;
  for (std::size_t j = 0; j < n_neighbors; j++) {
    auto w = nn_weights[i_nbrs + j];
    sumw += w;
    auto nbr = nn_index[i_nbrs + j];
    for (std::size_t k = 0; k < ndim; k++) {
      sumc[k] += train_embedding[ndim * nbr + k] * w;
    }
  }
}

void init_by_mean(std::size_t begin, std::size_t end, std::size_t ndim,
                  std::size_t n_neighbors, const std::vector<int> &nn_index,
                  const std::vector<float> &nn_weights,
                  std::size_t n_test_vertices,
                  const std::vector<float> &train_embedding,
                  std::size_t n_train_vertices, std::vector<float> &embedding) {
  bool weighted = nn_weights.size() > 0;

  std::vector<double> sumc(ndim);
  for (std::size_t i = begin; i < end; i++) {
    std::fill(sumc.begin(), sumc.end(), 0.0);

    double sumw = 0.0;
    // cost of checking this boolean N times is not going to be a bottleneck
    if (weighted) {
      sum_nbrs_weighted(i, train_embedding, nn_index, n_neighbors, ndim,
                        nn_weights, sumc, sumw);
    } else {
      sumw = static_cast<double>(n_neighbors);
      sum_nbrs(i, train_embedding, nn_index, n_neighbors, ndim, sumc);
    }

    for (std::size_t k = 0; k < ndim; k++) {
      embedding[ndim * i + k] = sumc[k] / sumw;
    }
  }
}

} // namespace uwot

#endif // UWOT_TRANSFORM_H
