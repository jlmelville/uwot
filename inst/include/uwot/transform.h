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

void initialize_by_mean(std::size_t begin, std::size_t end, std::size_t ndim,
                        std::size_t n_neighbors,
                        const std::vector<int> &nn_index,
                        std::size_t n_test_vertices,
                        const std::vector<float> &train_embedding,
                        std::size_t n_train_vertices,
                        std::vector<float> &embedding) {
  std::vector<double> sumc(ndim);
  for (std::size_t i = begin; i < end; i++) {
    std::fill(sumc.begin(), sumc.end(), 0.0);

    for (std::size_t j = 0; j < n_neighbors; j++) {
      std::size_t nbr = nn_index[i + j * n_test_vertices];
      for (std::size_t k = 0; k < ndim; k++) {
        sumc[k] += train_embedding[nbr + k * n_train_vertices];
      }
    }

    for (std::size_t k = 0; k < ndim; k++) {
      embedding[i + k * n_test_vertices] = sumc[k] / n_neighbors;
    }
  }
}

void initialize_by_weighted_mean(std::size_t begin, std::size_t end,
                                 std::size_t ndim, std::size_t n_neighbors,
                                 const std::vector<int> &nn_index,
                                 const std::vector<float> &nn_weights,
                                 std::size_t n_test_vertices,
                                 const std::vector<float> &train_embedding,
                                 std::size_t n_train_vertices,
                                 std::vector<float> &embedding) {
  std::vector<double> sumc(ndim);
  for (std::size_t i = begin; i < end; i++) {
    std::fill(sumc.begin(), sumc.end(), 0.0);

    double sumw = 0.0;

    for (std::size_t j = 0; j < n_neighbors; j++) {
      std::size_t nbr = nn_index[i + j * n_test_vertices];
      float w = nn_weights[i + j * n_test_vertices];
      sumw += w;
      for (std::size_t k = 0; k < ndim; k++) {
        sumc[k] += train_embedding[nbr + k * n_train_vertices] * w;
      }
    }

    for (std::size_t k = 0; k < ndim; k++) {
      embedding[i + k * n_test_vertices] = sumc[k] / sumw;
    }
  }
}

} // namespace uwot

#endif // UWOT_TRANSFORM_H
