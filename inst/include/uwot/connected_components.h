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

// Translated from the Python source code of:
//   scipy.sparse.csgraph.connected_components
//   Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
//   License: BSD, (C) 2012

#ifndef UWOT_CONNECTED_COMPONENTS_H
#define UWOT_CONNECTED_COMPONENTS_H

namespace uwot {

auto connected_components_undirected(std::size_t N,
                                     const std::vector<int> &indices1,
                                     const std::vector<int> &indptr1,
                                     const std::vector<int> &indices2,
                                     const std::vector<int> &indptr2)
    -> std::pair<unsigned int, std::vector<int>> {
  constexpr int VOID = -1;
  constexpr int END = -2;
  std::vector<int> labels(N, VOID);
  std::vector<int> SS(labels);
  unsigned int label = 0;
  auto SS_head = END;
  for (std::size_t v = 0; v < N; ++v) {
    auto vv = v;
    if (labels[vv] == VOID) {
      SS_head = vv;
      SS[vv] = END;
      while (SS_head != END) {
        vv = SS_head;
        SS_head = SS[vv];
        labels[vv] = label;
        for (auto jj = indptr1[vv]; jj < indptr1[vv + 1]; ++jj) {
          auto ww = indices1[jj];
          if (SS[ww] == VOID) {
            SS[ww] = SS_head;
            SS_head = ww;
          }
        }
        for (auto jj = indptr2[vv]; jj < indptr2[vv + 1]; ++jj) {
          auto ww = indices2[jj];
          if (SS[ww] == VOID) {
            SS[ww] = SS_head;
            SS_head = ww;
          }
        }
      }
      ++label;
    }
  }
  return {label, labels};
}
} // namespace uwot

#endif // UWOT_CONNECTED_COMPONENTS_H