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

#include <cmath>
#include <limits>
#include <vector>

#ifndef UWOT_SUPERVISED_H
#define UWOT_SUPERVISED_H

namespace uwot {

void fast_intersection(const std::vector<int> &rows,
                       const std::vector<int> &cols,
                       std::vector<double> &values,
                       const std::vector<int> &target,
                       double unknown_dist = 1.0, double far_dist = 5.0,
                       int na = (std::numeric_limits<int>::min)()

) {

  double ex_unknown = std::exp(-unknown_dist);
  double ex_far = std::exp(-far_dist);

  auto len = values.size();

  for (std::size_t nz = 0; nz < len; ++nz) {
    auto i = rows[nz];
    auto j = cols[nz];
    if (target[i] == na || target[j] == na) {
      values[nz] = values[nz] * ex_unknown;
    } else if (target[i] != target[j]) {
      values[nz] = values[nz] * ex_far;
    }
  }
}

} // namespace uwot
#endif // UWOT_SUPERVISED_H
