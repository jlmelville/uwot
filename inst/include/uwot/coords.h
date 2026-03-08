// BSD 2-Clause License
//
// Copyright 2021 James Melville
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

#ifndef UWOT_COORDS_H
#define UWOT_COORDS_H

#include <optional>
#include <utility>
#include <vector>

namespace uwot {

// For normal UMAP there is no separate tail embedding, so the head embedding
// is used anywhere a tail embedding is requested.
// When updating new values, tail_embedding stores the coordinates to optimize.
struct Coords {
  std::vector<float> head_embedding;
  std::optional<std::vector<float>> tail_embedding;

  explicit Coords(std::vector<float> head_embedding)
      : head_embedding(std::move(head_embedding)), tail_embedding(std::nullopt) {}

  Coords(std::vector<float> head_embedding, std::vector<float> tail_embedding)
      : head_embedding(std::move(head_embedding)),
        tail_embedding(std::move(tail_embedding)) {}

  auto get_tail_embedding() -> std::vector<float> & {
    return tail_embedding ? *tail_embedding : head_embedding;
  }

  auto get_head_embedding() -> std::vector<float> & { return head_embedding; }
};

} // namespace uwot

#endif
