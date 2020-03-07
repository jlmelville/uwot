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

// Three-component combined Tausworthe "taus88" PRNG from L'Ecuyer 1996.

#ifndef UWOT_TAUPRNG_H
#define UWOT_TAUPRNG_H

namespace uwot {

struct tau_prng {
  uint64_t state0;
  uint64_t state1; // technically this needs to always be > 7
  uint64_t state2; // and this should be > 15

  static constexpr uint64_t MAGIC0 = static_cast<uint64_t>(4294967294);
  static constexpr uint64_t MAGIC1 = static_cast<uint64_t>(4294967288);
  static constexpr uint64_t MAGIC2 = static_cast<uint64_t>(4294967280);

  tau_prng(uint64_t state0, uint64_t state1, uint64_t state2)
      : state0(state0), state1(state1 > 7 ? state1 : 8),
        state2(state2 > 15 ? state2 : 16) {}

  auto operator()() -> int32_t {
    state0 = (((state0 & MAGIC0) << 12) & 0xffffffff) ^
             ((((state0 << 13) & 0xffffffff) ^ state0) >> 19);
    state1 = (((state1 & MAGIC1) << 4) & 0xffffffff) ^
             ((((state1 << 2) & 0xffffffff) ^ state1) >> 25);
    state2 = (((state2 & MAGIC2) << 17) & 0xffffffff) ^
             ((((state2 << 3) & 0xffffffff) ^ state2) >> 11);

    return state0 ^ state1 ^ state2;
  }

  // return a value in (0, n]
  auto operator()(std::size_t n) -> std::size_t {
    std::size_t result = (*this)() % n;
    return result;
  }
};

} // namespace uwot

#endif // UWOT_TAUPRNG_H
