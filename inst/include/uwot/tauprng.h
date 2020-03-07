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

  int32_t operator()() {
    state0 = (((state0 & MAGIC0) << 12) & 0xffffffff) ^
             ((((state0 << 13) & 0xffffffff) ^ state0) >> 19);
    state1 = (((state1 & MAGIC1) << 4) & 0xffffffff) ^
             ((((state1 << 2) & 0xffffffff) ^ state1) >> 25);
    state2 = (((state2 & MAGIC2) << 17) & 0xffffffff) ^
             ((((state2 << 3) & 0xffffffff) ^ state2) >> 11);

    return state0 ^ state1 ^ state2;
  }

  // return a value in (0, n]
  std::size_t operator()(std::size_t n) {
    std::size_t result = (*this)() % n;
    return result;
  }
};

} // namespace uwot

#endif // UWOT_TAUPRNG_H
