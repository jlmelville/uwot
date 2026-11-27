#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>

// hnswlib's saveIndex writes a native header, a fixed-width level-zero table,
// then one length-prefixed upper-level link list per item (including empty
// lists). Check the same boundaries as loadIndex without allocating an index.
// [[Rcpp::export(rng = false)]]
bool hnsw_index_complete(const std::string &path) {
  std::ifstream input(path, std::ios::binary | std::ios::ate);
  constexpr std::streamoff header_size = 10 * sizeof(std::size_t) +
                                         sizeof(int) + sizeof(unsigned int) +
                                         sizeof(double);
  const std::streamoff end = input.tellg();
  if (end < header_size) {
    return false;
  }

  // The first four fields are level-zero offset, capacity, count and stride.
  input.seekg(2 * sizeof(std::size_t));
  std::size_t count = 0, stride = 0;
  input.read(reinterpret_cast<char *>(&count), sizeof(count));
  input.read(reinterpret_cast<char *>(&stride), sizeof(stride));
  const auto payload = static_cast<std::uintmax_t>(end - header_size);
  if (!input || count == 0 || stride == 0 || stride > payload / count) {
    return false;
  }
  const auto base_bytes = static_cast<std::uintmax_t>(count) * stride;
  if (count > (payload - base_bytes) / sizeof(unsigned int)) {
    return false;
  }
  input.seekg(header_size + static_cast<std::streamoff>(base_bytes));
  for (std::size_t i = 0; i < count; ++i) {
    unsigned int links = 0;
    if (!input.read(reinterpret_cast<char *>(&links), sizeof(links))) {
      return false;
    }
    if (links != 0) {
      input.ignore(links);
    }
  }
  return input.tellg() == end;
}
