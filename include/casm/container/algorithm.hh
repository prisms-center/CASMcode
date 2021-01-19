#ifndef CONTAINER_ALGO_HH
#define CONTAINER_ALGO_HH

#include <vector>

namespace CASM {

template <typename T>
std::vector<T> sequence(T first, T last) {
  std::vector<T> seq;
  if (0 < last - first) seq.reserve(last - first);
  while (first <= last) {
    seq.push_back(first);
    first += 1;
  }
  return seq;
}

template <typename T>
std::vector<T> sequence(T first, T inc, T last) {
  std::vector<T> seq;
  // Decreasing case
  if (first + inc < first) {
    while (last <= first) {
      seq.push_back(first);
      first += inc;
    }
  }
  // Increasing case
  else if (first < first + inc) {
    while (first <= last) {
      seq.push_back(first);
      first += inc;
    }
  } else {  // increment=0 case
    if (first == last) {
      seq.push_back(first);
    }
  }
  return seq;
}
}  // namespace CASM

#endif
