#ifndef CASM_clex_CompositionAxes_impl
#define CASM_clex_CompositionAxes_impl

#include "casm/clex/CompositionAxes.hh"

namespace CASM {

template <typename IterType>
void CompositionAxes::insert_enumerated(IterType begin, IterType end) {
  int i = 0;
  for (; begin != end; ++begin, ++i) {
    while (all_axes.count(std::to_string(i))) ++i;
    all_axes[std::to_string(i)] = *begin;
    enumerated.insert(std::to_string(i));
  }
}

}  // namespace CASM

#endif
