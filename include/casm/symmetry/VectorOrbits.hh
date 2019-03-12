#ifndef CASM_VectorOrbits
#define CASM_VectorOrbits

#include <vector>
#include <functional>
#include <utility>
#include <iostream>

#include "casm/symmetry/SymGroup.hh"

namespace CASM {

  /* -- Vector Orbit generating function declarations ------------------------------------- */

  /// \brief Generate Orbit<IntegralVector> using OrbitBranchSpecs
  template<typename GeneratorIterator, typename SymCompareType, typename OrbitOutputIterator>
  OrbitOutputIterator make_orbits(
    GeneratorIterator gen_begin,
    GeneratorIterator gen_end,
    const SymGroup &generating_group,
    const SymCompareType &sym_compare,
    OrbitOutputIterator result,
    std::ostream &status) {

    using OrbitType = OrbitOutpuIterator::container_type::value_type;
    for(; get_begin != gen_end; ++gen_begin) {
      *(result++) = OrbitType(*gen_begin, generating_group, sym_compare);
    }
    return result;
  }


}

#endif
