#ifndef CASM_PrimPeriodicDiffTransOrbitTraits
#define CASM_PrimPeriodicDiffTransOrbitTraits

#include <string>
#include "casm/kinetics/DiffusionTransformationTraits.hh"

namespace CASM {

  template<typename T> struct traits;

  template<>
  struct traits<PrimPeriodicDiffTransOrbit> {
    static const std::string name;
    static const std::string short_name;
    static const std::string orbit_type_name;
    static bool name_compare(std::string A, std::string B);
  };
}

#endif
