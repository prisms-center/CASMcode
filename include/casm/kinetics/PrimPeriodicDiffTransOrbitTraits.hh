#ifndef CASM_PrimPeriodicDiffTransOrbitTraits
#define CASM_PrimPeriodicDiffTransOrbitTraits

#include <string>
#include <vector>
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/symmetry/Orbit.hh"
namespace CASM {
  template<typename T> struct traits;

  namespace Kinetics {
    class DiffusionTransformation;
    class PrimPeriodicDiffTransSymCompare;
  }

  template<typename _Element, typename _SymCompare>
  class Orbit;

  typedef Orbit<Kinetics::DiffusionTransformation, Kinetics::PrimPeriodicDiffTransSymCompare> PrimPeriodicDiffTransOrbit;

  template<>
  struct traits<CASM::PrimPeriodicDiffTransOrbit> {
    static const std::string name;
    static const std::string short_name;
    static const std::string orbit_type_name;
    static bool name_compare(std::string A, std::string B);
  };
}

#endif