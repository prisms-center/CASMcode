#ifndef CASM_PrimPeriodicDiffTransOrbitTraits
#define CASM_PrimPeriodicDiffTransOrbitTraits

#include <string>

namespace CASM {
  template<typename T> struct traits;
  template<typename T> class PrimPeriodicSymCompare;
  template<typename _Element, typename _SymCompare> class Orbit;

  namespace Kinetics {
    class DiffusionTransformation;
  }
  template<> class PrimPeriodicSymCompare<Kinetics::DiffusionTransformation>;

  namespace Kinetics {
    typedef PrimPeriodicSymCompare<Kinetics::DiffusionTransformation> PrimPeriodicDiffTransSymCompare;
    typedef Orbit <
    Kinetics::DiffusionTransformation,
             Kinetics::PrimPeriodicDiffTransSymCompare > PrimPeriodicDiffTransOrbit;
  }

  template<>
  struct traits<Kinetics::PrimPeriodicDiffTransOrbit> {
    static const std::string name;
    static const std::string short_name;
    static const std::string orbit_type_name;
    static bool name_compare(std::string A, std::string B);
  };
}

#endif
