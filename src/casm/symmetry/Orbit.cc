#include "casm/symmetry/Orbit_impl.hh"
#include "casm/kinetics/PrimPeriodicDiffTransOrbitTraits.hh"
#include "casm/kinetics/DiffusionTransformation.hh"

namespace CASM {

  template<>
  std::string _generate_orbit_name(const Kinetics::PrimPeriodicDiffTransOrbit &orbit) {
    return traits<Kinetics::PrimPeriodicDiffTransOrbit>::orbit_type_name + "." + orbit.id();
  }

}
