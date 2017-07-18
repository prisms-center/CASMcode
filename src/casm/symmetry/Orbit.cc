#include "casm/symmetry/Orbit_impl.hh"
#include "casm/kinetics/PrimPeriodicDiffTransOrbitTraits.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/app/AppIO.hh"
#include "casm/clex/ChemicalReference.hh"
#include "casm/app/DirectoryStructure.hh"

namespace CASM {

  template<>
  std::string _generate_orbit_name(const Kinetics::PrimPeriodicDiffTransOrbit &orbit) {
    return traits<Kinetics::PrimPeriodicDiffTransOrbit>::orbit_type_name + "/" + orbit.id();
  }

  template<>
  std::ostream &write_pos(const Kinetics::PrimPeriodicDiffTransOrbit &orbit, std::ostream &sout) {
    std::vector<Kinetics::PrimPeriodicDiffTransOrbit> container;
    container.push_back(orbit);
    PrototypePrinter<Kinetics::DiffusionTransformation> printer;
    print_clust(container.begin(), container.end(), sout, printer);
    return sout;
  }

  template<>
  void _write_pos(const Kinetics::PrimPeriodicDiffTransOrbit &orbit) {

    const auto &dir = orbit.primclex().dir();
    try {
      fs::create_directories(dir.configuration_dir(orbit.name()));
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in Orbit::write_pos()." << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream file(dir.POS(orbit.name()));
    write_pos(orbit, file);
    return;
  };

}
