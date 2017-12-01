#include "casm/crystallography/BasicStructure_impl.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/basis_set/DoF.hh"

namespace CASM {

  template class BasicStructure<Site>;
  template Index BasicStructure<Site>::find<Coordinate>(const Coordinate &bsite) const;
  template Index BasicStructure<Site>::find<Site>(const Site &bsite) const;

}
