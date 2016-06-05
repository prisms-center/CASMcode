#include "casm/crystallography/BasicStructure_impl.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {

  template class BasicStructure<Site>;
  template Index BasicStructure<Site>::find<Coordinate>(const Coordinate &bsite, double tol) const;
  template Index BasicStructure<Site>::find<Site>(const Site &bsite, double tol) const;

}