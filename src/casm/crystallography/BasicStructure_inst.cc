#include "casm/crystallography/BasicStructure_impl.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {

  template class BasicStructure<Site>;
  template UnitCellCoord BasicStructure<Site>::get_unit_cell_coord<Coordinate>(const Coordinate &bsite, double tol) const;
  template UnitCellCoord BasicStructure<Site>::get_unit_cell_coord<Site>(const Site &bsite, double tol) const;
  template Index BasicStructure<Site>::find<Coordinate>(const Coordinate &bsite, double tol) const;
  template Index BasicStructure<Site>::find<Site>(const Site &bsite, double tol) const;

}