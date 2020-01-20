#ifndef BASICSTRUCTURETOOLS_HH
#define BASICSTRUCTURETOOLS_HH

#include "casm/global/definitions.hh"
#include <vector>

namespace CASM {
  namespace xtal {
    class Site;
    class Coordinate;

    /// return basis index of site that matches test_coord, if it is in basis
    /// otherwise, returns basis.size()
    Index find_index(const std::vector<Site> &basis, const Site &test_site);

    /// return basis index of site that matches test_site+shift, if it is in basis
    /// otherwise, returns basis.size()
    Index find_index(const std::vector<Site> &basis, const Site &test_site, const Coordinate &shift);

  } // namespace xtal
} // namespace CASM

#endif
