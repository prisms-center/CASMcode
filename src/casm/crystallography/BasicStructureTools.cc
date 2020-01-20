#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {
  namespace xtal {
    Index find_index(const std::vector<Site> &basis, const Site &test_site) {
      for(Index i = 0; i < basis.size(); i++) {
        if(basis[i].compare(test_site)) {
          return i;
        }
      }
      return basis.size();
    }

    Index find_index(const std::vector<Site> &basis, const Site &test_site, const Coordinate &shift) {
      for(Index i = 0; i < basis.size(); i++) {
        if(basis[i].compare(test_site, shift)) {
          return i;
        }
      }
      return basis.size();
    }

  } // namespace xtal
} // namespace CASM
