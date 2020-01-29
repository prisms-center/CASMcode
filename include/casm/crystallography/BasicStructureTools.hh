#ifndef BASICSTRUCTURETOOLS_HH
#define BASICSTRUCTURETOOLS_HH

#include "casm/global/definitions.hh"
#include <vector>

namespace CASM {
  namespace xtal {
    class Site;
    class Coordinate;
    class BasicStructure;
    class SymOp;

    /// return basis index of site that matches test_coord, if it is in basis
    /// otherwise, returns basis.size()
    Index find_index(const std::vector<Site> &basis, const Site &test_site);

    /// return basis index of site that matches test_site+shift, if it is in basis
    /// otherwise, returns basis.size()
    Index find_index(const std::vector<Site> &basis, const Site &test_site, const Coordinate &shift);

    /// Returns true if the structure describes a crystal primitive cell
    /// i.e., no translation smaller than a lattice vector can map the structure onto itself
    bool is_primitive(const BasicStructure &struc);

    /// Returns the smallest possible tiling unit of the given structure
    BasicStructure make_primitive(const BasicStructure &non_primitive_struc);

    /// Create the factor group of the given structure. If the structure has not degrees of freedom
    /// affected by time reversal, time reversal is ignored. Otherwise symmetry operations are checked
    /// for time reversal
    std::vector<SymOp> make_factor_group(const BasicStructure &struc);

  } // namespace xtal
} // namespace CASM

#endif
