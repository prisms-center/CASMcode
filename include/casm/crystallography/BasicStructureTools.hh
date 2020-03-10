#ifndef BASICSTRUCTURETOOLS_HH
#define BASICSTRUCTURETOOLS_HH

#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"
#include <vector>

namespace CASM {
  namespace xtal {
    class Site;
    class Coordinate;
    class BasicStructure;
    class SymOp;

    /// return basis index of site that matches test_site, if it is in basis
    /// otherwise, returns the size of the basis
    Index find_index(const std::vector<Site> &basis, const Site &test_site);

    /// return basis index of site that matches test_site+shift, if it is in basis
    /// otherwise, returns the size of the basis
    Index find_index_with_translation(const std::vector<Site> &basis, const Site &test_site, const Coordinate &shift);

    /// Returns true if the structure describes a crystal primitive cell
    /// i.e., no translation smaller than a lattice vector can map the structure onto itself
    bool is_primitive(const BasicStructure &struc);

    /// Returns the smallest possible tiling unit of the given structure
    BasicStructure make_primitive(const BasicStructure &non_primitive_struc);

    /// Create the factor group of the given structure. If the structure has no degrees of freedom
    /// affected by time reversal, time reversal is ignored. Otherwise symmetry operations are checked
    /// for time reversal
    std::vector<SymOp> make_factor_group(const BasicStructure &struc);

    /// Given a symmetry group, the basis of the structure will have
    /// each operation applied to it. The resulting set of basis
    /// from performing these operations will be averaged out,
    /// yielding a new average basis.
    BasicStructure symmetrize(const BasicStructure &structure, const std::vector<SymOp> &enforced_group);

    /// Given an integer transformation matrix, create a superstructure whose lattice is the product of
    /// the original lattice and the transformation matrix.
    BasicStructure make_superstructure(const BasicStructure &tiling_unit, const Eigen::Matrix3i transformation_matrix);

  } // namespace xtal
} // namespace CASM

#endif
