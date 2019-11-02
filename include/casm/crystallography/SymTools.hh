#ifndef XTALSYMTOOLS_HH
#define XTALSYMTOOLS_HH

#include <vector>;
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/CASM_global_definitions.hh"

namespace CASM {
  namespace xtal {
    class Lattice;

    /// \brief Construct indices of the subgroup that leaves a lattice unchanged
    std::vector<Index> invariant_subgroup_indices(const Lattice &lat, SymOpVector const &super_grp);

    template <typename ExternSymOpVector>
    std::vector<Index> invariant_subgroup_indices(const Lattice &lat, ExternSymOpVector const &super_grp) {
      return invariant_subgroup_indices(lat, adapter::Adapter<SymOpVector, ExternSymOpVector>()(super_grp));
    }

    /// \brief Construct indices of the subgroup for which this->is_equivalent(copy_apply(op, *this))
    template <typename OutputIt>
    OutputIt invariant_subgroup_indices(const Lattice &lat, const std::vector<SymOp> &super_group, OutputIt result) {
      LatticeIsEquivalent is_equiv(lat);

      Index ix = 0;
      for(auto it = super_group.begin(); it != super_group.end(); ++it) {
        if(is_equiv(*it)) {
          *result = ix;
          ++result;
        }
        ++ix;
      }
      return result;
    }

    /// Convert a Cartesian symmetry operation representation to fractional
    Eigen::Matrix3i symmetry_matrix_to_frac(const Lattice &lat, const Eigen::Matrix3d &cart_matrix);
    /// Convert a fractional symmetry operation representation to Cartesian
    Eigen::Matrix3d symmetry_matrix_to_cart(const Lattice &lat, const Eigen::Matrix3i &cart_matrix);

    /// \brief Return a copy of the given lattice, which obeys the symmetry of the given group \param enforced_group
    Lattice symmetrize(const Lattice &lat, const std::vector<SymOp> &enforced_group);

    //TODO
    //Why does this routine take a tolerance, if the lattice itself has a tolerance?
    //Should we get rid of the tolerance inside of lattice?
    /// \brief Return a copy of the given lattice, which obeys the symmetry of its point group, when generated
    /// within the tolerance \param point_group_tolerance
    Lattice symmetrize(const Lattice &lat, double point_group_tolerance);

    /// Relax the vectors of the given lattice such that it obeys the symmetry of the given group,
    /// where the symmetry operations are given in fractional representations
    /* Lattice symmetrized_with_fractional(const Lattice& lat, const std::vector<Eigen::Matrix3i> &fractional_point_group); */

    /// \brief Populate \param point_group with the point group of this lattice
    /// \param point_group should be empty
    /// \param pg_tol can be increased to find point group of lattice vectors
    /// that are slightly distorted due to numerical noise
    std::vector<SymOp> make_point_group(Lattice const &_lat);

    /// \brief Populate \param point_group with the point group of this lattice
    /// \param point_group should be empty
    /// \param pg_tol can be increased to find point group of lattice vectors
    /// that are slightly distorted due to numerical noise
    std::vector<SymOp> make_point_group(Lattice const &_lat, double _tol);

  } // namespace xtal
} // namespace CASM

#endif
