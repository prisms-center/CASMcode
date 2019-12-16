#ifndef XTALSYMTOOLS_HH
#define XTALSYMTOOLS_HH

#include <vector>
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/global/definitions.hh"

namespace CASM {
  namespace sym {
    //TODO: These could be template specializations of what's in symmetry/SymTools.hh
    //but I'm not sure how we're generalizing the apply/copy_apply stuff yet

    /// \brief Apply SymOp to a Lattice
    xtal::Lattice &apply(const xtal::SymOp &op, xtal::Lattice &lat);

    /// \brief Copy and apply SymOp to a Lattice
    xtal::Lattice copy_apply(const xtal::SymOp &op, xtal::Lattice lat_copy);

    template <typename ExternSymOp>
    xtal::Lattice copy_apply(const ExternSymOp &op, const xtal::Lattice &lat) {
      return copy_apply(adapter::Adapter<xtal::SymOp, ExternSymOp>()(op), lat);
    }
  }
  namespace xtal {

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

    // TODO
    // Unimplemented. Is this the same as cart2frac and frac2cart?
    /// Convert a Cartesian symmetry operation representation to fractional
    /* Eigen::Matrix3i symmetry_matrix_to_frac(const Lattice &lat, const Eigen::Matrix3d &cart_matrix); */
    /// Convert a fractional symmetry operation representation to Cartesian
    /* Eigen::Matrix3d symmetry_matrix_to_cart(const Lattice &lat, const Eigen::Matrix3i &cart_matrix); */

    /// \brief Return a copy of the given lattice, which obeys the symmetry of the given group \param enforced_group
    Lattice symmetrize(const Lattice &lat, const std::vector<SymOp> &enforced_group);

    // TODO
    // Why does this routine take a tolerance, if the lattice itself has a tolerance?
    // Should we get rid of the tolerance inside of lattice?
    /// \brief Return a copy of the given lattice, which obeys the symmetry of its point group, when generated
    /// within the tolerance \param point_group_tolerance
    Lattice symmetrize(const Lattice &lat, double point_group_tolerance);

    /// Relax the vectors of the given lattice such that it obeys the symmetry of the given group,
    /// where the symmetry operations are given in fractional representations
    /* Lattice symmetrized_with_fractional(const Lattice& lat, const std::vector<Eigen::Matrix3i> &fractional_point_group);
     */

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

    //************************************************************************************************************************//

    /// Check if there is a symmetry operation, op, and transformation matrix T,
    ///   such that scel is a superlattice of the result of applying op to unit
    ///
    /// \returns pair corresponding to first successful op and T, or with op=end if not successful
    template<typename Object, typename OpIterator>
    std::pair<OpIterator, Eigen::Matrix3d> is_equivalent_superlattice(
      const Object &scel,
      const Object &unit,
      OpIterator begin,
      OpIterator end,
      double tol) {

      std::pair<bool, Eigen::Matrix3d> res;
      for(auto it = begin; it != end; ++it) {
        res = is_superlattice(scel, sym::copy_apply(*it, unit), tol);
        if(res.first) {
          return std::make_pair(it, res.second);
        }
      }
      return std::make_pair(end, res.second);
    }

    ///\brief returns Lattice that is smallest possible superlattice of all input Lattice
    ///
    /// SymOpIterator are provided to apply to each Lattice in an attempt
    /// to find the smallest possible superduperlattice of all symmetrically transformed Lattice
    template<typename LatIterator, typename SymOpIterator>
    Lattice make_equivalent_superduperlattice(LatIterator begin,
                                              LatIterator end,
                                              SymOpIterator op_begin,
                                              SymOpIterator op_end) {

      Lattice best = *begin;
      for(auto it = ++begin; it != end; ++it) {
        Lattice tmp_best = make_superduperlattice(best, *it);
        for(auto op_it = op_begin; op_it != op_end; ++op_it) {
          Lattice test = make_superduperlattice(best, sym::copy_apply(*op_it, *it));
          if(std::abs(volume(test)) < std::abs(volume(tmp_best))) {
            tmp_best = test;
          }
        }
        best = tmp_best;
      }
      return best;
    }

  } // namespace xtal

} // namespace CASM

namespace CASM {
  namespace xtal {
    /// Given a symmetry group, the basis of the structure will have
    /// each operation applied to it. The resulting set of basis
    /// from performing these operations will be averaged out,
    /// yielding a new average basis.
    template <typename StructureType>
    StructureType symmetrize(const StructureType &structure, const std::vector<SymOp> &enforced_group) {
      // All your sites need to be within
      auto symmetrized_structure = structure;
      symmetrized_structure.reset();

      // First make a copy of your current basis
      // This copy will eventually become the new average basis.
      auto avg_basis = structure.basis();

      // Loop through given symmetry group an fill a temporary "operated basis"
      decltype(avg_basis) operbasis;

      // Loop through given symmetry group an fill a temporary "operated basis"
      for(Index rf = 0; rf < enforced_group.size(); rf++) {
        operbasis.clear();
        for(Index b = 0; b < symmetrized_structure.basis().size(); b++) {
          operbasis.push_back(enforced_group[rf] * symmetrized_structure.basis()[b]);
        }
        // Now that you have a transformed basis, find the closest mapping of atoms
        // Then average the distance and add it to the average basis
        for(Index b = 0; b < symmetrized_structure.basis().size(); b++) {
          double smallest = 1000000;
          Coordinate bshift(symmetrized_structure.lattice()), tshift(symmetrized_structure.lattice());
          for(Index ob = 0; ob < operbasis.size(); ob++) {
            double dist = operbasis[ob].min_dist(symmetrized_structure.basis()[b], tshift);
            if(dist < smallest) {
              bshift = tshift;
              smallest = dist;
            }
          }
          bshift.cart() *= (1.0 / enforced_group.size());
          avg_basis[b] += bshift;
        }
      }
      symmetrized_structure.set_basis(avg_basis);
      symmetrized_structure.update(); // TODO: Do we want this?
      return symmetrized_structure;
    }

    template <typename StructureType, typename ExternSymOpVector>
    StructureType symmetrize(const StructureType &structure, const ExternSymOpVector &enforced_group) {
      return xtal::symmetrize(structure, adapter::Adapter<SymOpVector, ExternSymOpVector>()(enforced_group));
    }
  } // namespace xtal
} // namespace CASM

#endif
