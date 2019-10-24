#ifndef CASM_LatticeCanonicalForm
#define CASM_LatticeCanonicalForm

#include "casm/CASM_global_definitions.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include <vector>

namespace CASM {
  class SymGroup;
  class SymOp;
  namespace xtal {

    class Lattice;

    namespace canonical {
      /// True if lat_column_mat is approximately equal to the lat_column_mat of the canonical_form, using
      /// the lattice point group to find the most canonical form
      bool check(const Lattice &lat);

      /// True if lat_column_mat is approximately equal to the lat_column_mat of the canonical_form, using
      /// the provided symmetry operations to find the most canonical form
      bool check(const Lattice &lat, std::vector<SymOp> const &g);


      /// Canonical equivalent lattice, using this lattice's point group
      Lattice equivalent(const Lattice &lat);

      /// Canonical equivalent lattice, using the provided group
      Lattice equivalent(const Lattice &lat, std::vector<SymOp> const &g);


      /// Uses provided group to find 'to_canonical' SymOp
      ///
      /// - Returns first SymOp for which canonical_form.is_equivalent(apply(op, *this))
      /// - Note that that copy_apply(this->to_canonical(), *this).is_canonical()
      ///   may be false because they may be equivalent, but without identical
      ///   lat_column_mat().
      Index operation_index(const Lattice &lat, std::vector<SymOp> const &g);
    }


    /// Check if canonical_form(ref_lattice) == canonical_form(other)
    /* bool lattices_are_symmetrically_equivalent(const Lattice& ref_lattice, const Lattice& other); */


    /// \brief Construct indices of the subgroup that leaves a lattice unchanged
    std::vector<Index> invariant_subgroup_indices(const Lattice &lat, std::vector<SymOp> const &super_grp);

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

  } // namespace xtal
} // namespace CASM

#endif
