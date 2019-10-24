#ifndef CASM_LatticeCanonicalForm_impl
#define CASM_LatticeCanonicalForm_impl

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/LatticeCanonicalForm.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {
  namespace xtal {
    namespace canonical {
      bool check(const Lattice &lat) {
        return canonical::check(lat, calc_point_group(lat));
      }

      bool check(const Lattice &lat, std::vector<SymOp> const &g) {
        return almost_equal(lat.lat_column_mat(), canonical::equivalent(lat, g).lat_column_mat(), lat.tol());
      }

      /* bool is_canonical(const Lattice& lat) { return is_canonical(lat, calc_point_group(lat)); } */

      Lattice equivalent(const Lattice &lat) {
        return canonical::equivalent(lat, calc_point_group(lat));
      }

      /// Uses provided group to find 'to_canonical' SymOp
      ///
      /// - Returns first SymOp for which canonical_form.is_equivalent(apply(op, *this))
      /// - Note that that copy_apply(this->to_canonical(), *this).is_canonical()
      ///   may be false because they may be equivalent, but without identical
      ///   lat_column_mat().
      Index operation_index(const Lattice &lat, std::vector<SymOp> const &g) {
        return canonical::operation_index(lat, g, lat.tol());
      }

      /// Canonical equivalent lattice, using the provided group
      Lattice equivalent_lattice(const Lattice &lat, std::vector<SymOp> const &g) {
        return canonical::equivalent(lat, g, lat.tol());
      }

    } // namespace canonical

    // --- Lattice canonical form finding ---

    /* bool lattices_are_symmetrically_equivalent(const Lattice& ref_lattice, const Lattice &other)  { */
    /*   return canonical_form(ref_lattice) == canonical_form(other); */
    /* } */

    /* bool is_canonical(const Lattice& lat, std::vector<SymOp> const &g) { */
    /*   return almost_equal(lat.lat_column_mat(), canonical_form(lat,g).lat_column_mat(), lat.tol()); */
    /* } */

    /// \brief Construct the subgroup that leaves a lattice unchanged
    std::vector<Index> invariant_subgroup_indices(const Lattice &lat, std::vector<SymOp> const &super_grp) {
      std::vector<Index> result;
      invariant_subgroup_indices(lat, super_grp, std::back_inserter(result));
      return result;
    }

    /// \brief Construct the subgroup for which this->is_equivalent(copy_apply(op, *this))
    /* template<typename Base> */
    /* std::vector<Index> LatticeCanonicalForm<Base>::invariant_subgroup_indices(std::vector<SymOp>::const_iterator begin,
     * std::vector<SymOp>::const_iterator end) const { */
    /*   std::vector<Index> result; */
    /*   invariant_subgroup_indices(begin, end, std::back_inserter(result)); */
    /*   return result; */
    /* } */

  } // namespace xtal
} // namespace CASM

#endif
