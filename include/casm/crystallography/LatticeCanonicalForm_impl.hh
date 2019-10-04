#ifndef CASM_LatticeCanonicalForm_impl
#define CASM_LatticeCanonicalForm_impl

#include "casm/crystallography/LatticeCanonicalForm.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/crystallography/Niggli.hh"

namespace CASM {
  namespace xtal {

    // --- Lattice canonical form finding ---

    template<typename Base>
    LatticeCanonicalForm<Base>::LatticeCanonicalForm() {}

    /// Uses Lattice point group
    template<typename Base>
    bool LatticeCanonicalForm<Base>::is_canonical() const {
      return is_canonical(calc_point_group(derived()));
    }

    /// Check if *this = other*U, where U is unimodular
    template<typename Base>
    bool LatticeCanonicalForm<Base>::is_equivalent(const Lattice &other) const {
      LatticeIsEquivalent is_equiv(derived());
      return is_equiv(other);
    }

    /// Check if this->canonical_form() == other.canonical_form()
    template<typename Base>
    bool LatticeCanonicalForm<Base>::is_sym_equivalent(const Lattice &other) const {
      return canonical_form() == other.canonical_form();
    }

    /// Canonical equivalent lattice, using this lattice's point group
    template<typename Base>
    Lattice LatticeCanonicalForm<Base>::canonical_form() const {
      return canonical_form(calc_point_group(derived()));
    }


    /// Uses provided group
    ///
    /// - True if lat_column_mat is approximately equal to the lat_column_mat of the canonical_form
    template<typename Base>
    bool LatticeCanonicalForm<Base>::is_canonical(std::vector<SymOp> const &g) const {
      return almost_equal(derived().lat_column_mat(), canonical_form(g).lat_column_mat(), _tol());
    }

    /// Uses provided group to find 'to_canonical' SymOp
    ///
    /// - Returns first SymOp for which canonical_form.is_equivalent(apply(op, *this))
    /// - Note that that copy_apply(this->to_canonical(), *this).is_canonical()
    ///   may be false because they may be equivalent, but without identical
    ///   lat_column_mat().
    template<typename Base>
    SymOp LatticeCanonicalForm<Base>::to_canonical(std::vector<SymOp> const &g) const {
      return _canonical_equivalent_lattice(derived(), g, _tol()).second;
    }

    /// Returns the inverse of to_canonical(g)
    /// - Note that that
    ///     copy_apply(this->from_canonical(), this->canonical_form()) == *this
    ///   may be false because they may be equivalent, but without identical
    ///   lat_column_mat().
    template<typename Base>
    SymOp LatticeCanonicalForm<Base>::from_canonical(std::vector<SymOp> const &g) const {
      return to_canonical(g).inverse();
    }

    /// Canonical equivalent lattice, using the provided group
    template<typename Base>
    Lattice LatticeCanonicalForm<Base>::canonical_form(std::vector<SymOp> const &g) const {
      return canonical_equivalent_lattice(derived(), g, _tol());
    }


    /// \brief Construct the subgroup that leaves a lattice unchanged
    template<typename Base>
    SymGroup LatticeCanonicalForm<Base>::invariant_subgroup(std::vector<SymOp> const &super_grp) const {
      return invariant_subgroup(super_grp.begin(), super_grp.end());
    }

    /// \brief Construct the subgroup for which this->is_equivalent(copy_apply(op, *this))
    template<typename Base>
    template<typename SymOpIt>
    SymGroup LatticeCanonicalForm<Base>::invariant_subgroup(SymOpIt begin, SymOpIt end) const {
      SymGroup result;
      invariant_subgroup(begin, end, std::back_inserter(result));
      result.set_lattice(derived());
      return result;
    }

    /// \brief Construct the subgroup for which this->is_equivalent(copy_apply(op, *this))
    template<typename Base>
    template<typename SymOpIt, typename OutputIt>
    OutputIt LatticeCanonicalForm<Base>::invariant_subgroup(SymOpIt begin, SymOpIt end, OutputIt result) const {
      LatticeIsEquivalent is_equiv(derived());
      return std::copy_if(begin, end, result, is_equiv);
    }


  }
}

#endif
