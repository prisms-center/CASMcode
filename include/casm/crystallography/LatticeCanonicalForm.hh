#ifndef CASM_LatticeCanonicalForm
#define CASM_LatticeCanonicalForm

#include <vector>
#include "casm/CASM_global_definitions.hh"

namespace CASM {
  class SymGroup;
  class SymOp;
  namespace xtal {

    class Lattice;

    /// Lattice canonical form finding
    template<typename Base>
    class LatticeCanonicalForm : public Base {
    public:

      using Base::derived;

      LatticeCanonicalForm();

      /// Uses Lattice point group
      bool is_canonical() const;

      /// Check if *this = other*U, where U is unimodular
      bool is_equivalent(const Lattice &other) const;

      /// Check if this->canonical_form() == other.canonical_form()
      bool is_sym_equivalent(const Lattice &other) const;

      /// Canonical equivalent lattice, using this lattice's point group
      Lattice canonical_form() const;


      /// Uses provided group
      ///
      /// - True if lat_column_mat is approximately equal to the lat_column_mat of the canonical_form
      bool is_canonical(std::vector<SymOp> const &g) const;

      /// Uses provided group to find 'to_canonical' SymOp
      ///
      /// - Returns first SymOp for which canonical_form.is_equivalent(apply(op, *this))
      /// - Note that that copy_apply(this->to_canonical(), *this).is_canonical()
      ///   may be false because they may be equivalent, but without identical
      ///   lat_column_mat().
      SymOp to_canonical(std::vector<SymOp> const &g) const;

      /// Returns the inverse of to_canonical(g)
      /// - Note that that
      ///     copy_apply(this->from_canonical(), this->canonical_form()) == *this
      ///   may be false because they may be equivalent, but without identical
      ///   lat_column_mat().
      SymOp from_canonical(std::vector<SymOp> const &g) const;

      /// Canonical equivalent lattice, using the provided group
      Lattice canonical_form(std::vector<SymOp> const &g) const;


      /// \brief Construct indices of the subgroup that leaves a lattice unchanged
      std::vector<Index> invariant_subgroup_indices(std::vector<SymOp> const &super_grp) const;

      /// \brief Construct indices of the subgroup for which this->is_equivalent(copy_apply(op, *this))
      std::vector<Index> invariant_subgroup_indices(std::vector<SymOp>::const_iterator begin, std::vector<SymOp>::const_iterator end) const;

      /// \brief Construct indices of the subgroup for which this->is_equivalent(copy_apply(op, *this))
      template<typename OutputIt>
      OutputIt invariant_subgroup_indices(std::vector<SymOp>::const_iterator begin, std::vector<SymOp>::const_iterator end, OutputIt result) const;

    private:

      double _tol() const {
        return derived().tol();
      }

    };
  }
}

#endif
