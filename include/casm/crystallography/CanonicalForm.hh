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
      bool check(const Lattice &lat, std::vector<CASM::SymOp> const &g);


      /**
       * Return a canonical Lattice by first converting the
       * given Lattice into the standard Niggli form,
       * followed by applying the point group of the Lattice
       * so that the one oriented in the most standard manner
       * is selected.
       */
      Lattice equivalent(Lattice const &in_lat, std::vector<CASM::SymOp> const &point_grp, double compare_tol);

      /// Canonical equivalent lattice, using this lattice's point group
      Lattice equivalent(const Lattice &lat);

      /// Canonical equivalent lattice, using the provided group
      Lattice equivalent(const Lattice &lat, std::vector<CASM::SymOp> const &g);


      /// Uses provided group to find 'to_canonical' SymOp
      ///
      /// - Returns first SymOp for which canonical_form.is_equivalent(apply(op, *this))
      /// - Note that that copy_apply(this->to_canonical(), *this).is_canonical()
      ///   may be false because they may be equivalent, but without identical
      ///   lat_column_mat().
      Index operation_index(const Lattice &lat, std::vector<CASM::SymOp> const &g);

      /// Return the index of the operation that makes the lattice canonical
      Index operation_index(Lattice const &in_lat, std::vector<CASM::SymOp> const &point_grp, double compare_tol);
    }

  } // namespace xtal
} // namespace CASM

#endif
