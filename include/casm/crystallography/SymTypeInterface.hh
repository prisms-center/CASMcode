#ifndef SYMTYPEINTERFACE_HH
#define SYMTYPEINTERFACE_HH

#include "casm/crystallography/SymType.hh"

// Specialize these functions for an automatic Adapter that
// can convert between your custom symmetry object, and
// xtal::SymOpType. The conversion allows you to pass your
// custom symmetry objects to the crystallography module,
// which only deals with xtal::SymOpType.

namespace CASM {
  namespace xtal {
    /// Accessor for SymOpType. Returns transformation matrix (Cartesian).
    template <typename ExternSymOpType>
    SymOpMatrixType &get_matrix(ExternSymOpType &op);
    template <typename ExternSymOpType>
    const SymOpMatrixType &get_matrix(const ExternSymOpType &op);

    /// Accessor for SymOpType. Returns translation vector (tau).
    template <typename ExternSymOpType>
    SymOpTranslationType &get_translation(ExternSymOpType &op);
    template <typename ExternSymOpType>
    const SymOpTranslationType &get_translation(const ExternSymOpType &op);

    /// Accessor for SymOpType. Returns whether the symmetry operation is time reversal active.
    template <typename ExternSymOpType>
    SymOpTimeReversalType get_time_reversal(const ExternSymOpType &op);
  } // namespace xtal
} // namespace CASM

#endif
