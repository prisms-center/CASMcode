#ifndef ADAPTERINTERFACE_HH
#define ADAPTERINTERFACE_HH

namespace CASM {
  namespace adapter {
    /// In order to convert types between modules (e.g. turn SymOp
    /// into SymOpType so the symmetry operation is primed for the crystallography
    /// module), this Adapter functor must be defined
    template <typename FromType, typename ToType>
    struct Adapter {
    };
  } // namespace adapter
} // namespace CASM

#endif
