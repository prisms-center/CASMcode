#ifndef ADAPTERINTERFACE_HH
#define ADAPTERINTERFACE_HH

namespace CASM {
  namespace adapter {
    /// In order to convert types between modules (e.g. turn SymOp
    /// into SymOpType so the symmetry operation is primed for the crystallography
    /// module), this cast function must be defined.
    template <typename FromType, typename ToType>
    struct Adapter {
      typedef FromType FromTypeIt;
      ToType operator()(const FromType &adaptable);
      ToType operator()(FromTypeIt begin, FromTypeIt end);
    };
  } // namespace adapter
} // namespace CASM

#endif
