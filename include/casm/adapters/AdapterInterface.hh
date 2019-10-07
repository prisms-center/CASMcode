#ifndef ADAPTERINTERFACE_HH
#define ADAPTERINTERFACE_HH

namespace CASM {
  namespace adapter {
    /// In order to convert types between modules (e.g. turn SymOp
    /// into SymOpType so the symmetry operation is primed for the crystallography
    /// module), this Adapter functor must be defined
    template <typename ToType, typename FromType>
    struct Adapter {
    };

    //TODO: Do we want this?
    /// Allows for conversion between types without having to create the Adapter
    /// functor. Instead, you can invoke the conversion with
    /// adapter::Adapt::cast<ToType,FromType>::cast(my_from_type)
    /*
    template <typename ToType, typename FromType>
    struct Adapt
    {
        template <typename... Ts>
        static ToType cast(Ts... args)
        {
            return Adapter<ToType,FromType>()(args...);
        }
    };
    */
  } // namespace adapter
} // namespace CASM

#endif
