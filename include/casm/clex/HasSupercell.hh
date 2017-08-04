#ifndef CASM_HasSupercell
#define CASM_HasSupercell

#include "casm/clex/HasPrimClex.hh"

namespace CASM {

  /// \brief Implements PrimClex dependent functions
  ///
  /// - Useful if MostDerived class has const Supercell& MostDerived::supercell() const;
  template<typename _Base>
  class HasSupercell : public HasPrimClex<_Base> {
  public:

    typedef HasPrimClex<_Base> Base;
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    const PrimClex &primclex() const;

  };

}

#endif
