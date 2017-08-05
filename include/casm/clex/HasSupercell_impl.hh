#ifndef CASM_HasSupercell_impl
#define CASM_HasSupercell_impl

#include "casm/clex/HasSupercell.hh"

namespace CASM {

  template<typename Base>
  const PrimClex &HasSupercell<Base>::primclex() const {
    return derived().supercell().primclex();
  }

}

#endif
