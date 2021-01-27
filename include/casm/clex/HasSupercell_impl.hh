#ifndef CASM_HasSupercell_impl
#define CASM_HasSupercell_impl

#include "casm/clex/HasSupercell.hh"

namespace CASM {

template <typename Base>
const Structure &HasSupercell<Base>::prim() const {
  return derived().supercell().prim();
}

template <typename Base>
std::shared_ptr<Structure const> const &HasSupercell<Base>::shared_prim()
    const {
  return derived().supercell().shared_prim();
}

template <typename Base>
double HasSupercell<Base>::crystallography_tol() const {
  return derived().supercell().crystallography_tol();
}

template <typename Base>
bool HasSupercell<Base>::has_primclex() const {
  return derived().supercell().has_primclex();
}

template <typename Base>
const PrimClex &HasSupercell<Base>::primclex() const {
  return derived().supercell().primclex();
}

}  // namespace CASM

#endif
