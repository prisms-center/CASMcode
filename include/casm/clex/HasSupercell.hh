#ifndef CASM_HasSupercell
#define CASM_HasSupercell

#include <memory>

namespace CASM {

class PrimClex;
class Structure;
class Supercell;

/// \brief Implements PrimClex dependent functions
///
/// - Useful if MostDerived class has const Supercell& MostDerived::supercell()
/// const;
template <typename Base>
class HasSupercell : public Base {
 public:
  typedef typename Base::MostDerived MostDerived;
  using Base::derived;

  const Structure &prim() const;

  std::shared_ptr<Structure const> const &shared_prim() const;

  double crystallography_tol() const;

  bool has_primclex() const;

  const PrimClex &primclex() const;
};

}  // namespace CASM

#endif
