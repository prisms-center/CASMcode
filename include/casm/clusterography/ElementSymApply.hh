#ifndef ELEMENTSYMAPPLY_HH
#define ELEMENTSYMAPPLY_HH

#include "casm/global/definitions.hh"
#include "casm/symmetry/SymTools.hh"

namespace CASM {

  class SymOp;

  /**
   * This is a CRTP mixin for the ClusterSymCompare classes.
   * Sometimes copy_apply must be called with a prim provided
   * (e.g. UnitCellCoord). The specialization of the ClusterSymCompare
   * classes must therefore pick this way of applying symmetry
   * to the Element if Element is of type UnitCellCoord.
   */

  template <typename Base>
  class CopyApplyWithPrim : public Base {

  public:
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    template<typename Element>
    Element copy_apply_impl(const SymOp &op, const Element &obj) const {
      return sym::copy_apply(op, obj, this->derived().prim());
    }
  };

  /**
   * This is a CRTP mixin for the ClusterSymCompare classes.
   * When applying symmetry only requires a SymOp and the
   * Element itself, this way of applying symmetry should
   * be invoked.
   */

  template <typename Base>
  class CopyApply : public Base {

  public:
    using MostDerived = typename Base::MostDerived;
    using Base::derived;

    template<typename Element>
    Element copy_apply_impl(const SymOp &op, const Element &obj) const {
      return CASM::copy_apply(op, obj);
    }
  };
} // namespace CASM

#endif
