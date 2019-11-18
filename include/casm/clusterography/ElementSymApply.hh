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
      //TODO: Is there a clever way to have the prim stored in here instead?
      return sym::copy_apply(op, obj, this->derived().prim());
    }
  };

  /**
   * This is a CRTP mixin for the ClusterSymCompare classes.
   * Sometimes copy_apply must be called with a prim provided
   * on each of the inner elements (e.g. CoordCluster<UnitCellCoord>)
   * The specialization of the ClusterSymCompare
   * classes must therefore pick this way of applying symmetry
   * to the Element if Element is of type CoordCluster<UnitCellCoord>.
   * This is done by iterating though each of the inner elements of
   * Element, and applying symmetry to them individually.
   */

  template <typename Base>
  class CopyApplyForEachWithPrim : public Base {

  public:
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    template<typename Element>
    Element copy_apply_impl(const SymOp &op, const Element &obj) const {
      auto result = obj;
      for(auto &e : result) {
        sym::apply(op, e, this->derived().prim());
      }
      return result;
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

  /**
   * This is a CRTP mixin for the ClusterSymCompare classes.
   * Sometimes copy_apply must be called
   * on each of the inner elements (e.g. OccPerturbation)
   * The specialization of the ClusterSymCompare
   * classes must therefore pick this way of applying symmetry
   * to the Element if Element is of type OccPerturbation.
   * This is done by iterating though each of the inner elements of
   * Element, and applying symmetry to them individually.
   */

  template <typename Base>
  class CopyApplyForEach: public Base {

  public:
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    template<typename Element>
    Element copy_apply_impl(const SymOp &op, const Element &obj) const {
      auto result = obj;
      for(auto &e : result) {
        CASM::apply(op, e);
      }
      return result;
    }
  };

} // namespace CASM

#endif
