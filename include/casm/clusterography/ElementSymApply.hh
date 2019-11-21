#ifndef ELEMENTSYMAPPLY_HH
#define ELEMENTSYMAPPLY_HH

#include "casm/global/definitions.hh"
#include "casm/misc/CRTPBase.hh"
#include "casm/symmetry/SymTools.hh"
#include "casm/kinetics/OccupationTransformation.hh"    //TODO: Forwad declary sym::copy_apply?
#include <memory>

namespace CASM {

  class SymOp;

  namespace xtal {
    class Structure;
  }

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

    template <typename Element>
    Element copy_apply_impl(const SymOp &op, const Element &obj) const {
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

    template <typename Element>
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

    template <typename Element>
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
  class CopyApplyForEach : public Base {

  public:
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    template <typename Element>
    Element copy_apply_impl(const SymOp &op, const Element &obj) const {
      auto result = obj;
      for(auto &e : result) {
        CASM::apply(op, e);
      }
      return result;
    }
  };

  namespace sym {

    /**
     * Functor class to apply sym via CopyApplyWithPrim.
     * A shared prim is given at construction, which is then forwarded to
     * sym::copy_apply(const SymOp&, const T&, const Structure&), applying
     * this function to the element.
     */

    class CopyApplyWithPrim : public CASM::CopyApplyWithPrim<CRTPBase<CopyApplyWithPrim>> {
      typedef xtal::Structure PrimType;
      typedef std::shared_ptr<const PrimType> PrimType_ptr;
      typedef CASM::CopyApplyWithPrim<CRTPBase<CopyApplyWithPrim>> Base;
      typedef typename Base::MostDerived MostDerived;
      using Base::derived;
      friend Base;

    public:
      explicit CopyApplyWithPrim(PrimType_ptr prim_ptr) : m_prim_ptr(prim_ptr) {}
      template <typename Element>
      Element operator()(const SymOp &op, const Element &obj) const {
        return Base::copy_apply_impl(op, obj);
      }

    private:
      const PrimType &prim() const {
        return *m_prim_ptr;
      }
      PrimType_ptr m_prim_ptr;
    };

    /**
     * Functor class to apply sym via CopyApplyForEachWithPrim.
     * A shared prim is given at construction, which is then forwarded to
     * sym::copy_apply(const SymOp&, const T&, const Structure&), applying
     * this function to each element of Element (e.g. if Element is a cluster
     * of UnitCellCoord, use this functior to call copy_apply on each UnitCellCoord),
     * while retaining the interface of copy_apply(const SymOp&, const T&)
     */

    class CopyApplyElementWiseWithPrim : public CopyApplyForEachWithPrim<CRTPBase<CopyApplyElementWiseWithPrim>> {
      typedef xtal::Structure PrimType;
      typedef std::shared_ptr<const PrimType> PrimType_ptr;
      typedef CopyApplyForEachWithPrim<CRTPBase<CopyApplyElementWiseWithPrim>> Base;
      typedef typename Base::MostDerived MostDerived;
      using Base::derived;
      friend Base;

    public:
      explicit CopyApplyElementWiseWithPrim(PrimType_ptr prim_ptr) : m_prim_ptr(prim_ptr) {}
      template <typename Element>
      Element operator()(const SymOp &op, const Element &obj) const {
        return Base::copy_apply_impl(op, obj);
      }

    private:
      const PrimType &prim() const {
        return *m_prim_ptr;
      }
      PrimType_ptr m_prim_ptr;
    };

    /**
     * Functor class to apply sym via CopyApply.
     * This functor requires no arguments at construction, and relies on
     * CASM::copy_apply(const SymOp&, const T&)
     */

    class CopyApplyDefault : public CASM::CopyApply<CRTPBase<CopyApplyDefault>> {
      typedef CASM::CopyApply<CRTPBase<CopyApplyDefault>> Base;
      typedef typename Base::MostDerived MostDerived;
      using Base::derived;

    public:

      CopyApplyDefault() {}
      template<typename DiscardedType>
      CopyApplyDefault(const DiscardedType &_) {}

      template <typename Element>
      Element operator()(const SymOp &op, const Element &obj) const {
        return Base::copy_apply_impl(op, obj);
      }

    private:
    };
  }

} // namespace CASM

#endif
