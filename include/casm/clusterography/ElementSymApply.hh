#ifndef ELEMENTSYMAPPLY_HH
#define ELEMENTSYMAPPLY_HH

#include "casm/global/definitions.hh"
#include "casm/misc/CRTPBase.hh"
//TODO: Fix this include insanity. This is only here because it has the
//sym::copy_apply routines that are needed for instantiation
#include "casm/symmetry/SymTools.hh"
#include "casm/clusterography/CoordCluster.hh"
#include "casm/kinetics/OccPerturbation.hh"

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
   * When applying symmetry only requires a SymOp and the
   * Element itself, this way of applying symmetry should
   * be invoked.
   */

  template <typename Base>
  class CopyApplyDefault : public Base {

  public:
    using MostDerived = typename Base::MostDerived;
    using Base::derived;

    template <typename Element>
    Element copy_apply_impl(const SymOp &op, const Element &obj) const {
      return CASM::copy_apply(op, obj);
    }
  };

  namespace sym {

    /**
     * Functor class to apply sym via CopyApplyWithPrim.
     * A shared prim is given at construction, which is then forwarded to
     * sym::copy_apply(const SymOp&, const T&, const Structure&), applying
     * this function to the element.
     */

    class CopyApplyWithPrim_f : public CASM::CopyApplyWithPrim<CRTPBase<CopyApplyWithPrim_f>> {
      typedef xtal::Structure PrimType;
      typedef std::shared_ptr<const PrimType> PrimType_ptr;
      typedef CASM::CopyApplyWithPrim<CRTPBase<CopyApplyWithPrim_f>> Base;
      typedef typename Base::MostDerived MostDerived;
      using Base::derived;
      friend Base;

    public:
      explicit CopyApplyWithPrim_f(PrimType_ptr prim_ptr) : m_prim_ptr(prim_ptr) {}
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

    class CopyApplyDefault_f : public CASM::CopyApplyDefault<CRTPBase<CopyApplyDefault_f>> {
      typedef CASM::CopyApplyDefault<CRTPBase<CopyApplyDefault_f>> Base;
      typedef typename Base::MostDerived MostDerived;
      using Base::derived;

    public:

      CopyApplyDefault_f() {}
      template<typename DiscardedType>
      CopyApplyDefault_f(const DiscardedType &_) {}

      template <typename Element>
      Element operator()(const SymOp &op, const Element &obj) const {
        return Base::copy_apply_impl(op, obj);
      }

    private:
    };
  }

} // namespace CASM

#endif
