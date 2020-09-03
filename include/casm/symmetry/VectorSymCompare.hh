#ifndef CASM_VectorSymCompare
#define CASM_VectorSymCompare

#include "casm/symmetry/SymCompare.hh"
#include "casm/symmetry/SymGroupRep.hh"
namespace CASM {

  template<typename Base>
  class VectorSymCompare;

  class VectorInvariants {
  public:
    VectorInvariants(Eigen::VectorXcd const &vector);

    double cols() const;
    double norm() const;

  private:
    double m_cols;
    double m_norm;
  };

  /// \brief Check if VectorInvariants are equal
  bool almost_equal(VectorInvariants const &A_invariants, VectorInvariants const &B_invariants, double tol);

  /// \brief Compare ClusterInvariants
  bool compare(VectorInvariants const &A_invariants, VectorInvariants const &B_invariants, double tol);

  template<typename Element, typename SymApply>
  class DirectionSymCompare;

  /// \brief Traits class for DirectionSymCompare
  template<typename _Element, typename _SymApply>
  struct traits<DirectionSymCompare<_Element, _SymApply>> {
    typedef _Element Element;
    typedef _SymApply SymApply;
    typedef DirectionSymCompare<_Element, _SymApply> MostDerived;
    typedef VectorInvariants InvariantsType;
  };

  template<typename Element, typename SymApply>
  class SubspaceSymCompare;

  /// \brief Traits class for SubspaceSymCompare
  template<typename _Element, typename _SymApply>
  struct traits<SubspaceSymCompare<_Element, _SymApply>> {
    typedef _Element Element;
    typedef _SymApply SymApply;
    typedef SubspaceSymCompare<_Element, _SymApply> MostDerived;
    typedef VectorInvariants InvariantsType;
  };

  /// \brief Traits class for VectorSymCompare
  ///
  template<typename Base>
  struct traits<VectorSymCompare<Base>> : public traits<Base> {};

  /// \brief CRTP Base class for Vector comparisons
  ///
  /// Implements:
  /// - 'invariants_compare_impl' using 'compare'
  /// - 'compare_impl'
  /// - 'canonical_transform_impl'
  ///
  /// Does not implement
  /// - 'representation_prepare_impl'
  /// - 'copy_apply_impl'
  ///
  /// traits<Element> requires:
  /// - typedef <ElementInvariants> InvariantsType;
  ///
  /// The VectorSymCompare hierarchy:
  /// - SymCompare
  ///   - EigenSymCompare
  ///     - DirectionSymCompare
  ///     - SubspaceSymCompare
  ///
  ///
  template<typename _Base>
  class EigenSymCompare : public _Base {

  public:

    /// Element refers to Vector, not element of Vector
    /*
    typedef typename traits<Derived>::MostDerived MostDerived;
    typedef typename traits<Derived>::Element Element;
    typedef Element VectorType;
    typedef typename traits<Element>::InvariantsType InvariantsType;
    */

    using Base = _Base;
    using MostDerived = typename Base::MostDerived;
    using Element = typename traits<MostDerived>::Element;
    using SymApply = typename traits<MostDerived>::SymApply;
    using typename Base::InvariantsType;
    using Base::derived;

    /// \brief Return tolerance
    double tol() const;

  protected:

    friend _Base;

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    template<typename ...Args>
    EigenSymCompare(double tol, Args &&... args) :
      m_tol(tol),
      m_apply(std::forward<Args>(args)...) {}


    /// \brief Make orbit invariants from one element in the orbit
    InvariantsType make_invariants_impl(Element const &obj) const;

    /// \brief Orders 'prepared' elements in the same orbit
    bool invariants_compare_impl(
      InvariantsType const &A_invariants, InvariantsType const &B_invariants) const;

    /// \brief Compares 'prepared' elements
    bool compare_impl(const Element &A, const Element &B) const;

    /// \brief Returns transformation that takes 'obj' to its prepared (canonical) form
    // For now, this is the the sorting permutation
    std::unique_ptr<SymOpRepresentation> canonical_transform_impl(Element const &obj)const;

    Element copy_apply_impl(SymOp const &op, Element obj) const {
      return m_apply(op, obj);
    }

    /// \brief Spatial prepare does not apply -- element is returned unchanged
    Element spatial_prepare_impl(Element obj) const;

    /// \brief Spatial prepare does not apply -- transform is always identity
    SymOp const &spatial_transform_impl() const;

  private:

    double m_tol;

    SymApply m_apply;

    /// Not applicable to VectorSymCompare -- always identity
    /// - Default SymOp constructor initializes to identity
    SymOp m_spatial_transform;
  };

  template<typename Element>
  class EigenSymRepApply {
  public:
    EigenSymRepApply(SymGroupRep const &_rep) :
      m_rep(_rep) {}

    Element operator()(SymOp const &op, Element obj) const {
      return (*m_rep.MatrixXd(op.index())) * obj;
    }
  private:
    SymGroupRep m_rep;
  };

  template<typename Element, typename SymApply>
  class DirectionSymCompare : public EigenSymCompare<SymCompare<CRTPBase<DirectionSymCompare<Element, SymApply>>>> {
  public:
    using Base = EigenSymCompare<SymCompare<CRTPBase<DirectionSymCompare<Element, SymApply>>>>;

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    template<typename ...Args>
    DirectionSymCompare(Args &&... args) :
      Base(std::forward<Args>(args)...) {}

  protected:
    friend class EigenSymCompare<SymCompare<CRTPBase<DirectionSymCompare<Element, SymApply>>>>;
    friend class SymCompare<CRTPBase<DirectionSymCompare<Element, SymApply>>>;


    /// \brief Prepare an element for comparison (by doing nothing)
    Element representation_prepare_impl(Element obj) const {
      return obj;
    }
  };

  template<typename Element, typename SymApply>
  class SubspaceSymCompare : public EigenSymCompare<SymCompare<CRTPBase<SubspaceSymCompare<Element, SymApply>>>> {
  public:
    using Base = EigenSymCompare<SymCompare<CRTPBase<SubspaceSymCompare<Element, SymApply>>>>;

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    template<typename ...Args>
    SubspaceSymCompare(Args &&... args) :
      Base(std::forward<Args>(args)...) {}

  protected:
    friend class EigenSymCompare<SymCompare<CRTPBase<SubspaceSymCompare<Element, SymApply>>>>;
    friend class SymCompare<CRTPBase<SubspaceSymCompare<Element, SymApply>>>;

    /// \brief Prepare an element for comparison by rotating, sorting, and changing sign of column vectors
    Element representation_prepare_impl(Element obj) const;

  };


}

#endif
