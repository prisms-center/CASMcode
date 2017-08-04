#ifndef CASM_ClusterSymCompare
#define CASM_ClusterSymCompare

#include "casm/symmetry/SymCompare.hh"

namespace CASM {

  template<typename Derived>
  class ClusterSymCompare;
  class PrimGrid;
  class PrimClex;
  class Supercell;

  /// \brief Traits class for ClusterSymCompare
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename Derived>
  struct traits<ClusterSymCompare<Derived>> : public traits<Derived> {};

  /// \brief CRTP Base class for Cluster comparisons
  ///
  /// Implements:
  /// - 'invariants_compare_impl' using 'compare'
  /// - 'compare_impl'
  ///
  /// Does not implement:
  /// - 'prepare_impl'
  ///
  /// traits<Element> requires:
  /// - static UnitCellCoord position(const Element& el) const;
  /// - typedef <ElementInvariants> InvariantsType;
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare
  ///   - ClusterSymCompare
  ///     - IntegralClusterSymCompare (implements 'compare_impl')
  ///       - LocalSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - PrimPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - ScelPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///
  /// \ingroup Clusterography
  ///
  template<typename Derived>
  class ClusterSymCompare : public SymCompare<ClusterSymCompare<Derived> > {

  public:

    /// Element refers to Cluster, not element of Cluster
    typedef typename traits<Derived>::MostDerived MostDerived;
    typedef typename traits<Derived>::Element Element;
    typedef Element ClusterType;
    typedef typename traits<Element>::InvariantsType InvariantsType;

    /// \brief Return tolerance
    double tol() const;

  protected:

    friend class SymCompare<ClusterSymCompare<Derived> >;
    using SymCompare<ClusterSymCompare<Derived> >::m_integral_tau;

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    ClusterSymCompare(double tol);

    // compare_impl : implement in Derived

    /// \brief Orders 'prepared' elements in the same orbit
    bool invariants_compare_impl(const InvariantsType &A, const InvariantsType &B) const;

    /// \brief Compares 'prepared' elements
    bool compare_impl(const Element &A, const Element &B) const;

  protected:

    /// \brief type-specific way to get position of element
    ///
    /// - Returns traits<Element>::position(el)
    static UnitCellCoord position(const Element &el);

  private:

    double m_tol;

  };

  template <
    typename ClusterType,
    typename = typename std::enable_if <
      std::is_same <
        UnitCellCoord,
        decltype(traits<ClusterType>::position(std::declval<ClusterType>()))
        >::value >::type >
  using enable_if_integral_position = ClusterType;


  /* -- LocalSymCompare<IntegralCluster> Declaration ------------------------------------- */

  /// \brief Comparisons of GenericCluster-derived types using aperiodic symmetry
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare<Derived>
  ///   - ClusterSymCompare<Derived> (implements 'compare_impl', 'invariants_compare_impl')
  ///     - AperiodicSymCompare<ClusterType> (implements 'prepare_impl')
  ///     - PrimPeriodicSymCompare<ClusterType> (implements 'prepare_impl')
  ///     - ScelPeriodicSymCompare<ClusterType> (implements 'prepare_impl')
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename Element>
  class AperiodicSymCompare<Element/*, enable_if_integral_position<Element>*/> :
    public ClusterSymCompare<AperiodicSymCompare<Element> > {

  public:

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    AperiodicSymCompare(double tol);

  protected:

    friend class SymCompare<ClusterSymCompare<AperiodicSymCompare<Element> > >;
    using ClusterSymCompare<AperiodicSymCompare<Element> >::m_integral_tau;

    /// \brief Prepare an element for comparison
    ///
    /// - Returns sorted
    Element prepare_impl(Element obj) const;

  };


  /// \brief Comparisons of GenericCluster-derived types using prim periodic symmetry
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare<Derived>
  ///   - ClusterSymCompare<Derived> (implements 'compare_impl', 'invariants_compare_impl')
  ///     - LocalSymCompare<ClusterType> (implements 'prepare_impl')
  ///     - PrimPeriodicSymCompare<ClusterType> (implements 'prepare_impl')
  ///     - ScelPeriodicSymCompare<ClusterType> (implements 'prepare_impl')
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename Element>
  class PrimPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/> :
    public ClusterSymCompare<PrimPeriodicSymCompare<Element> > {

  public:

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    PrimPeriodicSymCompare(double tol);

    PrimPeriodicSymCompare(const PrimClex &primclex);

  protected:

    friend class SymCompare<ClusterSymCompare<PrimPeriodicSymCompare<Element> > >;
    using ClusterSymCompare<PrimPeriodicSymCompare<Element> >::m_integral_tau;
    using ClusterSymCompare<PrimPeriodicSymCompare<Element> >::position;

    /// \brief Prepare an element for comparison
    ///
    /// - Sorts and translates so that obj[0] is in the origin unit cell
    Element prepare_impl(Element obj) const;

  };

  /// \brief Comparisons of GenericCluster-derived types using supercell periodic symmetry
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare<Derived>
  ///   - ClusterSymCompare<Derived> (implements 'compare_impl', 'invariants_compare_impl')
  ///     - LocalSymCompare<ClusterType> (implements 'prepare_impl')
  ///     - PrimPeriodicSymCompare<ClusterType> (implements 'prepare_impl')
  ///     - ScelPeriodicSymCompare<ClusterType> (implements 'prepare_impl')
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename Element>
  class ScelPeriodicSymCompare<Element/*, enable_if_integral_position<Element>*/> :
    public ClusterSymCompare<ScelPeriodicSymCompare<Element> > {

  public:

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    ScelPeriodicSymCompare(const PrimGrid &prim_grid, double tol);

    /// \brief Constructor
    ScelPeriodicSymCompare(const Supercell &scel);

  protected:

    friend class SymCompare<ClusterSymCompare<ScelPeriodicSymCompare<Element> > >;
    using ClusterSymCompare<ScelPeriodicSymCompare<Element> >::m_integral_tau;
    using ClusterSymCompare<ScelPeriodicSymCompare<Element> >::position;

    /// \brief Prepare an element for comparison
    ///
    /// - Sorts UnitCellCoord and translates so that obj[0] is within the supercell
    Element prepare_impl(Element obj) const;

    const PrimGrid *m_prim_grid;
  };

}

#endif
