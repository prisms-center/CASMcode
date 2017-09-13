#ifndef CASM_CoordCluster
#define CASM_CoordCluster

#include <vector>

#include "casm/symmetry/PermuteIterator.hh"
#include "casm/clusterography/GenericCluster.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/clex/HasCanonicalForm.hh"

namespace CASM {

  /* -- CoordCluster Declarations ------------------------------------- */

  class Structure;
  class SymOp;

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /** \defgroup CoordCluster

      \brief Functions and classes related to CoordCluster
      \ingroup Clusterography
  */

  /// \brief A cluster of Coordinate-like elements
  ///
  /// Beyond what a GenericCluster does, a CoordCluster:
  /// - has coordinates of elements
  /// - may be translated by UnitCell, by translating each element
  /// - has a list of displacements between coordinates, min_length, max_length
  ///
  /// Requires:
  ///   - Structure& MostDerived::prim_impl();
  ///   - Element& Element::operator+=(UnitCell trans);
  ///   - Coordinate MostDerived::coordinate_impl(size_type i) const;
  ///   - const std::vector<double>& MostDerived::InvariantsType::displacement() const;
  ///   - _Base inherits from GenericCluster
  ///
  /// \ingroup CoordCluster
  ///
  template<typename _Base>
  class GenericCoordCluster : public Translatable<GenericCluster<_Base>> {

  public:

    typedef Translatable<GenericCluster<_Base>> Base;
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    typedef typename traits<MostDerived>::Element Element;
    typedef typename traits<MostDerived>::InvariantsType InvariantsType;
    typedef typename traits<MostDerived>::size_type size_type;


    /// \brief Return the coordinate corresponding to element(i)
    Coordinate coordinate(size_type i) const {
      return derived().coordinate_impl(i);
    }

    const std::vector<double> &displacement() const {
      return this->invariants().displacement();
    }

    /// \brief Return the min pair distance, or 0.0 if size() <= 1
    double min_length() const {
      if(this->size() <= 1) {
        return 0.0;
      }
      return displacement().front();
    }

    /// \brief Return the max pair distance, or 0.0 if size() <= 1
    double max_length() const {
      if(this->size() <= 1) {
        return 0.0;
      }
      return displacement().back();
    }

    /// \brief Translate the cluster by a UnitCell translation
    MostDerived &operator+=(UnitCell trans) {
      for(auto it = this->begin(); it != this->end(); ++it) {
        *it += trans;
      }
      return this->derived();
    }

  };

  template<typename CoordType>
  struct traits<CoordCluster<CoordType>> {
    typedef CoordType Element;
    typedef ClusterInvariants<CoordCluster<CoordType>> InvariantsType;
    static CoordType position(const CoordCluster<CoordType> &clust);
    typedef unsigned int size_type;
  };

  template<typename CoordType>
  class CoordCluster : public
    CanonicalForm<ElementWiseSymApply<GenericCoordCluster<CRTPBase<CoordCluster<CoordType>>>>> {

  public:

    typedef Structure PrimType;
    typedef typename traits<CoordCluster<CoordType>>::Element Element;
    typedef typename traits<CoordCluster<CoordType>>::InvariantsType InvariantsType;
    typedef typename traits<CoordCluster<CoordType>>::size_type size_type;

    CoordCluster(const PrimType &prim) :
      m_prim_ptr(&prim) {}

    template<typename Iterator>
    CoordCluster(const PrimType &prim, Iterator begin, Iterator end) :
      m_prim_ptr(&prim),
      m_element(begin, end) {}

    const PrimType &prim() const {
      return *m_prim_ptr;
    }

    /// \brief Access vector of elements
    std::vector<Element> &elements() {
      this->reset_invariants();
      return m_element;
    }

    /// \brief const Access vector of elements
    const std::vector<Element> &elements() const {
      return m_element;
    }

  protected:

    friend GenericCoordCluster<CRTPBase<CoordCluster<CoordType>>>;

    Coordinate coordinate_impl(size_type i) const {
      return static_cast<Coordinate>(this->element(i));
    }

  private:

    std::vector<Element> m_element;
    const PrimType *m_prim_ptr;

  };

  template<typename CoordType>
  CoordType traits<CoordCluster<CoordType>>::position(const CoordCluster<CoordType> &clust) {
    return clust[0];
  }

}

#endif
