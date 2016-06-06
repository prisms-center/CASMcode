#ifndef CASM_CoordCluster
#define CASM_CoordCluster

#include <vector>

#include "casm/clusterography/GenericCluster.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {

  /* -- CoordCluster Declarations ------------------------------------- */

  class Structure;
  class SymOp;

  template<typename CoordType>
  class CoordCluster;

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /** \defgroup CoordCluster

      \brief Functions and classes related to CoordCluster
      \ingroup Clusterography
  */

  namespace CASM_TMP {

    /// \brief Traits class for CoordCluster<CoordType>
    ///
    /// \ingroup CoordCluster
    ///
    template<typename CoordType>
    struct traits<CoordCluster<CoordType> > {
      typedef CoordType Element;
      typedef ClusterInvariants<CoordCluster<CoordType> > InvariantsType;
    };
  }

  /// \brief A cluster of Coordinate-like elements
  ///
  /// Beyond what a GenericCluster does, a CoordCluster:
  /// - has a ptr to a primitive Structure
  /// - may be translated by UnitCell
  /// - expects static_cast<Coordinate>(CoordType) to be valid
  ///
  /// \ingroup CoordCluster
  ///
  template<typename CoordType>
  class CoordCluster : public ElementWiseSymCluster<CoordCluster<CoordType> > {

  public:

    typedef unsigned int size_type;
    typedef BasicStructure<Site> PrimType;

    /// \brief Construct an empty UnitCellCoordCluster
    explicit CoordCluster(const PrimType &_prim) :
      ElementWiseSymCluster<CoordCluster>(),
      m_prim_ptr(&_prim) {}

    /// \brief Construct a CoordCluster with a range of CoordType
    template<typename InputIterator>
    CoordCluster(const PrimType &_prim,
                 InputIterator _begin,
                 InputIterator _end) :
      ElementWiseSymCluster<CoordCluster>(_begin, _end),
      m_prim_ptr(&_prim) {}

    /// \brief Default copy constructor
    CoordCluster(const CoordCluster &other) = default;

    /// \brief Default assignment constructor
    CoordCluster &operator=(const CoordCluster &other) = default;

    /// \brief Default move constructor
    CoordCluster(CoordCluster &&other) = default;

    /// \brief Default move assignment
    CoordCluster &operator=(CoordCluster &&other) = default;


    /// \brief Return a reference to the primitive Structure
    const PrimType &prim() const {
      return *m_prim_ptr;
    }

    /// \brief Return the coordinate corresponding to element(i)
    Coordinate coordinate(size_type i) const {
      return static_cast<Coordinate>(this->element(i));
    }

    /// \brief Return the min pair distance
    double min_length() const {
      return this->invariants().displacement().front();
    }

    /// \brief Return the max pair distance
    double max_length() const {
      return this->invariants().displacement().back();
    }

    /// \brief Translate the cluster by a UnitCell translation
    CoordCluster &operator+=(UnitCell trans) {
      for(auto it = this->begin(); it != this->end(); ++it) {
        *it += trans;
      }
      return *this;
    }

    /// \brief Translate the UnitCellCoordCluster by a UnitCell translation
    CoordCluster &operator-=(UnitCell trans) {
      for(auto it = this->begin(); it != this->end(); ++it) {
        *it -= trans;
      }
      return *this;
    }


  private:

    const PrimType *m_prim_ptr;

  };

  /// \brief Translate a cluster
  ///
  /// \ingroup CoordCluster
  ///
  template<typename CoordType>
  CoordCluster<CoordType> operator+(CoordCluster<CoordType> cluster, UnitCell trans);

  /// \brief Translate a cluster
  ///
  /// \ingroup CoordCluster
  ///
  template<typename CoordType>
  CoordCluster<CoordType> operator-(CoordCluster<CoordType> cluster, UnitCell trans);



  /* -- GenericCluster Definitions ------------------------------------- */

  /// \brief Translate a cluster
  ///
  /// \ingroup CoordCluster
  ///
  template<typename CoordType>
  CoordCluster<CoordType> operator+(CoordCluster<CoordType> cluster, UnitCell trans) {
    return cluster += trans;
  }

  /// \brief Translate a cluster
  ///
  /// \ingroup CoordCluster
  ///
  template<typename CoordType>
  CoordCluster<CoordType> operator-(CoordCluster<CoordType> cluster, UnitCell trans) {
    return cluster -= trans;
  }

}

#endif