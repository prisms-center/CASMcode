#ifndef CASM_CoordCluster
#define CASM_CoordCluster

#include <vector>

#include "casm/clusterography/GenericCluster.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {

  /* -- CoordCluster Declarations ------------------------------------- */

  class Structure;
  class SymOp;

  /// \brief A cluster of Coordinate-like elements
  ///
  /// Beyond what a GenericCluster does, a CoordCluster:
  /// - has a ptr to a primitive Structure
  /// - may be translated by UnitCell
  ///
  /// \ingroup Clusterography
  ///
  template<typename CoordType>
  class CoordCluster : public GenericCluster<CoordType> {

  public:

    typedef unsigned int size_type;
    typedef BasicStructure<Site> PrimType;

    /// \brief Construct an empty UnitCellCoordCluster
    explicit CoordCluster(const PrimType &_prim) :
      GenericCluster<CoordType>(),
      m_prim_ptr(&_prim) {}

    /// \brief Construct a CoordCluster with a range of CoordType
    template<typename InputIterator>
    CoordCluster(const PrimType &_prim,
                 InputIterator _begin,
                 InputIterator _end) :
      GenericCluster<CoordType>(_begin, _end),
      m_prim_ptr(&_prim) {}

    /// \brief Default move constructor
    CoordCluster(CoordCluster &&other) = default;

    /// \brief Default move assignment
    CoordCluster &operator=(CoordCluster && other) = default;


    /// \brief Return a reference to the primitive Structure
    const Structure &prim() const {
      return *m_prim_ptr;
    }

    /// \brief Return the coordinate corresponding to element(i)
    Coordinate coordinate(size_type i) const {
      return prim().get_site(element(i));
    }

    /// \brief Translate the cluster by a UnitCell translation
    GenericCluster<CoordType> &operator+=(UnitCell trans) {
      for(auto it = begin(); it != end(); ++it) {
        *it += trans;
      }
      return *this;
    }

    /// \brief Translate the UnitCellCoordCluster by a UnitCell translation
    GenericCluster<CoordType> &operator-=(UnitCell trans) {
      for(auto it = begin(); it != end(); ++it) {
        *it -= trans;
      }
      return *this;
    }


  private:

    const PrimType *m_prim_ptr;

  };

  /// \brief Translate a cluster
  template<typename CoordType>
  CoordCluster<CoordType> operator+(CoordCluster<CoordType> cluster, UnitCell trans);

  /// \brief Translate a cluster
  ///
  /// \relates GenericCluster
  ///
  template<typename CoordType>
  CoordCluster<CoordType> operator-(CoordCluster<CoordType> cluster, UnitCell trans);



  /* -- GenericCluster Definitions ------------------------------------- */

  /// \brief Translate a cluster
  ///
  /// \relates GenericCluster
  ///
  template<typename CoordType>
  CoordCluster<CoordType> operator+(CoordCluster<CoordType> cluster, UnitCell trans) {
    return cluster += trans;
  }

  /// \brief Translate a cluster
  ///
  /// \relates GenericCluster
  ///
  template<typename CoordType>
  CoordCluster<CoordType> operator-(CoordCluster<CoordType> cluster, UnitCell trans) {
    return cluster -= trans;
  }

}

#endif