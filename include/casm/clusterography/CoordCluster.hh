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
      typedef CoordCluster<CoordType> MostDerived;
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
    typedef Structure PrimType;

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
    CoordCluster &operator=(CoordCluster && other) = default;


    /// \brief Return a reference to the primitive Structure
    const PrimType &prim() const {
      return *m_prim_ptr;
    }

    /// \brief Return the coordinate corresponding to element(i)
    Coordinate coordinate(size_type i) const {
      return static_cast<Coordinate>(this->element(i));
    }

    /// \brief Return the min pair distance, or 0.0 if size() <= 1
    double min_length() const {
      if(this->size() <= 1) {
        return 0.0;
      }
      return this->invariants().displacement().front();
    }

    /// \brief Return the max pair distance, or 0.0 if size() <= 1
    double max_length() const {
      if(this->size() <= 1) {
        return 0.0;
      }
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

    CoordCluster operator+(const UnitCell &trans) const {
      CoordCluster res(*this);
      return res += trans;
    }

    CoordCluster operator-(const UnitCell &trans) const {
      CoordCluster res(*this);
      return res -= trans;
    }


  private:

    const PrimType *m_prim_ptr;

  };

}

#endif