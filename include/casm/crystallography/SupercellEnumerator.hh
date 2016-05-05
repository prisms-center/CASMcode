#ifndef SupercellEnumerator_HH
#define SupercellEnumerator_HH

#include "casm/external/Eigen/Dense"

#include "casm/symmetry/SymGroup.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"

namespace CASM {
  /**
   * Given the dimensions of a square matrix and its determinant,
   * HermiteCounter will cycle through every possible matrix that
   * maintains it's Hermite normal form:
   *  -Upper triangular matrix
   *  -Determinant remains constant
   *  -row values to the right of the diagonal will always be smaller than
   *  the value of the diagonal
   *
   * For a determinant det, the initial value of the counter will be
   * a n x n identity matrix H with H(0,0)=det.
   * The final position will be a n x n identity matrix with H(n-1,n-1)=det.
   *
   * There are two main steps in the counter:
   *  -Incrementing the diagonal of the matrix such that its product remains
   *  equal to the determinant
   *  -Incrementing the upper triangular values such that they never exceed
   *  the diagonal
   *
   * The diagonal increments are achieved by working only with two adjacent
   * values at a time, distributing factors to the next diagonal element
   * only if it equals 1. If not, the next adjacent pair is selected.
   */

  class HermiteCounter {
  public:
    /// \brief constructor given the desired determinant and square matrix dimensions
    HermiteCounter(int init_determinant, int init_dim);

    //You probably will never need these. They're just here for testing more than anything.
    int pos() const;
    int att() const;
    Eigen::VectorXi diagonal() const;

  private:

    /// \brief Keeps track of the current adjacent diagonal element pair
    int m_pos;

    /// \brief The last value attempted to factorize a diagonal element
    int m_att;

    /// \brief Vector holding diagonal element values
    Eigen::VectorXi m_diagonal;

  };

  //******************************************************************************************************************//

  template <typename UnitType> class SupercellEnumerator;

  /// \brief Iterators used with SupercellEnumerator
  ///
  /// Allows iterating over unique supercells of a UnitType object (may be Lattice, eventually BasicStructure, or Structure)
  ///
  template <typename UnitType>
  class SupercellIterator {

  public:

    /// fixes alignment of m_deformation
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef std::forward_iterator_tag iterator_category;
    typedef int difference_type;
    typedef UnitType value_type;
    typedef const UnitType &reference;
    typedef const UnitType *pointer;

    SupercellIterator<UnitType>() {}

    SupercellIterator<UnitType>(const SupercellEnumerator<UnitType> &enumerator,
                                int volume);

    SupercellIterator<UnitType>(const SupercellIterator<UnitType> &B);

    SupercellIterator<UnitType> &operator=(const SupercellIterator<UnitType> &B);

    /// \brief Iterator comparison
    bool operator==(const SupercellIterator<UnitType> &B) const;

    /// \brief Iterator comparison
    bool operator!=(const SupercellIterator<UnitType> &B) const;

    /// \brief Access the supercell
    reference operator*() const;

    /// \brief Access the supercell
    pointer operator->() const;

    /// \brief const reference to the supercell matrix
    const Eigen::Matrix3i &matrix() const;

    /// \brief const reference to the SupercellEnumerator this is iterating with
    const SupercellEnumerator<UnitType> &enumerator() const;

    /// \brief Prefix increment operator. Increment to next unique supercell.
    SupercellIterator<UnitType> &operator++();

    /// \brief Postfix increment operator. Increment to next unique supercell.
    SupercellIterator<UnitType> operator++(int);


  private:

    /// \brief Uses _try_increment until the next unique supercell is found.
    void _increment();

    /// \brief Check if the current supercell matrix hermite normal form is in a canonical form
    bool _is_canonical() const;

    /// \brief Increment the supercell matrix by one (maintaining hermite normal form)
    void _try_increment();

    /// \brief Update m_super when required
    void _update_super();


    /// \brief Indicates if m_super reflects the current m_current supercell matrix
    mutable bool m_super_updated;

    /// \brief Pointer to SupercellEnumerator which holds the unit cell and point group
    const SupercellEnumerator<UnitType> *m_enum;

    /// \brief m_current supercell volume
    int m_vol;

    /// \brief Current supercell matrix
    Eigen::Matrix3i m_current;

    /// \brief A supercell, stored here so that iterator dereferencing will be OK. Only used when requested.
    mutable UnitType m_super;
  };


  /// \brief A fake container of supercell matrices
  ///
  /// Provides iterators over symmetrically unique supercells of some object of UnitType. Could be Lattice,
  /// eventually BasicStructure or Structure
  ///
  template<typename UnitType>
  class SupercellEnumerator {

  public:

    typedef unsigned long int size_type;

    typedef SupercellIterator<UnitType> const_iterator;

    /// \brief Construct a SupercellEnumerator
    ///
    /// \returns a SupercellEnumerator
    ///
    /// \param unit The thing that is tiled to form supercells. For now Lattice.
    /// \param tol Tolerance for generating the point group
    /// \param begin_volume The beginning volume to enumerate
    /// \param end_volume The past-the-last volume to enumerate
    ///
    SupercellEnumerator(UnitType unit,
                        double tol,
                        size_type begin_volume = 1,
                        size_type end_volume = ULONG_MAX);

    /// \brief Construct a SupercellEnumerator using custom point group operations
    ///
    /// \returns a SupercellEnumerator
    ///
    /// \param unit The thing that is tiled to form supercells. For now Lattice.
    /// \param point_grp Point group operations to use for checking supercell uniqueness.
    /// \param begin_volume The beginning volume to enumerate
    /// \param end_volume The past-the-last volume to enumerate
    ///
    SupercellEnumerator(UnitType unit,
                        const SymGroup &point_grp,
                        size_type begin_volume = 1,
                        size_type end_volume = ULONG_MAX);

    /// \brief Access the unit the is being made into supercells
    const UnitType &unit() const;

    /// \brief Access the unit lattice
    const Lattice &lattice() const;

    /// \brief Access the unit point group
    const SymGroup &point_group() const;

    /// \brief Set the beginning volume
    void begin_volume(size_type _begin_volume);

    /// \brief Get the beginning volume
    size_type begin_volume() const;

    /// \brief Set the end volume
    void end_volume(size_type _end_volume);

    /// \brief Get the end volume
    size_type end_volume() const;

    /// \brief A const iterator to the beginning volume
    const_iterator begin() const;

    /// \brief A const iterator to the past-the-last volume
    const_iterator end() const;

    /// \brief A const iterator to the beginning volume
    const_iterator cbegin() const;

    /// \brief A const iterator to the past-the-last volume
    const_iterator cend() const;

    /// \brief A const iterator to a specified volume
    const_iterator citerator(size_type volume) const;

  private:

    /// \brief The unit cell of the supercells
    UnitType m_unit;

    /// \brief The lattice of the unit cell
    Lattice m_lat;

    /// \brief The point group of the unit cell
    SymGroup m_point_group;

    /// \brief The first volume supercells to be iterated over (what cbegin uses)
    int m_begin_volume;

    /// \brief The past-the-last volume supercells to be iterated over (what cend uses)
    int m_end_volume;
  };


  /// \brief Return a transformation matrix that ensures a supercell of at least
  ///        some volume
  ///
  /// \param unit The thing that is tiled to form supercells. May not be the
  ///        primitive cell.
  /// \param T The transformation matrix of the unit, relative the prim. The
  ///        volume of the unit is determined from T.determinant().
  /// \param point_grp Point group operations to use for checking supercell uniqueness.
  /// \param volume The beginning volume to enumerate
  /// \param end_volume The past-the-last volume to enumerate
  /// \param fix_shape If true, enforce that S = T*m*I, where m is a scalar and
  ///        I is the identity matrix.
  ///
  /// \returns M, a transformation matrix, S = T*M, where M is an integer matrix
  ///          and S.determinant() >= volume
  ///
  ///
  template<typename UnitType>
  Eigen::Matrix3i enforce_min_volume(
    const UnitType &unit,
    const Eigen::Matrix3i &T,
    const SymGroup &point_grp,
    Index volume,
    bool fix_shape = false);



  template<typename UnitType>
  SupercellIterator<UnitType>::SupercellIterator(const SupercellEnumerator<UnitType> &enumerator,
                                                 int volume) :
    m_super_updated(false),
    m_enum(&enumerator) {

    (volume < 1) ? m_vol = 1 : m_vol = volume;
    if(volume < 1) { //Redundant if statement?
      m_vol = 1;
    }
    m_current = Eigen::Matrix3i::Identity();
    m_current(2, 2) = m_vol;

    if(!_is_canonical()) {
      _increment();
    }
  }

  template<typename UnitType>
  SupercellIterator<UnitType>::SupercellIterator(const SupercellIterator &B) {
    *this = B;
  };

  template<typename UnitType>
  SupercellIterator<UnitType> &SupercellIterator<UnitType>::operator=(const SupercellIterator &B) {
    m_enum = B.m_enum;
    m_current = B.m_current;
    m_vol = B.m_vol;
    m_super_updated = false;
    return *this;
  }

  template<typename UnitType>
  bool SupercellIterator<UnitType>::operator==(const SupercellIterator &B) const {
    return (m_enum == B.m_enum) && (m_vol == B.m_vol) && (m_current - B.m_current).isZero();
  }

  template<typename UnitType>
  bool SupercellIterator<UnitType>::operator!=(const SupercellIterator &B) const {
    return !(*this == B);
  }

  template<typename UnitType>
  typename SupercellIterator<UnitType>::reference SupercellIterator<UnitType>::operator*() const {
    if(!m_super_updated) {
      m_super = make_supercell(m_enum->unit(), m_current);
      m_super_updated = true;
    }
    return m_super;
  }

  template<typename UnitType>
  typename SupercellIterator<UnitType>::pointer SupercellIterator<UnitType>::operator->() const {
    if(!m_super_updated) {
      m_super = make_supercell(m_enum->unit(), m_current);
      m_super_updated = true;
    }
    return &m_super;
  }

  template<typename UnitType>
  const Eigen::Matrix3i &SupercellIterator<UnitType>::matrix() const {
    return m_current;
  }

  template<typename UnitType>
  const SupercellEnumerator<UnitType> &SupercellIterator<UnitType>::enumerator() const {
    return *m_enum;
  }

  // prefix
  template<typename UnitType>
  SupercellIterator<UnitType> &SupercellIterator<UnitType>::operator++() {
    _increment();
    return *this;
  }

  // postfix
  template<typename UnitType>
  SupercellIterator<UnitType> SupercellIterator<UnitType>::operator++(int) {
    SupercellIterator result(*this);
    _increment();
    return result;
  }


  template<typename UnitType>
  void SupercellIterator<UnitType>::_increment() {
    do {
      _try_increment();
    }
    while(!_is_canonical());
    m_super_updated = false;
  }

  template<typename UnitType>
  bool SupercellIterator<UnitType>::_is_canonical() const {

    Eigen::Matrix3i H, V;

    // apply point group operations to the supercell matrix to check if it is canonical form
    // S: supercell lattice column vectors,  U: unit cell lattice column vectors,  T: supercell transformation matrix
    //
    // S = U*T  and  op*S = U*T', solve for T' = f(T)
    //
    // op*U*T = U*T'
    // U.inv*op*U*T = T'

    for(int i = 0; i < m_enum->point_group().size(); i++) {

      Eigen::Matrix3d U = m_enum->lattice().lat_column_mat();
      Eigen::Matrix3d op = m_enum->point_group()[i].matrix();

      Eigen::Matrix3i transformed = iround(U.inverse() * op * U) * m_current;

      H = hermite_normal_form(transformed).first;

      // canonical only if m_current is '>' H, so if H '>' m_current, return false
      if(H(0, 0) > m_current(0, 0)) {
        return false;
      }
      if(H(0, 0) < m_current(0, 0)) {
        continue;
      }

      if(H(1, 1) > m_current(1, 1)) {
        return false;
      }
      if(H(1, 1) < m_current(1, 1)) {
        continue;
      }

      if(H(2, 2) > m_current(2, 2)) {
        return false;
      }
      if(H(2, 2) < m_current(2, 2)) {
        continue;
      }

      if(H(1, 2) > m_current(1, 2)) {
        return false;
      }
      if(H(1, 2) < m_current(1, 2)) {
        continue;
      }

      if(H(0, 2) > m_current(0, 2)) {
        return false;
      }
      if(H(0, 2) < m_current(0, 2)) {
        continue;
      }

      if(H(0, 1) > m_current(0, 1)) {
        return false;
      }
      if(H(0, 1) < m_current(0, 1)) {
        continue;
      }

    }
    return true;
  }

  template<typename UnitType>
  void SupercellIterator<UnitType>::_try_increment() {

    // order: try to increment m_current(1,2), (0,2), (0,1), (1,1), (0,0), m_vol
    //   but ensure still in valid hermite_normal_form with correct volume

    // try to increment m_current(1,2)
    if(m_current(1, 2) < m_current(1, 1) - 1) {
      m_current(1, 2)++;
      return;
    }
    m_current(1, 2) = 0;

    // try to increment m_current(0,2)
    if(m_current(0, 2) < m_current(0, 0) - 1) {
      m_current(0, 2)++;
      return;
    }
    m_current(0, 2) = 0;

    // try to increment m_current(0,1)
    if(m_current(0, 1) < m_current(0, 0) - 1) {
      m_current(0, 1)++;
      return;
    }
    m_current(0, 1) = 0;

    // try to increment m_current(1,1)
    do {
      m_current(1, 1)++;
    }
    while((m_vol / m_current(0, 0) % m_current(1, 1) != 0) && (m_current(1, 1) <= m_vol));
    if(m_current(1, 1) <= m_vol) {
      m_current(2, 2) = m_vol / (m_current(0, 0) * m_current(1, 1));
      return;
    }
    m_current(1, 1) = 1;

    // try to increment m_current(0,0)
    do {
      m_current(0, 0)++;
    }
    while((m_vol % m_current(0, 0) != 0) && (m_current(0, 0) <= m_vol));
    if(m_current(0, 0) <= m_vol) {
      m_current(2, 2) = m_vol / (m_current(0, 0) * m_current(1, 1));
      return;
    }
    m_current(0, 0) = 1;

    // increment m_vol
    m_vol++;
    m_current(2, 2) = m_vol;

    return;

  }

  template<typename UnitType>
  const UnitType &SupercellEnumerator<UnitType>::unit() const {
    return m_unit;
  }

  template<typename UnitType>
  const Lattice &SupercellEnumerator<UnitType>::lattice() const {
    return m_lat;
  }

  template<typename UnitType>
  const SymGroup &SupercellEnumerator<UnitType>::point_group() const {
    return m_point_group;
  }

  template<typename UnitType>
  void SupercellEnumerator<UnitType>::begin_volume(size_type _begin_volume) {
    m_begin_volume = _begin_volume;
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::size_type SupercellEnumerator<UnitType>::begin_volume() const {
    return m_begin_volume;
  }

  template<typename UnitType>
  void SupercellEnumerator<UnitType>::end_volume(size_type _end_volume) {
    m_end_volume = _end_volume;
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::size_type SupercellEnumerator<UnitType>::end_volume() const {
    return m_end_volume;
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::begin() const {
    return const_iterator(*this, m_begin_volume);
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::end() const {
    return const_iterator(*this, m_end_volume);
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::cbegin() const {
    return const_iterator(*this, m_begin_volume);
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::cend() const {
    return const_iterator(*this, m_end_volume);
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::citerator(size_type volume) const {
    return SupercellIterator<UnitType>(*this, volume);
  }


  // declare specializations for Lattice

  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    double tol,
                                                    size_type begin_volume,
                                                    size_type end_volume);

  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    const SymGroup &point_grp,
                                                    size_type begin_volume,
                                                    size_type end_volume);

  template<>
  Eigen::Matrix3i enforce_min_volume<Lattice>(
    const Lattice &unit,
    const Eigen::Matrix3i &T,
    const SymGroup &point_grp,
    Index volume,
    bool fix_shape);


  /// \brief Return canonical hermite normal form of the supercell matrix, and op used to find it
  std::pair<Eigen::MatrixXi, Eigen::MatrixXd>
  canonical_hnf(const Eigen::MatrixXi &T, const BasicStructure<Site> &unitcell);

}

#endif

