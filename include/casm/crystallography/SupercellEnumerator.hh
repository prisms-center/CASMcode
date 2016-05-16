#ifndef SupercellEnumerator_HH
#define SupercellEnumerator_HH

#include "casm/external/Eigen/Dense"

#include "casm/symmetry/SymGroup.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"

#include "casm/container/Counter.hh"

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
   * In addition, this class is limited to SQUARE matrices, and NONZERO
   * determinants. The idea is to use it for supcercell enumerations,
   * where these conditions are always met.
   *
   * For a determinant det, the initial value of the counter will be
   * a n x n identity matrix H with H(0,0)=det.
   * The final position will be a n x n identity matrix with H(n-1,n-1)=det.
   * Once the final position for a particular determinant is reached,
   * the counter starts over with the next integer determinant value.
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
    typedef Eigen::VectorXi::Scalar value_type;
    typedef CASM::Index Index;


    /// \brief constructor given the desired determinant and square matrix dimensions
    HermiteCounter(int init_determinant, int init_dim);

    //You probably will never need these. They're just here for testing more than anything.
    //Either way, they're safe to call.
    Index position() const;
    Eigen::VectorXi diagonal() const;
    //value_type low() const;
    //value_type high() const;

    /// \brief Get the current matrix the counter is on
    Eigen::MatrixXi current() const;

    /// \brief Get the current determinant
    value_type determinant() const;

    /// \brief Get the dimensions of *this
    Index dim() const;

    /// \brief reset the counter to the first iteration of the current determinant
    void reset_current();

    /// \brief Skip the remaining iterations and start at the next determinant value
    void next_determinant();

    /// \brief Reset the diagonal to the specified determinant and set the other elements to zero
    void jump_to_determinant(value_type new_det);

    /// \brief Jump to the next available HNF matrix.
    HermiteCounter &operator++();

    /// \brief Get the current matrix the counter is on
    Eigen::MatrixXi operator()() const;


  private:

    /// \brief Keeps track of the current diagonal element that needs to be factored
    Index m_pos;

    /// \brief Vector holding diagonal element values
    Eigen::VectorXi m_diagonal;

    /// \brief unrolled vector of the upper triangle (does not include diagonal elements)
    EigenVectorXiCounter m_upper_tri;

    /// \brief Go to the next values of diagonal elements that keep the same determinant
    Index _increment_diagonal();

  };

  namespace HermiteCounter_impl {
    /// \brief Find the next factor of the specified position and share with next element. Use attempt as starting point.
    HermiteCounter::Index _spill_factor(Eigen::VectorXi &diag, HermiteCounter::Index position, HermiteCounter::value_type attempt);

    /// \brief Spill the next factor of the specified element with its neighbor, and return new position
    HermiteCounter::Index next_spill_position(Eigen::VectorXi &diag, HermiteCounter::Index position);

    /// \brief Determine the number of elements in the upper triangular matrix (excluding diagonal)
    HermiteCounter::Index upper_size(HermiteCounter::Index init_dim);

    /// \brief Create a counter for the elements above the diagonal based on the current diagonal value
    EigenVectorXiCounter _upper_tri_counter(const Eigen::VectorXi &current_diag);

    /// \brief Assemble a matrix diagonal and unrolled upper triangle values into a matrix
    Eigen::MatrixXi _zip_matrix(const Eigen::VectorXi &current_diag, const Eigen::VectorXi &current_upper_tri);

    /// \brief Expand a n x n Hermite normal matrix into a m x m one (e.g. for 2D supercells)
    Eigen::MatrixXi _expand_dims_old(const Eigen::MatrixXi &hermit_mat, const Eigen::VectorXi &active_dims);

    /// \brief Expand a n x n Hermite normal matrix (H) into a m x m one through a m x m transformation matrix (T) (e.g. for arbitrary 2D supercells)
    Eigen::MatrixXi _expand_dims(const Eigen::MatrixXi &H, const Eigen::MatrixXi &T);

    /// \brief Unroll a Hermit normal form square matrix into a vector such that it's canonical form is easy to compare
    Eigen::VectorXi _canonical_unroll(const Eigen::MatrixXi &hermit_mat);
  }

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
                                int volume,
                                int dims);

    //SupercellIterator<UnitType>(const SupercellIterator<UnitType> &B); TODO

    //SupercellIterator<UnitType> &operator=(const SupercellIterator<UnitType> &B); TODO

    /// \brief Iterator comparison
    bool operator==(const SupercellIterator<UnitType> &B) const;

    /// \brief Iterator comparison
    bool operator!=(const SupercellIterator<UnitType> &B) const;

    /// \brief Access the supercell
    reference operator*() const;

    /// \brief Access the supercell
    pointer operator->() const;

    /// \brief constructed supercell matrix
    Eigen::Matrix3i matrix() const;

    /// \brief current volume
    HermiteCounter::value_type volume() const;

    /// \brief const reference to the SupercellEnumerator this is iterating with
    const SupercellEnumerator<UnitType> &enumerator() const;

    /// \brief Prefix increment operator. Increment to next unique supercell.
    SupercellIterator<UnitType> &operator++();

    /// \brief Postfix increment operator. Increment to next unique supercell.
    //SupercellIterator<UnitType> operator++(int);  TODO


  private:

    /// \brief Uses _try_increment until the next unique supercell is found.
    void _increment();

    /// \brief Check if the current supercell matrix hermite normal form is in a canonical form
    //bool _is_canonical() const;

    /// \brief Increment the supercell matrix by one (maintaining hermite normal form)
    void _try_increment();

    /// \brief Update m_super when required
    void _update_super();


    /// \brief Indicates if m_super reflects the current m_current supercell matrix
    mutable bool m_super_updated;

    /// \brief Pointer to SupercellEnumerator which holds the unit cell and point group
    const SupercellEnumerator<UnitType> *m_enum;

    /// \brief Current supercell matrix in HermitCounter form
    HermiteCounter m_current;

    /// \brief A supercell, stored here so that iterator dereferencing will be OK. Only used when requested.
    mutable UnitType m_super;

    /// \brief Keep track of the HNF matrices for the current determinant value
    std::vector<Eigen::Matrix3i> m_canon_hist;

  };


  /// \brief A fake container of supercell matrices
  ///
  /// Provides iterators over symmetrically unique supercells of some object of UnitType. Could be Lattice,
  /// eventually BasicStructure or Structure
  ///
  template<typename UnitType>
  class SupercellEnumerator {

  public:

    typedef long int size_type;

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
                        size_type begin_volume,
                        size_type end_volume,
                        int init_dims = 3,
                        Eigen::Matrix3i init_trans_mat = Eigen::Matrix3i::Identity());

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
                        size_type begin_volume,
                        size_type end_volume,
                        int init_dims = 3,
                        Eigen::Matrix3i init_trans_mat = Eigen::Matrix3i::Identity());

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

    /// \brief Get the transformation matrix that's being applied to the unit vectors
    const Eigen::Matrix3i &trans_mat() const;

    /// \brief Get the dimensions of the enumerator (1D, 2D or 3D)
    int dimension() const;

    /// \brief A const iterator to the beginning volume, specify here how the iterator should jump through the enumeration
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
    const UnitType m_unit;

    /// \brief The lattice of the unit cell
    Lattice m_lat;

    /// \brief The point group of the unit cell
    SymGroup m_point_group;         //factor group...?

    /// \brief The first volume supercells to be iterated over (what cbegin uses)
    const int m_begin_volume;

    /// \brief The past-the-last volume supercells to be iterated over (what cend uses)
    const int m_end_volume;

    /// \brief The transformation matrix to go from 1x1 or 2x2 Hermite normal form to 3x3 (for guiding enumeration). Vectors specified column-wise.
    const Eigen::Matrix3i m_trans_mat;

    /// \brief The number of lattice directions the enumeration is being done in
    const int m_dims;

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
                                                 int volume,
                                                 int dims):
    m_current(volume, dims),
    m_super_updated(false),
    m_enum(&enumerator) {
    if(enumerator.begin_volume() > enumerator.end_volume()) {
      throw std::runtime_error("The beginning volume of the SupercellEnumerator cannot be greater than the end volume!");
    }

    /*
    if(!_is_canonical())
    {
        _increment();
    }
    */
  }

  /*
  template<typename UnitType>
  SupercellIterator<UnitType>::SupercellIterator(const SupercellIterator &B)
  {
      *this = B;
  };

  template<typename UnitType>
  SupercellIterator<UnitType> &SupercellIterator<UnitType>::operator=(const SupercellIterator &B)
  {
      m_enum = B.m_enum;
      m_current = B.m_current;
      //m_vol = B.m_vol;
      m_super_updated = false;
      return *this;
  }
  */

  template<typename UnitType>
  bool SupercellIterator<UnitType>::operator==(const SupercellIterator &B) const {
    return (m_enum == B.m_enum) && (m_current() - B.m_current()).isZero();
  }

  template<typename UnitType>
  bool SupercellIterator<UnitType>::operator!=(const SupercellIterator &B) const {
    return !(*this == B);
  }

  template<typename UnitType>
  typename SupercellIterator<UnitType>::reference SupercellIterator<UnitType>::operator*() const {
    if(!m_super_updated) {
      m_super = make_supercell(m_enum->unit(), matrix());
      m_super_updated = true;
    }
    return m_super;
  }

  template<typename UnitType>
  typename SupercellIterator<UnitType>::pointer SupercellIterator<UnitType>::operator->() const {
    if(!m_super_updated) {
      m_super = make_supercell(m_enum->unit(), matrix());
      m_super_updated = true;
    }
    return &m_super;
  }

  template<typename UnitType>
  HermiteCounter::value_type SupercellIterator<UnitType>::volume() const {
    return m_current.determinant();
  }

  template<typename UnitType>
  Eigen::Matrix3i SupercellIterator<UnitType>::matrix() const {
    Eigen::Matrix3i expanded = HermiteCounter_impl::_expand_dims(m_current(), m_enum->trans_mat());
    return canonical_hnf(expanded, m_enum->point_group(), m_enum->lattice()).first;
    //return hermite_normal_form(expanded).first;
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

  /*
  // postfix
  template<typename UnitType>
  SupercellIterator<UnitType> SupercellIterator<UnitType>::operator++(int)
  {
      SupercellIterator result(*this);
      _increment();
      return result;
  }
  */


  template<typename UnitType>
  void SupercellIterator<UnitType>::_increment() {
    m_canon_hist.push_back(matrix());
    HermiteCounter::value_type last_determinant = m_current.determinant();
    ++m_current;

    if(last_determinant != m_current.determinant()) {
      m_canon_hist.clear();
    }

    while(std::find(m_canon_hist.begin(), m_canon_hist.end(), matrix()) != m_canon_hist.end()) {
      ++m_current;
    }

    std::cout << matrix() << std::endl << std::endl;
    m_super_updated = false;
  }

  /*
  template<typename UnitType>
  bool SupercellIterator<UnitType>::_is_canonical() const
  {
      // apply point group operations to the supercell matrix to check if it is canonical form
      // S: supercell lattice column vectors,  U: unit cell lattice column vectors,  T: supercell transformation matrix
      //
      // S = U*T  and  op*S = U*T', solve for T' = f(T)
      //
      // op*U*T = U*T'
      // U.inv*op*U*T = T'

      const Eigen::Matrix3i curr_mat = matrix();
      const Eigen::VectorXi unrolled_current = HermiteCounter_impl::_canonical_unroll(curr_mat);

      for(Index i = 0; i < m_enum->point_group().size(); i++)
      {
          const Eigen::Matrix3d U = m_enum->lattice().lat_column_mat();
          const Eigen::Matrix3d op = m_enum->point_group()[i].matrix();
          const Eigen::Matrix3i transformed = iround(U.inverse() * op * U) * curr_mat;
          //const Eigen::Matrix3i transformed = iround(U.inverse() * op * U) * m_current;
          const Eigen::Matrix3i H = hermite_normal_form(transformed).first;
          const Eigen::VectorXi unrolled_H = HermiteCounter_impl::_canonical_unroll(H);

          // canonical only if m_current is '>' H, so if H '>' m_current, return false
          for(Index i = 0; i < unrolled_current.size(); i++)
          {
              if(unrolled_H(i) > unrolled_current(i))
              {
                  return false;
              }
              else if(unrolled_H(i) < unrolled_current(i))
              {
                  break;
              }
          }
      }
      return true;
  }
  */

  //********************************************************************************************************//

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
  const Eigen::Matrix3i &SupercellEnumerator<UnitType>::trans_mat() const {
    return m_trans_mat;
  }

  template<typename UnitType>
  int SupercellEnumerator<UnitType>::dimension() const {
    return m_dims;
  }

  /*
  template<typename UnitType>
  void SupercellEnumerator<UnitType>::begin_volume(size_type _begin_volume)
  {
      m_begin_volume = _begin_volume;
  }
  */

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::size_type SupercellEnumerator<UnitType>::begin_volume() const {
    return m_begin_volume;
  }

  /*
  template<typename UnitType>
  void SupercellEnumerator<UnitType>::end_volume(size_type _end_volume)
  {
      m_end_volume = _end_volume;
  }
  */

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::size_type SupercellEnumerator<UnitType>::end_volume() const {
    return m_end_volume;
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::begin() const {
    return const_iterator(*this, m_begin_volume, dimension());
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::end() const {
    return const_iterator(*this, m_end_volume, dimension());
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::cbegin() const {
    return const_iterator(*this, m_begin_volume, dimension());
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::cend() const {
    return const_iterator(*this, m_end_volume, dimension());
  }

  template<typename UnitType>
  typename SupercellEnumerator<UnitType>::const_iterator SupercellEnumerator<UnitType>::citerator(size_type volume) const {
    return SupercellIterator<UnitType>(*this, volume, dimension());
  }


  // declare specializations for Lattice

  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    double tol,
                                                    size_type begin_volume,
                                                    size_type end_volume,
                                                    int init_dims,
                                                    Eigen::Matrix3i init_trans_mat);

  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    const SymGroup &point_grp,
                                                    size_type begin_volume,
                                                    size_type end_volume,
                                                    int init_dims,
                                                    Eigen::Matrix3i init_trans_mat);

  template<>
  Eigen::Matrix3i enforce_min_volume<Lattice>(
    const Lattice &unit,
    const Eigen::Matrix3i &T,
    const SymGroup &point_grp,
    Index volume,
    bool fix_shape);


  /// \brief Return canonical hermite normal form of the supercell matrix, and op used to find it
  std::pair<Eigen::MatrixXi, Eigen::MatrixXd>
  canonical_hnf(const Eigen::MatrixXi &T, const SymGroup &effective_pg, const Lattice &ref_lattice);

}

#endif

