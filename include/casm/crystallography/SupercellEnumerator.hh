#ifndef SupercellEnumerator_HH
#define SupercellEnumerator_HH

#include "casm/external/Eigen/Dense"

#include "casm/misc/cloneable_ptr.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/casm_io/Log.hh"

#include "casm/container/Counter.hh"

namespace CASM {

  class PrimClex;

  /** \defgroup LatticeEnum Lattice Enumerators
   *
   *  \ingroup Enumerator
   *  \ingroup Lattice
   *
   *  \brief Enumerates Lattice
   *
   *  @{
   */

  /// Data structure for holding supercell enumeration properties
  class ScelEnumProps {

  public:

    typedef long size_type;

    /// \brief Constructor
    ///
    /// \param begin_volume The beginning volume to enumerate
    /// \param end_volume The past-the-last volume to enumerate
    /// \param dirs String indicating which lattice vectors to enumerate
    ///        over. Some combination of 'a', 'b', and 'c', where 'a' indicates
    ///        the first lattice vector of the unit cell, 'b' the second, and 'c'
    ///        the third.
    /// \param generating_matrix This matrix, G, transforms the primitive lattice
    ///         vectors into the unit cell lattice vectors which are used to generate
    ///        supercells. So the generated supercells, S = P*G*T, where S and P
    ///        are column vector matrices of the supercell and primitive cell,
    ///        respectively, and G and T are integer tranformation matrices.
    ScelEnumProps(size_type begin_volume,
                  size_type end_volume,
                  std::string dirs = "abc",
                  Eigen::Matrix3i generating_matrix = Eigen::Matrix3i::Identity()) :
      m_begin_volume(begin_volume),
      m_end_volume(end_volume),
      m_dims(dirs.size()),
      m_dirs(dirs) {

      for(int i = 0; i < m_dirs.size(); i++) {
        if(m_dirs[i] != 'a' && m_dirs[i] != 'b' && m_dirs[i] != 'c') {
          std::string msg = "Error constructing ScelEnumProps: an element of dirs != 'a', 'b', or 'c'";
          default_err_log().error("Constructing ScelEnumProps");
          default_err_log() << msg << "\n" << std::endl;
          throw std::invalid_argument(msg);
        }
      }

      // add missing directions to 'dirs',
      // then permute 'G' so that the specified 'dirs' are first
      while(m_dirs.size() != 3) {
        if(std::find(m_dirs.begin(), m_dirs.end(), 'a') == m_dirs.end()) {
          m_dirs.push_back('a');
        }
        if(std::find(m_dirs.begin(), m_dirs.end(), 'b') == m_dirs.end()) {
          m_dirs.push_back('b');
        }
        if(std::find(m_dirs.begin(), m_dirs.end(), 'c') == m_dirs.end()) {
          m_dirs.push_back('c');
        }
      }
      Eigen::Matrix3i P = Eigen::Matrix3i::Zero();
      P(m_dirs[0] - 'a', 0) = 1;
      P(m_dirs[1] - 'a', 1) = 1;
      P(m_dirs[2] - 'a', 2) = 1;

      m_gen_mat = generating_matrix * P;
    }


    size_type begin_volume() const {
      return m_begin_volume;
    }

    size_type end_volume() const {
      return m_end_volume;
    }

    int dims() const {
      return m_dims;
    }

    std::string dirs() const {
      return m_dirs;
    }

    Eigen::Matrix3i generating_matrix() const {
      return m_gen_mat;
    }

  private:

    size_type m_begin_volume;
    size_type m_end_volume;
    int m_dims;
    std::string m_dirs;
    Eigen::Matrix3i m_gen_mat;

  };

  /// \brief Read unit cell transformation matrix from JSON input
  Eigen::Matrix3i make_unit_cell(PrimClex &primclex, const jsonParser &json);

  /// \brief Make a ScelEnumProps object from JSON input
  ScelEnumProps make_scel_enum_props(PrimClex &primclex, const jsonParser &input);

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

    /// \brief constructor to satisfy iterator requirements. Do not recommend.
    //HermiteCounter() {};

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

    /// \brief Expand a n x n Hermite normal matrix (H) into a m x m one through a m x m generating matrix (G) (e.g. for arbitrary 2D supercells)
    Eigen::MatrixXi _expand_dims(const Eigen::MatrixXi &H, const Eigen::MatrixXi &G);

    /// \brief Unroll a Hermit normal form square matrix into a vector such that it's canonical form is easy to compare
    Eigen::VectorXi _canonical_unroll(const Eigen::MatrixXi &hermit_mat);

    /// \brief Compare two integer matrices and see which one is lexicographically greatest. Returns true if H0<H1
    bool _canonical_compare(const Eigen::MatrixXi &H0, const Eigen::MatrixXi &H1);
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

    //required for all iterators
    //SupercellIterator<UnitType>(const SupercellIterator<UnitType> &B);

    //required for all iterators?
    SupercellIterator<UnitType> &operator=(const SupercellIterator<UnitType> &B);

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
    notstd::cloneable_ptr<HermiteCounter> m_current;

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
    /// \param enum_props Data structure specifying how to enumerate supercells
    /// \param tol Tolerance for generating the point group
    ///
    SupercellEnumerator(UnitType unit,
                        const ScelEnumProps &enum_props,
                        double tol);

    /// \brief Construct a SupercellEnumerator using custom point group operations
    ///
    /// \returns a SupercellEnumerator
    ///
    /// \param unit The thing that is tiled to form supercells. For now Lattice.
    /// \param point_grp Point group operations to use for checking supercell uniqueness.
    /// \param enum_props Data structure specifying how to enumerate supercells
    ///
    SupercellEnumerator(UnitType unit,
                        const SymGroup &point_grp,
                        const ScelEnumProps &enum_props);


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
    const Eigen::Matrix3i &gen_mat() const;

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

    /// \brief This matrix (G) specifies new lattice vectors to enumerate over column-wise, such that the resulting transformation matrix
    /// is G*B, where B is the block matrix constructed out of the HermiteCounter matrix H.
    const Eigen::Matrix3i m_gen_mat;

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
    m_super_updated(false),
    m_enum(&enumerator),
    m_current(notstd::make_cloneable<HermiteCounter>(volume, dims)) {
    if(enumerator.begin_volume() > enumerator.end_volume()) {
      throw std::runtime_error("The beginning volume of the SupercellEnumerator cannot be greater than the end volume!");
    }

    if(dims < 1) {
      throw std::runtime_error("Dimensions to count over must be greater than 0!");
    }
  }

  /*
  template<typename UnitType>
  SupercellIterator<UnitType>::SupercellIterator(const SupercellIterator &B)
  {
      *this = B;
  };
  */

  template<typename UnitType>
  SupercellIterator<UnitType> &SupercellIterator<UnitType>::operator=(const SupercellIterator &B) {
    m_enum = B.m_enum;
    m_current = B.m_current;
    m_super_updated = false;

    m_canon_hist = B.m_canon_hist;
    return *this;
  }

  template<typename UnitType>
  bool SupercellIterator<UnitType>::operator==(const SupercellIterator &B) const {
    return (m_enum == B.m_enum) && (matrix() - B.matrix()).isZero();
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
    return m_current->determinant();
  }

  template<typename UnitType>
  Eigen::Matrix3i SupercellIterator<UnitType>::matrix() const {
    Eigen::Matrix3i expanded = HermiteCounter_impl::_expand_dims((*m_current)(), m_enum->gen_mat());
    return canonical_hnf(expanded, m_enum->point_group(), m_enum->lattice());
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
    HermiteCounter::value_type last_determinant = m_current->determinant();
    ++(*m_current);

    if(last_determinant != m_current->determinant()) {
      m_canon_hist.clear();
    }

    while(std::find(m_canon_hist.begin(), m_canon_hist.end(), matrix()) != m_canon_hist.end()) {
      ++(*m_current);
    }

    m_super_updated = false;
  }

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
  const Eigen::Matrix3i &SupercellEnumerator<UnitType>::gen_mat() const {
    return m_gen_mat;
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
                                                    const ScelEnumProps &enum_props,
                                                    double tol);

  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    const SymGroup &point_grp,
                                                    const ScelEnumProps &enum_props);

  template<>
  Eigen::Matrix3i enforce_min_volume<Lattice>(
    const Lattice &unit,
    const Eigen::Matrix3i &T,
    const SymGroup &point_grp,
    Index volume,
    bool fix_shape);


  /// \brief Return canonical hermite normal form of the supercell matrix, and op used to find it
  Eigen::Matrix3i canonical_hnf(const Eigen::Matrix3i &T, const SymGroup &effective_pg, const Lattice &ref_lattice);

  /** @} */
}

#endif
