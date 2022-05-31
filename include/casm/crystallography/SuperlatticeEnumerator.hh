#ifndef SuperlatticeEnumerator_HH
#define SuperlatticeEnumerator_HH

#include "casm/crystallography/HermiteCounter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/Dense"
#include "casm/misc/cloneable_ptr.hh"

// Including this file allows passing arbitrary
// symmetry types from the outside without requiring
// conversion.
#include "casm/crystallography/Adapter.hh"

namespace CASM {
namespace xtal {

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
  ///         vectors into the unit cell lattice vectors which are used to
  ///         generate
  ///        supercells. So the generated supercells, S = P*G*T, where S and P
  ///        are column vector matrices of the supercell and primitive cell,
  ///        respectively, and G and T are integer tranformation matrices.
  /// \param diagonal_only If true, restrict T to diagonal matrices
  /// \param fixed_shape If true, restrict T to diagonal matrices with
  ///     diagonal coefficients [m, 1, 1] (1d), [m, m, 1] (2d),
  ///     or [m, m, m] (3d), where the dimension == dirs.size().
  ScelEnumProps(size_type begin_volume, size_type end_volume,
                std::string dirs = "abc",
                Eigen::Matrix3i generating_matrix = Eigen::Matrix3i::Identity(),
                bool diagonal_only = false, bool fixed_shape = false)
      : m_begin_volume(begin_volume),
        m_end_volume(end_volume),
        m_dims(dirs.size()),
        m_dirs(dirs),
        m_diagonal_only(diagonal_only),
        m_fixed_shape(fixed_shape) {
    if (begin_volume < 1) {
      std::string msg = "Error constructing ScelEnumProps: begin_volume < 1";
      throw std::invalid_argument(msg);
    }

    for (int i = 0; i < m_dirs.size(); i++) {
      if (m_dirs[i] != 'a' && m_dirs[i] != 'b' && m_dirs[i] != 'c') {
        std::string msg =
            "Error constructing ScelEnumProps: an element of dirs != 'a', 'b', "
            "or 'c'";
        throw std::invalid_argument(msg);
      }
    }

    if (fixed_shape) {
      m_diagonal_only = true;
    }

    // add missing directions to 'dirs',
    // then permute 'G' so that the specified 'dirs' are first
    while (m_dirs.size() != 3) {
      if (std::find(m_dirs.begin(), m_dirs.end(), 'a') == m_dirs.end()) {
        m_dirs.push_back('a');
      }
      if (std::find(m_dirs.begin(), m_dirs.end(), 'b') == m_dirs.end()) {
        m_dirs.push_back('b');
      }
      if (std::find(m_dirs.begin(), m_dirs.end(), 'c') == m_dirs.end()) {
        m_dirs.push_back('c');
      }
    }
    Eigen::Matrix3i P = Eigen::Matrix3i::Zero();
    P(m_dirs[0] - 'a', 0) = 1;
    P(m_dirs[1] - 'a', 1) = 1;
    P(m_dirs[2] - 'a', 2) = 1;

    m_gen_mat = generating_matrix * P;
  }

  size_type begin_volume() const { return m_begin_volume; }

  size_type end_volume() const { return m_end_volume; }

  int dims() const { return m_dims; }

  std::string dirs() const { return m_dirs; }

  Eigen::Matrix3i generating_matrix() const { return m_gen_mat; }

  bool diagonal_only() const { return m_diagonal_only; }

  bool fixed_shape() const { return m_fixed_shape; }

 private:
  size_type m_begin_volume;
  size_type m_end_volume;
  int m_dims;
  std::string m_dirs;
  Eigen::Matrix3i m_gen_mat;
  bool m_diagonal_only;
  bool m_fixed_shape;
};

//******************************************************************************************************************//

class SuperlatticeEnumerator;

/// \brief Iterators used with SuperlatticeEnumerator
///
/// Allows iterating over unique superlattices
///

class SuperlatticeIterator {
 public:
  /// fixes alignment of m_deformation
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef std::forward_iterator_tag iterator_category;
  typedef int difference_type;
  typedef Lattice value_type;
  typedef const Lattice &reference;
  typedef const Lattice *pointer;

  SuperlatticeIterator() {}

  SuperlatticeIterator(const SuperlatticeEnumerator &enumerator, int volume,
                       int dims);

  // required for all iterators
  // SuperlatticeIterator<UnitType>(const SuperlatticeIterator<UnitType> &B);

  // required for all iterators?
  SuperlatticeIterator &operator=(const SuperlatticeIterator &B);

  /// \brief Iterator comparison
  bool operator==(const SuperlatticeIterator &B) const;

  /// \brief Iterator comparison
  bool operator!=(const SuperlatticeIterator &B) const;

  /// \brief Access the supercell
  reference operator*() const;

  /// \brief Access the supercell
  pointer operator->() const;

  /// \brief constructed supercell matrix
  Eigen::Matrix3i const &matrix() const;

  /// \brief current volume
  HermiteCounter::value_type volume() const;

  /// \brief const reference to the SuperlatticeEnumerator this is iterating
  /// with
  const SuperlatticeEnumerator &enumerator() const;

  /// \brief Prefix increment operator. Increment to next unique supercell.
  SuperlatticeIterator &operator++();

  /// \brief Postfix increment operator. Increment to next unique supercell.
  // SuperlatticeIterator<UnitType> operator++(int);  TODO

 private:
  /// \brief Advance m_current until the next unique, valid supercell is found.
  void _increment();

  /// \brief Check m_current for uniqueness and diagonal_only, fixed_shape
  /// validity
  bool _current_is_valid_and_unique() const;

  /// \brief Advance m_current by one, updating flags and history
  void _advance_one();

  /// \brief Advance m_current if it is invalid, updating flags and history
  void _advance_if_invalid();

  /// \brief Pointer to SuperlatticeEnumerator which holds the unit cell and
  /// point group
  const SuperlatticeEnumerator *m_enum;

  /// \brief Counter over hermite normal form matrices in the dimensions being
  /// enumerated
  notstd::cloneable_ptr<HermiteCounter> m_current;

  /// \brief Indicates if m_super reflects the current m_current matrix
  mutable bool m_super_updated;

  /// \brief A supercell, stored here so that iterator dereferencing will be OK.
  /// Only used when requested.
  mutable Lattice m_super;

  /// \brief Keep track of the HNF matrices for the current determinant value
  std::vector<Eigen::Matrix3i> m_canon_hist;

  /// \brief Indicates if m_matrix reflects the current m_current matrix
  mutable bool m_matrix_updated;

  /// \brief The transformation matrix to m_super; m_super = m_enum->unit() *
  /// m_matrix
  mutable Eigen::Matrix3i m_matrix;
};

/// \brief A fake container of supercell matrices
///
/// Provides iterators over symmetrically unique superlattices.

class SuperlatticeEnumerator {
 public:
  typedef long int size_type;

  typedef SuperlatticeIterator const_iterator;

  /// \brief Construct a SuperlatticeEnumerator using custom point group
  /// operations
  ///
  /// \returns a SuperlatticeEnumerator
  ///
  /// \param unit The thing that is tiled to form supercells. For now Lattice.
  /// \param point_grp Point group operations to use for checking supercell
  /// uniqueness. \param enum_props Data structure specifying how to enumerate
  /// supercells
  ///
  SuperlatticeEnumerator(const Lattice &unit, const SymOpVector &point_grp,
                         const ScelEnumProps &enum_props);

  /// Initialize with the unit lattice to be tiled, and specification for how to
  /// enumerate the superlattices. The symmetry operations are given with an
  /// iterator range, whose type must either have an Adapter functor defined, or
  /// the same external accessors as SymOpType
  template <typename ExternSymGroupTypeIt>
  SuperlatticeEnumerator(ExternSymGroupTypeIt begin, ExternSymGroupTypeIt end,
                         const Lattice &unit, const ScelEnumProps &enum_props)
      : SuperlatticeEnumerator(
            unit,
            adapter::Adapter<SymOpVector, ExternSymGroupTypeIt>()(begin, end),
            enum_props) {}

  /// \brief Access the unit the is being made into superlattices
  const Lattice &unit() const;

  /// \brief Access the unit point group
  const SymOpVector &point_group() const;

  /// \brief Set the beginning volume
  void begin_volume(size_type _begin_volume);

  /// \brief Get the beginning volume
  size_type begin_volume() const;

  /// \brief Set the end volume
  void end_volume(size_type _end_volume);

  /// \brief Get the end volume
  size_type end_volume() const;

  /// \brief Get the transformation matrix that's being applied to the unit
  /// vectors
  const Eigen::Matrix3i &gen_mat() const;

  /// \brief Get the dimensions of the enumerator (1D, 2D or 3D)
  int dimension() const;

  /// \brief If true, T, of S=P*G*T is restricted to diagonal matrices
  bool diagonal_only() const;

  /// \brief If true, T, of S=P*G*T is restricted to diagonal matrices with
  ///     diagonal coefficients [m, 1, 1] (1d), [m, m, 1] (2d),
  ///     or [m, m, m] (3d).
  bool fixed_shape() const;

  /// \brief A const iterator to the beginning volume, specify here how the
  /// iterator should jump through the enumeration
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
  const Lattice m_unit;

  /// \brief The lattice of the unit cell
  /* Lattice m_lat; */

  /// \brief The point group of the unit cell
  SymOpVector m_point_group;

  /// \brief The first volume supercells to be iterated over (what cbegin uses)
  const int m_begin_volume;

  /// \brief The past-the-last volume supercells to be iterated over (what cend
  /// uses)
  const int m_end_volume;

  /// \brief This matrix (G) specifies new lattice vectors to enumerate over
  /// column-wise, such that the resulting transformation matrix is G*B, where B
  /// is the block matrix constructed out of the HermiteCounter matrix H.
  const Eigen::Matrix3i m_gen_mat;

  /// \brief The number of lattice directions the enumeration is being done in
  const int m_dims;

  /// If true, restrict T to diagonal matrices, in S=P*G*T
  bool m_diagonal_only;

  /// \param fixed_shape If true, restrict T, in S=P*G*T, to diagonal matrices
  /// with
  ///     diagonal coefficients [m, 1, 1] (1d), [m, m, 1] (2d),
  ///     or [m, m, m] (3d).
  bool m_fixed_shape;
};

//********************************************************************************************************//

/// \brief Return a transformation matrix that ensures a supercell of at least
///        some volume
///
/// \param unit The thing that is tiled to form supercells. May not be the
///        primitive cell.
/// \param T The transformation matrix of the unit, relative the prim. The
///        volume of the unit is determined from T.determinant().
/// \param point_grp Point group operations to use for checking supercell
/// uniqueness. \param volume The beginning volume to enumerate \param
/// end_volume The past-the-last volume to enumerate \param fix_shape If true,
/// enforce that S = T*m*I, where m is a scalar and
///        I is the identity matrix.
///
/// \returns M, a transformation matrix, S = T*M, where M is an integer matrix
///          and S.determinant() >= volume
///
///
Eigen::Matrix3i enforce_min_volume(const Lattice &unit,
                                   const Eigen::Matrix3i &T,
                                   const SymOpVector &point_grp, Index volume,
                                   bool fix_shape = false);

template <typename ExternSymGroupTypeIt>
Eigen::Matrix3i enforce_min_volume(ExternSymGroupTypeIt begin,
                                   ExternSymGroupTypeIt end,
                                   const Lattice &unit,
                                   const Eigen::Matrix3i &T, Index volume,
                                   bool fix_shape = false) {
  return enforce_min_volume(
      unit, T,
      adapter::Adapter<SymOpVector, ExternSymGroupTypeIt>()(begin, end), volume,
      fix_shape);
}

/// \brief Return canonical hermite normal form of the supercell matrix
///
/// \returns Eigen::Matrix3i of H in canonical form
///
/// \param T A supercell matrix (Eigen::Matrix3i), such that S = U*T,
///          where S is the superlattice and U the unit lattice, as column
///          vector matrices
/// \param ref_lattice The lattice the transformation matrix T is meant to be
/// acting on
///
/// \param effective_pg Group of symmetry operations to use for determining
/// whether a canonical
///        HNF matrix has been found. This is probably the point group of a
///        structure or configuration
//         that has a lattice ref_lat (not to be confused with the point group
//         of ref_lat itself).
///
/// Canonical form is such that T is in Hermite normal form (as from
/// CASM::hermite_normal_form),
///  and the unrolled coefficients
/// \code
/// [a f e]
/// [0 b d] -> abcdef
/// [0 0 c]
/// \endcode
/// form the highest lexicographic order when considering equivalent
/// superlattices by point group operations.
///
/// - Equivalent superlattices can be obtained using point group operations: S'
/// = op*S = U*H*V,
///   where V is integer and has determinant +/- 1
/// - Substituting S = U*T, we have op*U*T = U*H*V.
/// - Or H*V = U.inverse*op*U*T, which is the hermite normal form of
/// U.inverse*op*U*T
/// - So T is canonical if it is in hermite normal form and for all point group
/// operations
///   it has a higher lexicographic order than the resulting H
///
///
/// \relatesalso Lattice
///
Eigen::Matrix3i canonical_hnf(const Eigen::Matrix3i &T,
                              const SymOpVector &effective_pg,
                              const Lattice &ref_lattice);
template <typename ExternSymGroupTypeIt>
Eigen::Matrix3i canonical_hnf(ExternSymGroupTypeIt begin,
                              ExternSymGroupTypeIt end,
                              const Eigen::Matrix3i &T,
                              const Lattice &ref_lattice) {
  return canonical_hnf(
      T, adapter::Adapter<SymOpVector, ExternSymGroupTypeIt>()(begin, end),
      ref_lattice);
}

}  // namespace xtal
}  // namespace CASM

#endif
