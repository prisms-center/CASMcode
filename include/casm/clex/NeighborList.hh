#ifndef CASM_NeighborList_HH
#define CASM_NeighborList_HH

#include <memory>
#include <set>
#include <vector>

#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {

namespace xtal {
class UnitCell;
class UnitCellCoord;
class Superlattice;
}  // namespace xtal
using xtal::UnitCell;
using xtal::UnitCellCoord;

/** \ingroup ClexClex
 *  @{
 */

/// \brief The PrimNeighborList gives the coordinates of UnitCell that are
/// neighbors of the origin UnitCell
///
/// - The PrimNeighborList is constructed with a weighting matrix, W, that
/// defines the shape of neighborhood
/// - The canonical order of neighboring UnitCell is obtained by
/// lexicographically sorting [r, i, j, k],
///   where r = (i,j,k).transpose() * W * (i,j,k).
/// - The canonical order of UnitCellCoord is obtained by lexicographically
/// sorting [r, i, j, k, b],
///   where r = (i,j,k).transpose() * W * (i,j,k).
/// - The PrimNeighborList can be expanded as needed to increase the range of
/// the neighborhood without affecting
///   the order of the neighbors.
///
class PrimNeighborList {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef Index Scalar;
  typedef Eigen::Matrix3l Matrix3Type;
  typedef Eigen::Matrix<long, Eigen::Dynamic, 1> VectorXType;
  typedef std::set<UnitCell, std::function<bool(UnitCell, UnitCell)> >
      NeighborSet;
  typedef NeighborSet::size_type size_type;
  typedef NeighborSet::const_iterator const_iterator;
  typedef std::set<int> SublatIndices;

  /// \brief Constructor, specifying the weighting matrix to use, the indices
  /// of sublattices to include, and the total number of sublattices
  template <typename SublatIterator>
  PrimNeighborList(const Matrix3Type &W, SublatIterator begin,
                   SublatIterator end, size_type n_sublattices);

  /// \brief Return the weighting matrix W used to define the canonical order
  Matrix3Type weight_matrix() const;

  /// \brief Expand the neighbor list to include the given UnitCellCoord
  void expand(UnitCell const &uc);

  /// \brief Expand the neighbor list to include the given UnitCellCoord
  void expand(UnitCellCoord const &uccoord);

  /// \brief Expand the neighbor list to include the given range of
  /// UnitCellCoord
  template <typename UnitCellCoordIterator>
  void expand(UnitCellCoordIterator begin, UnitCellCoordIterator end);

  /// \brief size of the neighborhood of unit cells
  size_type size() const;

  /// \brief const_iterator over the neighborhood of unit cells
  const_iterator begin() const;

  /// \brief const_iterator over the neighborhood of unit cells
  const_iterator end() const;

  /// \brief const_iterator over the neighborhood of unit cells
  const_iterator cbegin() const;

  /// \brief const_iterator over the neighborhood of unit cells
  const_iterator cend() const;

  /// \brief pair of const_iterators over a range of indices of sublattices to
  /// include
  const SublatIndices &sublat_indices() const;

  /// \brief The total number of sublattices
  size_type n_sublattices() const;

  /// \brief Returns a NeighborList weighting matrix appropriate for a
  /// particular lattice
  static Matrix3Type make_weight_matrix(const Eigen::Matrix3d lat_column_mat,
                                        Index max_element_value, double tol);

  /// \brief Get neighborlist index of UnitCellCoord @param _ucc, expanding
  /// neighborhood if necessary
  Scalar neighbor_index(UnitCellCoord const &_ucc);

  /// \brief Get neighborlist indices of a collection of UnitCells, stored in
  /// @param _uc_container
  template <typename UnitCellCoordIterator>
  std::vector<Scalar> neighbor_indices(UnitCellCoordIterator _begin,
                                       UnitCellCoordIterator _end);

  /// \brief Clone
  std::unique_ptr<PrimNeighborList> clone() const;

 private:
  /// \brief Get neighborlist index of UnitCellCoord @param _ucc, without
  /// expanding neighborhood
  Scalar _neighbor_index(UnitCellCoord const &_ucc) const;

  /// \brief Expand the neighbor list to include the given UnitCell, but do not
  /// do additional updates returns true if added UnitCell is new
  bool _expand(UnitCell const &uc);

  /// \brief Expand the neighbor list to include the given UnitCellCoord, but do
  /// not do additional updates returns true if added UnitCellCoord is new
  bool _expand(UnitCellCoord const &uccoord);

  /// \brief Ensure that all intermediate UnitCell are included in our
  /// neighborhood
  void _expand(Scalar prev_range);

  /// \brief Calculate A.transpose()*M*A
  Scalar _dist(const UnitCell &A) const;

  /// \brief Return [r, i, j, k], where r = _dist(A)
  VectorXType _add_dist(const UnitCell &A) const;

  /// \brief Convert [i,j,k] -> [r,i,j,k] and then lexicographically compare
  bool _compare_unitcell(const UnitCell &A, const UnitCell &B) const;

  /// \brief Lexicographical comparison
  static bool _compare_vec(const VectorXType &A, const VectorXType &B);

  /// \brief Weighting matrix
  ///
  /// - The canonical order of UnitCellCoord is obtained by lexicographically
  /// sorting [r, i, j, k, b],
  ///   where r = (i,j,k).transpose() * W * (i,j,k).
  ///
  Matrix3Type m_W;

  /// \brief Cholesky decomposition: W = U.transpose()*U, useful for getting
  /// bounding box
  Eigen::MatrixXd m_Uinv;

  /// \brief the neighborhood of [i, j, k]
  NeighborSet m_neighborhood;

  /// \brief the neighborhood, m_neighborhood, contains all UnitCell with r <=
  /// m_range
  Scalar m_range;

  /// \brief the indices of sublattices that should be included
  SublatIndices m_sublat_indices;

  /// \brief the total number of sublattices
  size_type m_n_sublattices;
};

/// SuperNeighborList, linear indices of neighboring sites and unit cells
///
/// The SuperNeighborList takes the ordering of unit cells neighboring the
/// origin unit cell from the PrimNeighborList and uses it to construct lists of
/// linear indices of neighboring unit cells and neighboring sites for each unit
/// cell in a particular supercell.
///
/// Neighbor list index naming conventions:
/// - "unitcell_index", a linear index given to all unit cell "within" a
/// supercell. Convert
///   between unitcell_index and xtal::UnitCell with
///   xtal::UnitCellIndexConverter.
/// - "site_index", a linear index given to all sites "within" a supercell.
/// Convert between
///   site_index and xtal::UnitCellCoord with xtal::UnitCellCoordIndexConverter.
///
class SuperNeighborList {
 public:
  typedef Index size_type;

  /// Constructor
  SuperNeighborList(Eigen::Matrix3l const &transformation_matrix_to_super,
                    PrimNeighborList const &prim_nlist);

  /// Constructor
  SuperNeighborList(const xtal::Superlattice &prim_grid,
                    const PrimNeighborList &prim_nlist);

  // --- Keep inlined functions inline for most efficient use  ---

  size_type n_unitcells() const { return m_prim_grid_size; }

  /// Get unitcell_index from site_index
  size_type unitcell_index(size_type site_index) const {
    return site_index % m_prim_grid_size;
  }

  /// Get sublattice index from site_index
  size_type sublat_index(size_type site_index) const {
    return site_index / m_prim_grid_size;
  }

  /// \brief Get neighbor index from site_index (use for periodic Clexulator
  /// point corr evalulations `neighbor_index` argument)
  ///
  /// - Returns -1 if site is on a sublattice that is not included in the
  /// neighbor list
  int neighbor_index(size_type site_index) const {
    return m_site_index_to_neighbor_index[site_index];
  }

  /// \brief const Access the list of sites neighboring a particular unit cell
  const std::vector<size_type> &sites(size_type unitcell_index) const {
    return m_site[unitcell_index];
  }

  /// \brief const Access the list of unitcells neighboring a particular unit
  /// cell
  const std::vector<size_type> &unitcells(size_type unitcell_index) const {
    return m_unitcell[unitcell_index];
  }

  /// \brief Returns true if periodic images of the neighbor list overlap
  ///
  /// If periodic images of the neighborhood overlap, Clexulator 'delta' values
  /// will be incorrect.
  bool overlaps() const { return m_overlaps; }

  /// \brief Clone
  std::unique_ptr<SuperNeighborList> clone() const;

 private:
  /// \brief store prim grid size for site index -> unitcell index conversion
  /// unitcell_index = site_index % m_prim_grid_size
  size_type m_prim_grid_size;

  /// \brief m_site[unitcell_index][neighbor site index]
  ///
  /// - Configuration sites are ordered in blocks corresponding to each
  /// sublattice, b. Neighbors are
  ///   specific only to a unitcell, so the neighbor list for site index s and s
  ///   + n*scel.volume() are identical.
  /// - So m_site.size() == m_prim_grid_size
  ///
  std::vector<std::vector<size_type> > m_site;

  /// \brief m_unitcell[unitcell_index][neighbor unitcell index]
  std::vector<std::vector<size_type> > m_unitcell;

  /// \brief neighbor index =
  ///     m_linear_site_index_to_neighbor_index[linear site index]
  std::vector<int> m_site_index_to_neighbor_index;

  /// \brief True if periodic images of the neighbor list overlap
  bool m_overlaps;
};

/// \brief Constructor
///
/// \param W A positive-definite weighting matrix for calculating the weighted
///          squared L2 norm of the unitcell indices.
///
/// - The canonical order of UnitCell is obtained by lexicographically sorting
/// [r, i, j, k],
///   where r = (i,j,k).transpose() * W * (i,j,k).
/// - The canonical order of UnitCellCoord is obtained by lexicographically
/// sorting [r, i, j, k, b],
///   where r = (i,j,k).transpose() * W * (i,j,k).
///
template <typename SublatIterator>
PrimNeighborList::PrimNeighborList(const Matrix3Type &W, SublatIterator begin,
                                   SublatIterator end, size_type n_sublattices)
    : m_W(W),
      m_neighborhood(std::bind(&PrimNeighborList::_compare_unitcell, this,
                               std::placeholders::_1, std::placeholders::_2)),
      m_range(0),
      m_sublat_indices(begin, end),
      m_n_sublattices(n_sublattices) {
  // throw if W is not positive definite
  Eigen::LLT<Eigen::MatrixXd> llt(W.cast<double>());
  if (llt.info() != Eigen::Success) {
    std::cerr << "Error constructing PrimNeighborList: weight matrix is not "
                 "positive definite."
              << std::endl;
    std::cerr << "weight matrix: \n" << W << std::endl;
    throw std::runtime_error(
        "Error constructing PrimNeighborList: weight matrix is not positive "
        "definite.");
  }

  // store U.inverse()
  Eigen::MatrixXd U = llt.matrixU();
  m_Uinv = U.inverse();

  // add origin unit cell to neighborhood
  m_neighborhood.insert(UnitCell(0, 0, 0));
}

/// \brief Expand the neighbor list to include the given range of UnitCellCoord
///
/// - The canonical order of UnitCellCoord is obtained by lexicographically
/// sorting [r, i, j, k, b],
///   where r = (i,j,k).transpose() * W * (i,j,k).
/// - The canonical order of UnitCell is obtained by lexicographically sorting
/// [r, i, j, k]
/// - After determining the sorted order of the requested UnitCell, all
/// intermediate UnitCell
///   will also be added so that indices into the neighbor list are always
///   constant and in canonical order
///
template <typename UnitCellCoordIterator>
void PrimNeighborList::expand(UnitCellCoordIterator begin,
                              UnitCellCoordIterator end) {
  // save the old range
  Scalar prev_range = m_range;

  bool any_new = false;
  for (auto it = begin; it != end; ++it) {
    any_new = _expand(*it) || any_new;
  }

  if (!any_new) {
    return;
  }

  // otherwise, ensure all intermediate UnitCell are included
  _expand(prev_range);
}

/// \brief Get neighborlist indices of a of UnitCells, passed by begin and end
/// iterator
template <typename UnitCellCoordIterator>
std::vector<PrimNeighborList::Scalar> PrimNeighborList::neighbor_indices(
    UnitCellCoordIterator _begin, UnitCellCoordIterator _end) {
  expand(_begin, _end);
  std::vector<Scalar> result;
  std::transform(_begin, _end, std::back_inserter(result),
                 [this](UnitCellCoord const &A) -> Scalar {
                   return this->_neighbor_index(A);
                 });
  return result;
}

/** @} */
}  // namespace CASM

#endif
