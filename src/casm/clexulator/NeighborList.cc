#include "casm/clexulator/NeighborList.hh"

#include "casm/container/Counter.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace clexulator {

/// \brief Return the weighting matrix used to define the canonical order
PrimNeighborList::Matrix3Type PrimNeighborList::weight_matrix() const {
  return m_W;
}

/// \brief Expand the neighbor list to include the given UnitCellCoord
void PrimNeighborList::expand(xtal::UnitCellCoord const &uccoord) {
  expand(uccoord.unitcell());
}

/// \brief Expand the neighbor list to include the given UnitCellCoord
void PrimNeighborList::expand(xtal::UnitCell const &uc) {
  // save the old range
  Scalar prev_range = m_range;

  auto result = m_neighborhood.insert(uc);

  // if a new UnitCell
  if (result.second) {
    // ensure all intermediate UnitCell are included
    _expand(prev_range);
  }
}

/// \brief Expand the neighbor list to include the given UnitCellCoord
bool PrimNeighborList::_expand(xtal::UnitCellCoord const &uccoord) {
  return _expand(uccoord.unitcell());
}

/// \brief Expand the neighbor list to include the given UnitCellCoord
bool PrimNeighborList::_expand(xtal::UnitCell const &uc) {
  return m_neighborhood.insert(uc).second;
}

/// \brief Ensure that all intermediate UnitCell are included in our
/// neighborhood
void PrimNeighborList::_expand(Scalar prev_range) {
  // if some are new, we need to make sure we including all intermediate
  // UnitCell in our neighborhood get the score of the last element in the set
  m_range = _dist(*m_neighborhood.rbegin());

  // count over possible UnitCell, adding to neighborhood if prev_range < dist
  // <= m_range we're using a std::set for m_neighborhood, so it gets sorted

  // We want the bounding box that contains the ellipsoid defined by
  //   m_range = (i,j,k).transpose() * W * (i,j,k)
  //           = F.t * U.t * U * F
  //           = F.t * (M.inv).t * (M.inv) * F
  //         1 = F.t * (sqrt(m_range)*M.inv).t * (sqrt(m_range)*M.inv) * F
  // has maximum bb_end
  //   bb_end(i) = M.row(i).norm() (from considerations of an affine
  //   transformation M of unit sphere) where M = sqrt(m_range)*U.inv
  // and minimum
  //   bb_begin = -bb_end

  Eigen::MatrixXd M = sqrt(1.0 * m_range) * m_Uinv;

  auto d = [&](int i) { return std::lround(std::ceil(M.row(i).norm())); };

  VectorXType bb_end(3);
  bb_end << d(0), d(1), d(2);
  VectorXType bb_begin = -bb_end;
  VectorXType bb_incr = VectorXType::Constant(3, 1);

  typedef Counter<VectorXType, Scalar, Index,
                  CASM_TMP::ParenthesesAccess<VectorXType, Scalar, Index> >
      VectorXCounter;

  VectorXCounter counter(bb_begin, bb_end, bb_incr);

  Scalar dist;
  for (; counter.valid(); ++counter) {
    dist = _dist(counter.current());
    if (prev_range < dist && dist <= m_range) {
      m_neighborhood.insert(counter.current());
    }
  }
}

/// \brief size of the neighborhood of unit cells
PrimNeighborList::size_type PrimNeighborList::size() const {
  return m_neighborhood.size();
}

/// \brief const_iterator over the neighborhood of unit cells
PrimNeighborList::const_iterator PrimNeighborList::begin() const {
  return m_neighborhood.cbegin();
}

/// \brief const_iterator over the neighborhood of unit cells
PrimNeighborList::const_iterator PrimNeighborList::end() const {
  return m_neighborhood.cend();
}

/// \brief const_iterator over the neighborhood of unit cells
PrimNeighborList::const_iterator PrimNeighborList::cbegin() const {
  return m_neighborhood.cbegin();
}

/// \brief const_iterator over the neighborhood of unit cells
PrimNeighborList::const_iterator PrimNeighborList::cend() const {
  return m_neighborhood.cend();
}

/// \brief const_iterator over the neighborhood of unit cells
const PrimNeighborList::SublatIndices &PrimNeighborList::sublat_indices()
    const {
  return m_sublat_indices;
}

/// \brief The total number of sublattices
PrimNeighborList::size_type PrimNeighborList::n_sublattices() const {
  return m_n_sublattices;
}

/// \brief Get neighborlist index of UnitCellCoord @param _ucc, expanding
/// neighborhood if necessary
PrimNeighborList::Scalar PrimNeighborList::neighbor_index(
    xtal::UnitCellCoord const &_ucc) {
  expand(_ucc);
  return _neighbor_index(_ucc);
}

/// \brief Get neighborlist index of UnitCellCoord @param _ucc, without
/// expanding neighborhood
PrimNeighborList::Scalar PrimNeighborList::_neighbor_index(
    xtal::UnitCellCoord const &_ucc) const {
  Scalar uc_ind(find_index(m_neighborhood, _ucc.unitcell()));
  Scalar sublat_dist(find_index(sublat_indices(), _ucc.sublattice()));

  return uc_ind * sublat_indices().size() + sublat_dist;
}

/// \brief Calculate A.transpose()*M*A
PrimNeighborList::Scalar PrimNeighborList::_dist(
    xtal::UnitCell const &A) const {
  return A.transpose() * m_W * A;
}

/// \brief Return [r, i, j, k], where r = _dist(A)
PrimNeighborList::VectorXType PrimNeighborList::_add_dist(
    xtal::UnitCell const &A) const {
  VectorXType vec(4);
  vec(0) = _dist(A);
  vec.segment(1, 3) = A;
  return vec;
}

/// \brief Convert [i,j,k] -> [r,i,j,k] and then lexicographically compare
bool PrimNeighborList::_compare_unitcell(xtal::UnitCell const &A,
                                         xtal::UnitCell const &B) const {
  return _compare_vec(_add_dist(A), _add_dist(B));
}

/// \brief Lexicographical comparison
bool PrimNeighborList::_compare_vec(VectorXType const &A,
                                    VectorXType const &B) {
  return std::lexicographical_compare(A.data(), A.data() + A.size(), B.data(),
                                      B.data() + B.size());
}

/// \brief Returns a NeighborList weighting matrix appropriate for a particular
/// lattice
///
/// Returns the integer weight matrix, M, closest to lat.transpose()*lat with
/// all elements less than max_element_value.
///
/// - C = L*F, where L is the lattice vectors as a column matrix
/// - Want r = F.transpose*M*F such that r = C.transpose*C, where F is
///   fractional coordinates, and C is cartesian coordinates.
/// - Find F.t*M*F = F.t*L.t*L*F -> W ~ L.tranpose*L
/// - Approch, calculate L.tranpose*L, then divide by minimum element -> M'
/// - consider integer n=1,2,...,n_last for all n such that all elements of n*M'
/// are
///   less than 'max_element_value'
/// - return the first n*M' that is integer. If none are integer, return
/// n_last*M'
PrimNeighborList::Matrix3Type PrimNeighborList::make_weight_matrix(
    const Eigen::Matrix3d lat_column_mat, Index max_element_value, double tol) {
  Eigen::Matrix3d W = lat_column_mat.transpose() * lat_column_mat;

  double min = std::numeric_limits<double>::max();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (!almost_zero(W(i, j), tol) && std::abs(W(i, j)) < std::abs(min)) {
        min = W(i, j);
      }
    }
  }
  W /= min;

  // make positive definite
  Eigen::LLT<Eigen::MatrixXd> llt(W);
  if (llt.info() != Eigen::Success) {
    W = -W;
  }

  double n = 1.0;
  while (!is_integer(n * W, tol) &&
         ((n + 1.0) * W).cwiseAbs().maxCoeff() < max_element_value) {
    n += 1.0;
  }
  return lround(n * W);
}

/// \brief Clone
std::unique_ptr<PrimNeighborList> PrimNeighborList::clone() const {
  return std::unique_ptr<PrimNeighborList>(new PrimNeighborList(*this));
}

/// \brief Constructor
///
/// \param transformation_matrix_to_super The transformation matrix, T, from
/// primitive lattice
///        vectors to supercell lattice vectors: `supercell_lat_column_mat = T *
///        prim_lat_column_mat`.
/// \param prim_nlist A reference to a PrimNeighborList defining the
/// neighborhood for the origin
///        unit cell. It is only used in the constructor, any changes to the
///        PrimNeighborList will not be reflected in this SuperNeighborList.
///
/// - The canonical order of UnitCellCoord is obtained by lexicographically
/// sorting [r, i, j, k, b],
///   where r = (i,j,k).transpose() * W * (i,j,k).
/// - The canonical order of UnitCell is obtained by lexicographically sorting
/// [r, i, j, k]
/// - The sublattice iterators enable restricting the neighbor list to only
/// sites that have degrees of freedom
///
SuperNeighborList::SuperNeighborList(
    Eigen::Matrix3l const &transformation_matrix_to_super,
    const PrimNeighborList &prim_nlist) {
  xtal::UnitCellIndexConverter ijk_index_converter{
      transformation_matrix_to_super};
  // confusingly, `ijk_index_converter.total_sites()` is number of unitcells in
  // the supercell
  m_prim_grid_size = ijk_index_converter.total_sites();
  m_site.resize(m_prim_grid_size);
  m_unitcell.resize(m_prim_grid_size);

  // use the PrimNeighborList to generate the UnitCell and Site indices for
  //   the neighbors of each UnitCell in the supercell

  // for each unit cell in the supercell
  for (Index i = 0; i < m_prim_grid_size; ++i) {
    // for each neighbor unitcell
    for (auto it = prim_nlist.begin(); it != prim_nlist.end(); ++it) {
      // get the neighbor unitcell's index
      xtal::UnitCell neighbor_unitcell = ijk_index_converter(i) + *it;
      size_type neighbor_unitcell_index =
          ijk_index_converter(neighbor_unitcell);

      // store the unitcell index
      m_unitcell[i].push_back(neighbor_unitcell_index);

      // calculate and store the site indices for all sites in the neighbor
      // unitcell that are requested
      // - Depends on Configuration sites being stored in blocks by sublattice
      // and unitcell indices
      //   determined by the UnitCellCoordIndexConverter ordering
      for (auto b_it = prim_nlist.sublat_indices().begin();
           b_it != prim_nlist.sublat_indices().end(); ++b_it) {
        m_site[i].push_back((*b_it) * m_prim_grid_size +
                            neighbor_unitcell_index);
      }
    }
  }

  // populate m_site_index_to_neighbor_index
  m_site_index_to_neighbor_index.clear();
  auto sublat_indices_begin = prim_nlist.sublat_indices().begin();
  auto sublat_indices_end = prim_nlist.sublat_indices().end();
  for (Index b = 0; b < prim_nlist.n_sublattices(); ++b) {
    auto it = prim_nlist.sublat_indices().find(b);
    int neighbor_index = -1;
    if (it != sublat_indices_end) {
      neighbor_index = std::distance(sublat_indices_begin, it);
    }
    for (Index i = 0; i < m_prim_grid_size; ++i) {
      m_site_index_to_neighbor_index.push_back(neighbor_index);
    }
  }

  // check for overlap of periodic images
  //
  // there is an overlap if any of the neighboring unitcell indices is repeated
  // so sort and check if any two neighboring indices are the same
  std::vector<size_type> nlist = m_unitcell[0];
  std::sort(nlist.begin(), nlist.end());
  m_overlaps = std::adjacent_find(nlist.begin(), nlist.end()) != nlist.end();
}

/// \brief Clone
std::unique_ptr<SuperNeighborList> SuperNeighborList::clone() const {
  return std::unique_ptr<SuperNeighborList>(new SuperNeighborList(*this));
}

}  // namespace clexulator

}  // namespace CASM
