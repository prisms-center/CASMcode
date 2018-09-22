#include "casm/clex/NeighborList.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/PrimGrid.hh"

namespace CASM {

  /// \brief Return the weighting matrix used to define the canonical order
  PrimNeighborList::Matrix3Type PrimNeighborList::weight_matrix() const {
    return m_W;
  }

  /// \brief Expand the neighbor list to include the given UnitCellCoord
  void PrimNeighborList::expand(UnitCellCoord uccoord) {

    // save the old range
    Scalar prev_range = m_range;

    auto result = m_neighborhood.insert(uccoord.unitcell());

    // if a new UnitCell
    if(result.second) {
      // ensure all intermediate UnitCell are included
      _expand(prev_range);
    }
  }

  /// \brief Ensure that all intermediate UnitCell are included in our neighborhood
  void PrimNeighborList::_expand(Scalar prev_range) {

    // if some are new, we need to make sure we including all intermediate UnitCell in our neighborhood
    // get the score of the last element in the set
    m_range = _dist(*m_neighborhood.rbegin());

    // count over possible UnitCell, adding to neighborhood if prev_range < dist <= m_range
    // we're using a std::set for m_neighborhood, so it gets sorted

    // We want the bounding box that contains the ellipsoid defined by
    //   m_range = (i,j,k).transpose() * W * (i,j,k)
    //           = F.t * U.t * U * F
    //           = F.t * (M.inv).t * (M.inv) * F
    //         1 = F.t * (sqrt(m_range)*M.inv).t * (sqrt(m_range)*M.inv) * F
    // has maximum bb_end
    //   bb_end(i) = M.row(i).norm() (from considerations of an affine transformation M of unit sphere)
    //   where M = sqrt(m_range)*U.inv
    // and minimum
    //   bb_begin = -bb_end

    Eigen::MatrixXd M = sqrt(1.0 * m_range) * m_Uinv;

    auto d = [&](int i) {
      return std::lround(std::ceil(M.row(i).norm()));
    };

    VectorXType bb_end(3);
    bb_end << d(0), d(1), d(2);
    VectorXType bb_begin = -bb_end;
    VectorXType bb_incr = VectorXType::Constant(3, 1);

    typedef Counter< VectorXType,
            Scalar,
            Index,
            CASM_TMP::ParenthesesAccess<VectorXType, Scalar, Index> > VectorXCounter;

    VectorXCounter counter(bb_begin, bb_end, bb_incr);

    Scalar dist;
    for(; counter.valid(); ++counter) {
      dist = _dist(counter.current());
      if(prev_range < dist && dist <= m_range) {
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
  const PrimNeighborList::SublatIndices &PrimNeighborList::sublat_indices() const {
    return m_sublat_indices;
  }

  /// \brief Calculate A.transpose()*M*A
  PrimNeighborList::Scalar PrimNeighborList::_dist(const UnitCell &A) const {
    return A.transpose() * m_W * A;
  }

  /// \brief Return [r, i, j, k], where r = _dist(A)
  PrimNeighborList::VectorXType PrimNeighborList::_add_dist(const UnitCell &A) const {
    VectorXType vec(4);
    vec(0) = _dist(A);
    vec.segment(1, 3) = A;
    return vec;
  }

  /// \brief Convert [i,j,k] -> [r,i,j,k] and then lexicographically compare
  bool PrimNeighborList::_compare_unitcell(const UnitCell &A, const UnitCell &B) const {
    return _compare_vec(_add_dist(A), _add_dist(B));
  }

  /// \brief Lexicographical comparison
  bool PrimNeighborList::_compare_vec(const VectorXType &A, const VectorXType &B) {
    return std::lexicographical_compare(A.data(), A.data() + A.size(), B.data(), B.data() + B.size());
  }

  /// \brief Returns a NeighborList weighting matrix appropriate for a particular lattice
  ///
  /// Returns the integer weight matrix, M, closest to lat.transpose()*lat with
  /// all elements less than max_element_value.
  ///
  /// - C = L*F, where L is the lattice vectors as a column matrix
  /// - Want r = F.transpose*M*F such that r = C.transpose*C, where F is
  ///   fractional coordinates, and C is cartesian coordinates.
  /// - Find F.t*M*F = F.t*L.t*L*F -> W ~ L.tranpose*L
  /// - Approch, calculate L.tranpose*L, then divide by minimum element -> M'
  /// - consider integer n=1,2,...,n_last for all n such that all elements of n*M' are
  ///   less than 'max_element_value'
  /// - return the first n*M' that is integer. If none are integer, return n_last*M'
  PrimNeighborList::Matrix3Type PrimNeighborList::make_weight_matrix(const Eigen::Matrix3d lat_column_mat, Index max_element_value, double tol) {

    Eigen::Matrix3d W = lat_column_mat.transpose() * lat_column_mat;

    double min = std::numeric_limits<double>::max();
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        if(!almost_zero(W(i, j), tol) && std::abs(W(i, j)) < std::abs(min)) {
          min = W(i, j);
        }
      }
    }
    W /= min;

    // make positive definite
    Eigen::LLT<Eigen::MatrixXd> llt(W);
    if(llt.info() != Eigen::Success) {
      W = -W;
    }

    double n = 1.0;
    while(!is_integer(n * W, tol) &&
          ((n + 1.0)*W).cwiseAbs().maxCoeff() < max_element_value) {
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
  /// \param prim_grid A grid of unit cells describing the supercell this neighbor list pertains to
  /// \param prim_nlist_begin, prim_nlist_end Iterators over range [begin, end) of UnitCell that are
  ///        the neighbors of the origin UnitCell
  /// \param sublat_begin, sublat_end  Iterators over range [begin, end) of basis sites that should be
  ///        included as UnitCellCoord neighbors
  ///
  /// - The canonical order of UnitCellCoord is obtained by lexicographically sorting [r, i, j, k, b],
  ///   where r = (i,j,k).transpose() * W * (i,j,k).
  /// - The canonical order of UnitCell is obtained by lexicographically sorting [r, i, j, k]
  /// - The sublattice iterators enable restricting the neighbor list to only sites that have degrees of freedom
  ///
  SuperNeighborList::SuperNeighborList(const PrimGrid &prim_grid,
                                       const PrimNeighborList &prim_nlist) :
    m_prim_grid_size(prim_grid.size()),
    m_site(prim_grid.size()),
    m_unitcell(prim_grid.size()) {

    // use the PrimNeighborList to generate the UnitCell and Site indices for
    //   the neighbors of each UnitCell in the supercell

    UnitCellCoord bijk;

    // for each unit cell in the supercell
    for(Index i = 0; i < m_prim_grid_size; ++i) {

      // for each neighbor unitcell
      for(auto it = prim_nlist.begin(); it != prim_nlist.end(); ++it) {

        // get the neighbor unitcell's index
        size_type unitcell_index = prim_grid.find(prim_grid.unitcell(i) + *it);

        // store the unitcell index
        m_unitcell[i].push_back(unitcell_index);

        // calculate and store the site indices for all sites in the neighbor unitcell that are requested
        // - Depends on Configuration sites being stored in blocks by sublattice and unitcell indices
        //   determined by the PrimGrid ordering
        for(auto b_it = prim_nlist.sublat_indices().begin(); b_it != prim_nlist.sublat_indices().end(); ++b_it) {
          m_site[i].push_back((*b_it)*m_prim_grid_size + unitcell_index);
        }

      }
    }
  }

  /// \brief const Access the list of sites neighboring a particular unitcell
  const std::vector<SuperNeighborList::size_type> &SuperNeighborList::sites(size_type unitcell_index) const {
    return m_site[unitcell_index];
  }

  /// \brief const Access the list of unitcells neighboring a particular unitcell
  const std::vector<SuperNeighborList::size_type> &SuperNeighborList::unitcells(size_type unitcell_index) const {
    return m_unitcell[unitcell_index];
  }

  /// \brief Returns true if periodic images of the neighbor list overlap
  ///
  /// If periodic images of the neighborhood overlap, Clexulator 'delta' values will be incorrect.
  bool SuperNeighborList::overlaps() const {

    // there is an overlap if any of the neighboring unitcell indices is repeated
    // so sort and check if any two neighboring indices are the same

    std::vector<size_type> nlist = m_unitcell[0];
    std::sort(nlist.begin(), nlist.end());
    return std::adjacent_find(nlist.begin(), nlist.end()) != nlist.end();
  }

  /// \brief Clone
  std::unique_ptr<SuperNeighborList> SuperNeighborList::clone() const {
    return std::unique_ptr<SuperNeighborList>(new SuperNeighborList(*this));
  }

}
