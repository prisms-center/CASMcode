#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  //********************************************************************
  /**
   * This function generates a grid of points between max_radius and
   * min_radius. Additionally, it also fills up the points with a basis
   */
  //********************************************************************

  template<typename CoordType, typename CoordType2>
  Array<CoordType> Lattice::gridstruc_build(double max_radius, double min_radius, Array<CoordType> basis, CoordType2 lat_point) {
    Eigen::Vector3i dim;
    dim = enclose_sphere(max_radius);
    EigenCounter<Eigen::Vector3i > grid_count(-dim, dim, Eigen::Vector3i(1));
    double min_dist, dist;
    Array<CoordType> gridstruc;
    Eigen::Vector3i temp;

    do {
      lat_point(FRAC) = grid_count();

      for(Index i = 0; i < basis.size(); i++) {
        CoordType tatom(basis[i] + lat_point);
        //get distance to closest basis site in the unit cell at the origin

        min_dist = 1e20;
        for(Index j = 0; j < basis.size(); j++) {
          dist = tatom.dist(basis[j]);
          if(dist < min_dist)
            min_dist = dist;
        }
        if(min_dist < min_radius) {
          continue;
        }
        if(min_dist < max_radius) {
          gridstruc.push_back(tatom);
          //          std::cout<<"tatom"<<tatom<<"\t Min Dist"<<min_dist<<"\n";
        }
      }
    }
    while(++grid_count);

    return gridstruc;
  }


  ///\brief returns Lattice that is smallest possible supercell of all input Lattice
  ///
  /// If SymOpIterator are provided they are applied to each Lattice in an attempt
  /// to find the smallest possible superdupercell of all symmetrically transformed Lattice
  template<typename LatIterator, typename SymOpIterator>
  Lattice superdupercell(LatIterator begin,
                         LatIterator end,
                         SymOpIterator op_begin,
                         SymOpIterator op_end) {

    Lattice best = *begin;
    for(auto it = ++begin; it != end; ++it) {
      Lattice tmp_best = superdupercell(best, *it);
      for(auto op_it = op_begin; op_it != op_end; ++op_it) {
        Lattice test = superdupercell(best, copy_apply(*op_it, *it));
        if(std::abs(volume(test)) < std::abs(volume(tmp_best))) {
          tmp_best = test;
        }
      }
      best = tmp_best;
    }
    return best;
  }

  /// \brief Output the SymOp that leave this lattice invariant
  template<typename SymOpIterator, typename SymOpOutputIterator>
  SymOpOutputIterator Lattice::find_invariant_subgroup(
    SymOpIterator begin,
    SymOpIterator end,
    SymOpOutputIterator result,
    double pg_tol) const {

    LatticeIsEquivalent is_equiv(*this, pg_tol);
    return std::copy_if(begin, end, result, is_equiv);
  }


  /// Check if there is a symmetry operation, op, and transformation matrix T,
  ///   such that scel is a supercell of the result of applying op to unit
  ///
  /// \returns pair corresponding to first successful op and T, or with op=end if not successful
  template<typename Object, typename OpIterator>
  std::pair<OpIterator, Eigen::MatrixXi> is_supercell(
    const Object &scel,
    const Object &unit,
    OpIterator begin,
    OpIterator end,
    double tol) {

    std::pair<bool, Eigen::MatrixXi> res;
    for(auto it = begin; it != end; ++it) {
      res = is_supercell(scel, copy_apply(*it, unit), tol);
      if(res.first) {
        return std::make_pair(it, res.second);
      }
    }
    return std::make_pair(end, res.second);
  }
}

