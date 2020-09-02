#include <numeric>
#include "casm/misc/CASM_math.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/external/Eigen/CASM_AddOns"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/symmetry/VectorSymCompare_impl.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/misc/algorithm.hh"
#include "casm/container/Counter.hh"

#include "casm/casm_io/Log.hh" // TODO: replace error log with exceptions

namespace CASM {
  namespace Local {

    struct IrrepCompare {
      bool operator()(SymRepTools::IrrepInfo const &irrep_a, SymRepTools::IrrepInfo const &irrep_b) const {
        Eigen::VectorXcd const &a = irrep_a.characters;
        Eigen::VectorXcd const &b = irrep_b.characters;

        typedef std::complex<double> C;

        // Check if a and b are identity irrep, identity always compares less than other irreps
        bool a_id = almost_equal(a[0], C(1., 0.)) && almost_equal(C(a.sum()), C(a.size(), 0.));
        bool b_id = almost_equal(b[0], C(1., 0.)) && almost_equal(C(b.sum()), C(b.size(), 0.));

        if(a_id != b_id)
          return a_id;
        else if(a_id) {
          // Index is used to break ties if they are both identity
          return irrep_a.index < irrep_b.index;
        }

        // Low-dimensional irreps come before higher dimensional
        if(!almost_equal(a[0], b[0]))
          return a[0].real() < b[0].real();

        // 'gerade' irreps come before 'ungerade' irreps
        // This check may need to be improved to know whether inversion is actually present
        bool a_g = almost_equal(a[0], a[a.size() - 1]);
        bool b_g = almost_equal(b[0], b[b.size() - 1]);
        if(a_g != b_g)
          return a_g;

        // Finally, compare lexicographically (real first, then imag)
        for(Index i = 0; i < a.size(); ++i) {
          if(!almost_equal(a[i].real(), b[i].real()))
            return a[i].real() > b[i].real();
        }

        if(irrep_a.complex || irrep_b.complex) {
          for(Index i = 0; i < a.size(); ++i) {
            if(!almost_equal(a[i].imag(), b[i].imag()))
              return a[i].imag() > b[i].imag();
          }
        }

        // Index is used to break ties
        return irrep_a.index < irrep_b.index;
      }
    };

    //*******************************************************************************************
    template<typename T>
    struct _RealType;

    template<typename T>
    using _Real = typename _RealType<T>::Type;

    template<typename T>
    struct _RealType<std::vector<T> > {
      using Type = std::vector<_Real<T> >;
    };

    template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
    struct _RealType<Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime> > {
      using Type = Eigen::Matrix<double, RowsAtCompileTime, ColsAtCompileTime>;
    };


    template<typename Derived>
    _Real<Derived> _real(Eigen::MatrixBase<Derived> const &mat) {
      return mat.real().template cast<double>();
    }

    template<typename T>
    _Real<std::vector<T> > _real(std::vector<T> const &vec) {

      std::vector<_Real<T> > result;
      result.reserve(vec.size());
      for(T const &el : vec) {
        result.emplace_back(_real(el));
      }

      return result;
    }

    //*******************************************************************************************
    /// \brief Matrix such that result(i,j) is sum of squares over 'p' of  op_rep[p].matrix(i,j). 'p' only spans operations in head_group
    /// The resulting matrix reveals the block_diagonalization of SymGroupRep _rep for the subset of operations contained in head_group
    static Eigen::MatrixXd _block_shape_matrix(SymGroupRep const &_rep, SymGroup const &head_group) {
      if(!_rep.size() || !head_group.size() || !_rep.MatrixXd(head_group[0]))
        return Eigen::MatrixXd();

      Eigen::MatrixXd block_shape(_rep.MatrixXd(head_group[0])->cwiseProduct(*_rep.MatrixXd(head_group[0])));

      for(Index i = 1; i < head_group.size(); i++) {
        block_shape += _rep.MatrixXd(head_group[i])->cwiseProduct(*_rep.MatrixXd(head_group[i]));
      }

      return block_shape;
    }

    //*******************************************************************************************

    //assumes that representation is real-valued and irreducible
    static Eigen::MatrixXcd _irrep_symmetrizer_from_directions(multivector<Eigen::VectorXcd>::X<2> const &special_directions,
                                                               Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                                                               double vec_compare_tol) {

      // Four strategies, in order of desparation
      // 1) find a spanning set of orthogonal axes within a single orbit of special directions
      // 2) find a spanning set of orthogonal axes within the total set of special directions
      // 3) perform qr decomposition on lowest-multiplicity orbit to find spanning set of axes
      // 4) find column-echelon form of _subspace matrix to get a sparse/pretty set of axes
      Eigen::MatrixXcd result;
      Index dim = _subspace.cols();
      Index min_mult = 10000;
      bool orb_orthog = false;
      bool tot_orthog = false;

      Eigen::MatrixXcd axes, orb_axes, tot_axes;
      Index tot_col(0);
      tot_axes.setZero(_subspace.rows(), dim);

      for(auto const &orbit : special_directions) {
        // Strategy 1
        if((orb_orthog && orbit.size() < min_mult) || !orb_orthog) {
          orb_axes.setZero(_subspace.rows(), dim);
          Index col = 0;
          for(auto const &el : orbit) {
            if(almost_zero((el.adjoint()*orb_axes).eval(), vec_compare_tol)) {
              if(col < orb_axes.cols())
                orb_axes.col(col++) = el;
              else {
                std::stringstream errstr;
                errstr << "Error in irrep_symmtrizer_from_directions(). Constructing coordinate axes from special directions of space spanned by row vectors:\n "
                       << _subspace.transpose()
                       << "\nAxes collected thus far are the row vectors:\n" << orb_axes.transpose()
                       << "\nBut an additional orthogonal row vector has been found:\n" << el.transpose()
                       << "\nindicating that subspace matrix is malformed.";
                throw std::runtime_error(errstr.str());
              }
            }
          }
          if(col == dim) {
            //std::cout << "PLAN A-- col: " << col << "; dim: " << dim << "; min_mult: " << min_mult << "; orthog: " << orb_orthog << ";\naxes: \n" << axes << "\norb_axes: \n" << orb_axes << "\n\n";
            orb_orthog = true;
            min_mult = orbit.size();
            axes = orb_axes;
          }
        }

        // Greedy(ish) implementation of strategy 2 -- may not find a solution, even if it exists
        if(!orb_orthog && !tot_orthog) {
          for(auto const &el : orbit) {
            if(almost_zero((el.adjoint()*tot_axes).eval(), vec_compare_tol)) {
              if(tot_col < tot_axes.cols())
                tot_axes.col(tot_col++) = el;
              else {
                std::stringstream errstr;
                errstr << "Error in irrep_symmtrizer_from_directions(). Constructing coordinate axes from special directions of space spanned by row vectors:\n "
                       << _subspace.transpose()
                       << "\nAxes collected thus far are the row vectors:\n" << tot_axes.transpose()
                       << "\nBut an additional orthogonal row vector has been found:\n" << el.transpose()
                       << "\nindicating that subspace matrix is malformed.";
                throw std::runtime_error(errstr.str());
              }
            }
          }
          if(tot_col == dim) {
            //std::cout << "PLAN B-- col: " << tot_col << "; dim: " << dim << "; min_mult: " << min_mult << "; orthog: " << tot_orthog << ";\naxes: \n" << axes << "\ntot_axes: \n" << tot_axes << "\n\n";
            tot_orthog = true;
            axes = tot_axes;
          }
        }

        // Strategy 3
        if(!orb_orthog && !tot_orthog && orbit.size() < min_mult) {
          orb_axes.setZero(_subspace.rows(), orbit.size());
          for(Index col = 0; col < orbit.size(); ++col) {
            orb_axes.col(col) = orbit[col];
          }
          //std::cout << "PLAN C--  dim: " << dim << "; min_mult: " << min_mult << "; orthog: " << orb_orthog << ";\naxes: \n" << axes;
          min_mult = orbit.size();
          axes = Eigen::MatrixXcd(orb_axes.colPivHouseholderQr().matrixQ()).leftCols(dim);
          //std::cout << "\norb_axes: \n" << orb_axes << "\n\n";
        }
      }
      //std::cout << "axes: \n" << axes << "\n";
      if(axes.cols() == 0) {
        //std::cout << "Plan D\n";
        // SubspaceSymCompare doesn't actually use _rep
        result = _subspace.colPivHouseholderQr().solve(representation_prepare_impl(_subspace, vec_compare_tol));
      }
      else {
        // Strategy 4
        result = _subspace.colPivHouseholderQr().solve(axes);
      }

      //std::cout << "result: \n" << result << "\n"
      //<< "_subspace*result: \n" << _subspace *result << "\n";
      return result;//.transpose();
    }

    //*******************************************************************************************

    static SymRepTools::IrrepWedge _wedge_from_pseudo_irrep(SymRepTools::IrrepInfo const &irrep,
                                                            SymGroupRep const &_rep,
                                                            SymGroup const &head_group) {
      Eigen::MatrixXd t_axes = irrep.trans_mat.transpose().real();
      Eigen::MatrixXd axes = representation_prepare_impl(t_axes, TOL);
      Eigen::VectorXd v = axes.col(0);
      Eigen::VectorXd vbest;

      axes.setZero(irrep.vector_dim(), irrep.irrep_dim());

      axes.col(0) = v;

      for(Index i = 1; i < axes.cols(); ++i) {
        double bestproj = -1;
        for(SymOp const &op : head_group) {
          v = (*_rep.MatrixXd(op)) * axes.col(0);
          //std::cout << "v: " << v.transpose() << std::endl;
          bool skip_op = false;
          for(Index j = 0; j < i; ++j) {
            if(almost_equal(v, axes.col(j))) {
              skip_op = true;
              break;
            }
          }
          if(skip_op)
            continue;

          //std::cout << "bproj: " << bestproj << "  proj: " << (v.transpose()*axes).transpose() << std::endl;
          if(bestproj < (v.transpose()*axes).sum()) {
            bestproj = (v.transpose() * axes).sum();
            vbest = v;
          }
        }
        axes.col(i) = vbest;
      }
      return SymRepTools::IrrepWedge(irrep, axes);
    }

    //*******************************************************************************************

    bool _rep_check(SymGroupRep const &_rep, SymGroup const &head_group, bool verbose) {
      bool passed = true;
      for(Index ns = 0; ns < head_group.size(); ns++) {
        Eigen::MatrixXd tmat = *(_rep.MatrixXd(head_group[ns]));
        if(!almost_equal<double>((tmat * tmat.transpose()).trace(), tmat.cols())) {
          passed = false;
          if(verbose)
            std::cout << "  Representation failed full-rank/orthogonality check!\n"
                      << "  ns: " << ns << "\n"
                      << "  Matrix: \n" << tmat << "\n\n"
                      << "  Matrix*transpose: \n" << tmat.transpose()*tmat << "\n\n";
        }
        for(Index ns2 = ns; ns2 < head_group.size(); ns2++) {
          auto prod = (*(_rep.MatrixXd(head_group[ns]))) * (*(_rep.MatrixXd(head_group[ns2])));
          Index iprod = head_group[ns].ind_prod(head_group[ns2]);
          double norm = (prod - * (_rep.MatrixXd(iprod))).norm();
          if(!almost_zero(norm)) {
            passed = false;
            if(verbose)
              std::cout << "  Representation failed multiplication check!\n"
                        << "  ns: " << ns << "  ns2: " << ns2 << " iprod: " << iprod << " NO MATCH\n"
                        << "  prod: \n" << prod << "\n\n"
                        << "  mat(ns): \n" << *(_rep.MatrixXd(ns)) << "\n\n"
                        << "  mat(ns2): \n" << *(_rep.MatrixXd(ns2)) << "\n\n"
                        << "  mat(iprod): \n" << *(_rep.MatrixXd(iprod)) << "\n\n";
          }
        }
      }
      return passed;
    }

    template<typename T>
    SymRepTools::IrrepInfo _subspace_to_full_space(SymRepTools::IrrepInfo const &irrep,
                                                   Eigen::MatrixBase<T> const &subspace) {
      SymRepTools::IrrepInfo result(irrep);
      result.trans_mat = irrep.trans_mat * subspace.adjoint().template cast<std::complex<double> >();

      result.directions.clear();
      for (const auto& direction_orbit : irrep.directions){
          std::vector<Eigen::VectorXd> new_orbit;
          new_orbit.reserve(direction_orbit.size());
          for (const auto& directions: direction_orbit){
              new_orbit.push_back(subspace*directions);
          }
          result.directions.push_back(std::move(new_orbit));
      }
      return result;
    }

  }
}

namespace CASM {
  namespace SymRepTools {
    IrrepInfo::IrrepInfo(Eigen::MatrixXcd _trans_mat, Eigen::VectorXcd _characters) :
      index(0),
      trans_mat(std::move(_trans_mat)),
      characters(std::move(_characters)) {
      complex = !almost_zero(trans_mat.imag());
    }


  }

  //*******************************************************************************************

  Index num_blocks(SymGroupRep const &_rep, SymGroup const &head_group) {
    Eigen::MatrixXd bmat(Local::_block_shape_matrix(_rep, head_group));

    if(bmat.cols() == 0)
      return 0;

    Index Nb = 1;

    //count zeros on first superdiagonal
    for(EigenIndex i = 0; i + 1 < bmat.rows(); i++) {
      if(almost_zero(bmat(i, i + 1)))
        Nb++;
    }

    return Nb;
  }

  //*******************************************************************************************

  //assumes that representation is real-valued and irreducible
  SymRepTools::Symmetrizer irrep_symmetrizer(SymGroupRep const &_rep,
                                             SymGroup const &head_group,
                                             double vec_compare_tol) {
    Index dim = _rep.MatrixXd(0)->cols();
    return irrep_symmetrizer(_rep, head_group, Eigen::MatrixXcd::Identity(dim, dim), vec_compare_tol);
  }

  //*******************************************************************************************

  //assumes that representation is real-valued and irreducible
  SymRepTools::Symmetrizer irrep_symmetrizer(SymGroupRep const &_rep,
                                             SymGroup const &head_group,
                                             Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                                             double vec_compare_tol) {

    // Result of sdirs is ordered by symmetry breaking properties of different directions.
    // We won't do any element-based comparisons after this point, in order to ensure that
    // our choice of axes is totally determined by symmetry
    multivector<Eigen::VectorXcd>::X<2>  sdirs = special_irrep_directions(_rep, head_group, _subspace, vec_compare_tol);
    return std::make_pair(Local::_irrep_symmetrizer_from_directions(sdirs, _subspace, vec_compare_tol), std::move(sdirs));
  }

  //*******************************************************************************************

  Eigen::MatrixXd full_trans_mat(std::vector<SymRepTools::IrrepInfo> const &irreps) {
    Index row = 0;
    Index col = 0;
    for(auto const &irrep : irreps) {
      col = irrep.vector_dim();
      row += irrep.irrep_dim();
    }
    Eigen::MatrixXd trans_mat(row, col);
    row = 0;
    for(auto const &irrep : irreps) {
      trans_mat.block(row, 0, irrep.irrep_dim(), irrep.vector_dim()) = irrep.trans_mat.real();
      row += irrep.irrep_dim();
    }
    return trans_mat;
  }

  //*******************************************************************************************

  multivector< Eigen::VectorXd >::X<3>
  special_total_directions(SymGroupRep const &_rep,
                           SymGroup const &head_group) {

    std::vector<SymRepTools::IrrepInfo> irreps = irrep_decomposition(_rep, head_group, false);

    multivector< Eigen::VectorXd >::X<3> result;
    result.reserve(irreps.size());

    for(auto const &irrep : irreps) {
      result.emplace_back(Local::_real(irrep.directions));
    }
    return result;
  }

  //*******************************************************************************************
  /// Returns array of orbits of high-symmetry directions in the vector space on
  /// which this representation is defined. This routine is different from special_total_directions
  /// in that it does not rely on the character tables to generate the special directions. The
  /// routine also adopts the faster method to generating the irreducible transformation matrix

  multivector< Eigen::VectorXcd >::X<2>
  special_irrep_directions(SymGroupRep const &_rep,
                           SymGroup const &head_group,
                           double vec_compare_tol) {
    Index dim = _rep.dim();
    return special_irrep_directions(_rep, head_group, Eigen::MatrixXcd::Identity(dim, dim), vec_compare_tol);
  }

  //*******************************************************************************************
  /// Returns array of orbits of high-symmetry directions in the vector space on
  /// which this representation is defined. This routine is different from special_total_directions
  /// in that it does not rely on the character tables to generate the special directions. The
  /// routine also adopts the faster method to generating the irreducible transformation matrix

  multivector< Eigen::VectorXcd >::X<2>
  special_irrep_directions(SymGroupRep const &_rep,
                           SymGroup const &head_group,
                           Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                           double vec_compare_tol,
                           bool all_subgroups) {
    auto sgroups = (all_subgroups ? head_group.subgroups() : head_group.small_subgroups());

    std::vector< Eigen::VectorXcd > tdirs;
    Eigen::MatrixXd R;
    Index dim = _rep.dim();

    // Loop over small (i.e., cyclic) subgroups and hope that each special
    // direction is invariant to at least one small subgroup
    for(auto const &orbit : sgroups) {
      // Reynolds for small subgroup *(orbit.begin()) in irrep subspace i
      R.setZero(dim, dim);

      for(Index op : * (orbit.begin())) {
        R += *(_rep.MatrixXd(head_group[op]));
      }

      if((R * _subspace).norm() < TOL)
        continue;

      // Find spanning vectors of column space of R*_subspace, which is projection of _subspace into its invariant component
      auto QR = (R * _subspace).colPivHouseholderQr();
      QR.setThreshold(TOL);
      //std::cout<< " ---------------------------------------------------------- " << std::endl;
      //std::cout << "Identity projector: \n" << R << "\n\n"
      //          << "Projected coordinates:\n" << R*_subspace << "\n\n";
      //std::cout << "Q matrix: \n" << Eigen::MatrixXd(QR.matrixQ()) << "\n";
      //std::cout << "R matrix: \n" << Eigen::MatrixXd(QR.matrixR().template triangularView<Eigen::Upper>()) << "\n";
      // // If only one spanning vector, it is special direction
      if(QR.rank() > 1)
        continue;
      Eigen::MatrixXcd Q = QR.matrixQ();
      //std::cout << "Special direction identified:\n" << Q.col(0).transpose() << "\n\n";
      //std::cout << "MatrixR:a\n" << Eigen::MatrixXd(QR.matrixR().template triangularView<Eigen::Upper>()) << "\n";
      // // PRINT OUT  the Reynolds matrix if its going to be added
      //std::cout<< " ---------------------------------------------------------- " << std::endl;
      //std::cout<< R <<std::endl;
      // Convert from irrep subspace back to total space and push_back
      tdirs.push_back(Q.col(0));
      tdirs.push_back(-Q.col(0));
      //std::cout<< " The added vector is : "<<t_result.back().back().transpose()<<std::endl;
      // //Calculate the rank using an SVD decomposition
      // Eigen::JacobiSVD<Eigen::MatrixXd> svd(R);
      // svd.setThreshold(1e-5);
      //std::cout<<"The SVD rank is          : "<<svd.rank()<<std::endl;
      //std::cout<<"The singular values are  : "<<svd.singularValues().transpose()<<std::endl;
      // //Loop over the SVD and count the number of non-zero singular values
      // int svd_rank = 0;
      // for(auto singular_index = 0 ; singular_index<svd.singularValues().size() ; singular_index++)
      //   svd_rank += (std::abs(svd.singularValues()[singular_index]) > 1e-5);
      //std::cout<<"The calculated rank is   : "<<svd_rank<<std::endl;
      //std::cout<<"The threshold is         : "<<svd.threshold()<<std::endl;
      //std::cout<<"The maxPivot  is         : "<<QR.maxPivot()<<std::endl;
      //std::cout<< " ---------------------------------------------------------- " << std::endl;
    }

    // t_result may contain duplicates, or elements that are equivalent by
    // symmetry. To discern more info, we need to exclude duplicates and find
    // the orbit of the directions. this should also
    // reveal the invariant subgroups.

    using SymCompareType = DirectionSymCompare<Eigen::VectorXcd, EigenSymRepApply<Eigen::VectorXcd>>;
    using VectorOrbit = Orbit<SymCompareType>;
    std::set<VectorOrbit> orbit_result;
    make_orbits(tdirs.begin(),
                tdirs.end(),
                head_group,
                _SymCompareType(vec_compare_tol, _rep),
                std::inserter(orbit_result, orbit_result.begin()));

    multivector<Eigen::VectorXcd>::X<2> result;
    //std::cout << "subspace: \n" << _subspace << std::endl;
    //std::cout << "Found " << orbit_result.size() << " direction orbits:\n";
    for(VectorOrbit const &orbit : orbit_result) {
      //std::cout << orbit.begin()->transpose() << std::endl;
      //std::cout << "---------\n";
      result.emplace_back(orbit.begin(), orbit.end());
    }

    if(all_subgroups || result.size() >= _subspace.cols()) {
      //std::cout << "special_irrep_direcions RETURNING RESULT\n";
      return result;
    }
    else {
      //std::cout << "result size: " << result.size() << "; _subspace cols: " << _subspace.cols() << "\n";
      //std::cout << "special_irrep_direcions REDOING CALCULATION\n";
      return special_irrep_directions(_rep, head_group, _subspace, vec_compare_tol, true);
    }
  }


  //*******************************************************************************************

  // Similar to algorithm for finding special directions, except we also find special planes,
  // n-planes, etc. The matrices that are returned have column vectors that span the special
  // subspaces.  The subspaces are arranged in orbits of subspaces that are equivalent by symmetry

  std::vector<std::vector< Eigen::MatrixXd> > special_subspaces(SymGroupRep const &_rep, SymGroup const &head_group) {
    if(!_rep.size() || !_rep.MatrixXd(0)) {
      default_err_log() << "CRITICAL ERROR: In special_subspaces() called on imcompatible SymGroupRep.\n Exiting...\n";
      exit(1);
    }

    Index i, j, k;

    std::vector<std::vector<Eigen::MatrixXd> > ssubs; //std::vector of orbits of special directions

    Eigen::MatrixXd tsub, ttrans;

    // Operation indices of subgroups
    auto const &sg_op_inds(head_group.subgroups());

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(*(_rep.MatrixXd(0)));

    QR.setThreshold(TOL);

    Eigen::MatrixXd Reynolds(*(_rep.MatrixXd(0)));

    // 'i' loops over subgroup "orbits". There are 0 to 'dim' special directions associated with each orbit.
    for(i = 0; i < sg_op_inds.size(); i++) {
      if(sg_op_inds[i].size() == 0 || sg_op_inds[i].begin()->size() == 0) {
        default_err_log() << "CRITICAL ERROR: In special_subspaces(), attempting to use a zero-size subgroup.\n Exiting...\n";
        exit(1);
      }

      Reynolds.setZero();
      // j loops over elements of "prototype" subgroup to get the Reynold's operator for the vector
      // space on which the representation is defined
      for(Index op : * (sg_op_inds[i].begin())) {
        Reynolds += *(_rep.MatrixXd(head_group[op]));
      }

      Reynolds /= double(sg_op_inds[i].begin()->size());

      // Column space of Reynold's matrix is invariant under the subgroup
      QR.compute(Reynolds);

      // If there are no special subspaces, the Reynold's operator has zero
      if(QR.rank() < 1)
        continue;

      // The first (QR.rank()) columns of QR.matrixQ() span a special subspace --
      // this can be accessed as QR.matrixQ().leftCols(QR.rank())
      tsub = Eigen::MatrixXd(QR.householderQ()).leftCols(QR.rank());

      // See if tsub has been found. Because we only keep orthonormal spanning vectors, this check is significantly simplified
      //  -- the matrix (tsub.transpose()*ssubs[j][k]) must have a determinant of +/-1 if the two matrices span the same column space
      for(j = 0; j < ssubs.size(); j++) {

        if(ssubs[j][0].cols() != tsub.cols())
          continue;

        for(k = 0; k < ssubs[j].size(); k++) {
          double tdet = (tsub.transpose() * ssubs[j][k]).determinant();
          if(almost_equal(std::abs(tdet), 1.0))
            break;
        }
        if(k < ssubs[j].size())
          break;
      }
      if(j < ssubs.size())
        continue;

      //tsub is new -- get equivalents
      ssubs.push_back(std::vector<Eigen::MatrixXd>());
      for(j = 0; j < head_group.size(); j++) {
        ttrans = (*(_rep.MatrixXd(head_group[j]))) * tsub;
        for(k = 0; k < ssubs.back().size(); k++) {
          double tdet = (ttrans.transpose() * ssubs.back()[k]).determinant();
          if(almost_equal(std::abs(tdet), 1.0))
            break;
        }
        if(k == ssubs.back().size())
          ssubs.back().push_back(ttrans);
      }

      // Unlike directions, the negative orientation of a subspace matrix is not considered
      // distinct. This is because we are only concerned about coordinate axes.

    }

    // Order subspaces by dimension, small to large; within blocks of equal dimension, order from
    // lowest multiplicity to highest
    for(i = 0; i < ssubs.size(); i++) {
      for(j = i + 1; j < ssubs.size(); j++) {
        if(ssubs[i][0].cols() > ssubs[j][0].cols())
          ssubs[i].swap(ssubs[j]);
        if(ssubs[i][0].cols() == ssubs[j][0].cols() && ssubs[i].size() > ssubs[j].size())
          ssubs[i].swap(ssubs[j]);
      }
    }

    return ssubs;
  }
  //*******************************************************************************************

  bool is_irrep(SymGroupRep const &_rep, SymGroup const &head_group) {
    double tvalue = 0;

    for(Index i = 0; i < head_group.size(); i++) {
      tvalue += (_rep[head_group[i].index()]->character()) * (_rep[head_group[i].index()]->character());
    }

    if(Index(round(tvalue)) == head_group.size()) {
      return true;
    }
    else {
      return false;
    }
  }

  //*******************************************************************************************

  VectorSpaceSymReport vector_space_sym_report(SymGroupRep const &_rep,
                                               SymGroup const &head_group,
                                               Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                                               bool calc_wedge) {

    VectorSpaceSymReport result;


    if(calc_wedge) {
      result.irreducible_wedge = SymRepTools::symrep_subwedges(_rep, head_group, _subspace);

      if(!result.irreducible_wedge.empty()) {
        result.irreps.reserve(result.irreducible_wedge[0].irrep_wedges().size());
        for(auto const &wedge : result.irreducible_wedge[0].irrep_wedges()) {
          result.irreps.push_back(wedge.irrep_info);
        }
      }
    }
    else {
      result.irreps = irrep_decomposition(_rep, head_group, _subspace, false);
    }
    result.symmetry_adapted_dof_subspace = full_trans_mat(result.irreps).adjoint();

    result.axis_glossary = std::vector<std::string>(result.symmetry_adapted_dof_subspace.rows(), "x");
    Index i = 0;
    for(std::string &x : result.axis_glossary) {
      x += std::to_string(++i);
    }

    for(SymOp const &op : head_group) {
      result.symgroup_rep.push_back(*(_rep.MatrixXd(op)));
    }
    return result;
  }



  //*******************************************************************************************

  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          bool allow_complex) {
    return irrep_decomposition(_rep,
                               head_group,
                               Eigen::MatrixXd::Identity(_rep.dim(), _rep.dim()),
                               allow_complex);
  }

  //*******************************************************************************************

  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                                                          bool allow_complex) {
    return irrep_decomposition(_rep,
                               head_group,
    [&](Eigen::Ref<const Eigen::MatrixXcd> const & f_subspace) {
      return irrep_symmetrizer(_rep, head_group, f_subspace, TOL);
    },
    _subspace,
    allow_complex);
  }

  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          SymRepTools::SymmetrizerFunction const &symmetrizer_func,
                                                          Eigen::MatrixXd subspace,
                                                          bool allow_complex) {

    if(_rep.dim() != subspace.rows()) {
      throw std::runtime_error("In irrep_decomposition, subspace matrix does not have proper number of rows (should have "
                               + std::to_string(_rep.dim()) + ", but has " + std::to_string(subspace.rows()));
    }


    Eigen::MatrixXd symspace(subspace.rows(), subspace.cols()*head_group.size());
    Index l = 0;
    if(!subspace.isIdentity()) {
      //Expand subpsace if it isn't already an invariant subspace
      for(auto const &op : head_group) {
        symspace.block(0, l, subspace.rows(), subspace.cols()) = (*(_rep[op.index()]->MatrixXd())) * subspace;
        l += subspace.cols();
      }
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(symspace);
      qr.setThreshold(TOL);
      subspace = Eigen::MatrixXd(qr.householderQ()).leftCols(qr.rank());
    }

    //std::cout << "subspace after expansion: \n" << subspace << std::endl;

    SymGroupRep sub_rep = coord_transformed_copy(_rep, subspace.transpose());

    auto subspace_symmetrizer = [&](const Eigen::Ref<const Eigen::MatrixXcd>& _subspace){
        return irrep_symmetrizer(sub_rep, head_group, _subspace, TOL);
    };

    std::vector<SymRepTools::IrrepInfo> irreps = irrep_decomposition(sub_rep, head_group, subspace_symmetrizer, allow_complex);

    std::vector<SymRepTools::IrrepInfo> result;
    result.reserve(irreps.size());
    l = 0;
    for(auto const &irrep : irreps) {
      //std::cout << "irrep " << ++l << " index: " << irrep.index << "\n";
      

      result.push_back(Local::_subspace_to_full_space(irrep, subspace));
    }
    return result;
  }
  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  std::vector<SymRepTools::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                          SymGroup const &head_group,
                                                          SymRepTools::SymmetrizerFunction const &symmetrizer_func,
                                                          bool allow_complex) {
    std::set<SymRepTools::IrrepInfo, Local::IrrepCompare > result{Local::IrrepCompare()};
    if(!_rep.size() || !head_group.size() || !_rep.MatrixXd(head_group[0])) {
      default_err_log() << "WARNING:  In irrep_decomposition(), size of representation is " << _rep.size() << " and MatrixXd address is " << _rep.MatrixXd(head_group[0]) << std::endl
                        << "          No valid irreps will be returned.\n";
      return {};
    }
    assert(Local::_rep_check(_rep, head_group, true) && "REPRESENTATION IS ILL-DEFINED!!");
    int dim(_rep.dim());
    std::vector<Eigen::MatrixXcd> commuters;
    std::cout.precision(8);
    std::cout.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    if(!is_irrep(_rep, head_group)) {
      // Identity always commutes, and by putting it first we prevent some accidental degeneracies
      commuters.push_back(Eigen::MatrixXcd::Identity(dim, dim) / sqrt(double(dim)));
    }
    typedef std::complex<double> cplx;
    std::vector<cplx> phase;
    phase.push_back(cplx(1.0, 0.0)); // 1+0i
    phase.push_back(cplx(0.0, 1.0)); // 0+1i

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esolve;
    Eigen::HouseholderQR<Eigen::MatrixXcd> qr;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> colqr;
    colqr.setThreshold(0.001);
    int Nfound(0);
    Eigen::MatrixXcd tmat(Eigen::MatrixXcd::Zero(dim, dim)), tcommute(Eigen::MatrixXcd::Zero(dim, dim));
    Eigen::MatrixXcd adapted_subspace(Eigen::MatrixXcd::Zero(dim, dim));

    // Initialize kernel as a random orthogonal matrix
    Eigen::MatrixXcd kernel(((1000 * Eigen::MatrixXcd::Random(dim, dim)) / 1000.).householderQr().householderQ());
    kernel.setIdentity(dim, dim);
    //std::cout << "kernel matrix:\n" << kernel << std::endl;
    // Build 'commuters', which span space of real matrices that commute with SymOpReps
    // 'kernel' is the kernel of adapted_subspace, and as the loop progresses, 'kernel' shrinks and the rank of adapted_subspace increases
    for(Index nph = 0; nph < 2; nph++) {
      for(Index kci = 0; kci < kernel.cols(); kci++) {
        bool found_new_irreps(false);
        for(Index kcj = kernel.cols() - 1; kci <= kcj;  kcj--) {
          //std::cout << "nph: " << nph << "  kci: " << kci << "  kcj: " << kcj << std::endl;
          if(kci == kcj && nph > 0) continue;
          tcommute.setZero();
          //std::cout << "kernel imag:\n" << kernel.imag() << std::endl;
          // form commuters by taking outer product of i'th column of the kernel with the j'th column
          // commuters are constructed to be self-adjoint, which assures eigenvalues are real
          tmat = phase[nph] * kernel.col(kci) * kernel.col(kcj).adjoint() // outer product
                 + std::conj(phase[nph]) * kernel.col(kcj) * kernel.col(kci).adjoint(); // adjoint of outer product
          //std::cout << "tmat:\n" << tmat << std::endl;
          //apply reynolds operator

          for(SymOp const &op : head_group) {
            //std::cout << "Op " << ns++ << ":\n" << (*(_rep.MatrixXd(op))) << std::endl;
            tcommute += (*(_rep.MatrixXd(op))) * tmat * (*(_rep.MatrixXd(op))).transpose();
          }


          //std::cout << "Raw commuter: \n" << tcommute << std::endl;

          //Do Gram-Shmidt while building 'commuters'

          for(Index nc = 0; nc < commuters.size(); nc++) {
            cplx tproj((commuters[nc].array()*tcommute.array().conjugate()).sum()); //Frobenius product
            //std::cout << "tproj" << tproj << std::endl;
            tcommute -= tproj * commuters[nc];
          }

          double tnorm((tcommute.array().conjugate()*tcommute.array()).sum().real());

          //std::cout << "Commuter: \n" << tcommute << std::endl;
          //std::cout << "norm: " << tnorm << std::endl;
          if(tnorm > TOL) {
            commuters.push_back(tcommute / sqrt(tnorm));
          }
          else continue;  // Attempt to construct the next commuter...


          //Finished building commuter now we can try to harvest irreps from it

          // construct adapted_subspace from the non-degenerate irreps obtained from the commuter

          Index nc = commuters.size() - 1;

          // magnify the range of eigenvalues to be (I think) independent of matrix dimension by multiplying by dim^{3/2}
          esolve.compute(double(dim)*sqrt(double(dim))*kernel.adjoint()*commuters[nc]*kernel);

          std::vector<Index> subspace_dims = partition_distinct_values(esolve.eigenvalues());

          //std::cout << "My big matrix: \n" << commuters[nc] << std::endl;
          //std::cout << "My little matrix: \n" << double(dim)*sqrt(double(dim))*kernel.adjoint()*commuters[nc]*kernel << std::endl;

          // Columns of tmat are orthonormal eigenvectors of commuter in terms of natural basis
          // (they were calculated in terms of kernel as basis)
          //tmat = kernel/*.conjugate()*/ * (esolve.eigenvectors().householderQr().householderQ());

          tmat = kernel * esolve.eigenvectors();

          //std::cout << "Raw eigenvecs: \n" << esolve.eigenvectors() << std::endl;
          //std::cout << "Rotated eigenvecs: \n" << tmat << std::endl;

          //std::cout << "Eigenvec unitarity: \n" << esolve.eigenvectors().adjoint()* esolve.eigenvectors() << std::endl;

          //std::cout << "rebuild little 1:\n"
          //<< esolve.eigenvectors()*esolve.eigenvalues().asDiagonal()*esolve.eigenvectors().inverse()
          //<< std::endl;

          //std::cout << "rebuild big 1:\n"
          //<< kernel*esolve.eigenvectors()*esolve.eigenvalues().asDiagonal()*esolve.eigenvectors().inverse()*kernel.adjoint()
          //<< std::endl;

          //std::cout << "rebuild big 2:\n"
          //<< tmat*esolve.eigenvalues().asDiagonal()*tmat.adjoint()
          //<< std::endl;



          // make transformed copy of the representation
          std::vector<Eigen::MatrixXcd> trans_rep(head_group.size());

          for(Index i = 0; i < head_group.size(); i++) {
            trans_rep[i] = tmat.adjoint() * (*_rep.MatrixXd(head_group[i])) * tmat;
          }
          //std::cout << "Eigenvals: " << esolve.eigenvalues().transpose() << " \ndims: " << subspace_dims << std::endl;

          Index last_i = 0;
          for(Index ns = 0; ns < subspace_dims.size(); ns++) {
            double sqnorm(0);
            Eigen::VectorXcd char_vec(head_group.size());

            // 'ns' indexes an invariant subspace
            // Loop over group operation and calculate character vector for the group representation on this
            // invariant subspace. If the squared norm of the character vector is equal to the group order,
            // the invariant subspace is also irreducible
            for(Index ng = 0; ng < trans_rep.size(); ng++) {
              char_vec[ng] = cplx(0, 0);
              for(Index i = last_i; i < last_i + subspace_dims[ns]; i++)
                char_vec[ng] += trans_rep[ng](i, i);
              //std::norm is squared norm
              sqnorm += std::norm(char_vec[ng]);
            }
            //std::cout << "char_vec:\n" << char_vec << std::endl;
            //std::cout << "subspace: " << ns << "; sqnorm: " << sqnorm << "\n";
            if(almost_equal(sqnorm, double(head_group.size()))) { // this representation is irreducible
              Eigen::MatrixXcd t_irrep_subspace;

              if(allow_complex) {
                t_irrep_subspace = tmat.block(0, last_i, dim, subspace_dims[ns]);
              }
              else {
                t_irrep_subspace.resize(dim, 2 * subspace_dims[ns]);

                t_irrep_subspace.leftCols(subspace_dims[ns])
                  = sqrt(2.0) * tmat.block(0, last_i, dim, subspace_dims[ns]).real().cast<std::complex<double> >();
                t_irrep_subspace.rightCols(subspace_dims[ns])
                  = sqrt(2.0) * tmat.block(0, last_i, dim, subspace_dims[ns]).imag().cast<std::complex<double> >();
              }
              //Only append to adapted_subspace if the new columns extend the space (i.e., are orthogonal)

              if(almost_zero((t_irrep_subspace.adjoint()*adapted_subspace).norm(), 0.001)) {
                qr.compute(t_irrep_subspace);
                //it seems stupid to use two different decompositions that do almost the same thing, but
                // HouseholderQR is not rank-revealing, and colPivHouseholder mixes up the columns of the Q matrix.
                colqr.compute(t_irrep_subspace);
                Index rnk = colqr.rank();
                //std::cout << "t_irrep_subspace with rank: " << rnk << " --\n" << t_irrep_subspace << std::endl << std::endl;
                t_irrep_subspace = Eigen::MatrixXcd(qr.householderQ()).leftCols(rnk);
                //std::cout << "MatrixR: \n" << colqr.matrixR() << std::endl;
                auto symmetrizer = symmetrizer_func(t_irrep_subspace);
                SymRepTools::IrrepInfo t_irrep((t_irrep_subspace * symmetrizer.first).adjoint(), char_vec);
                t_irrep.directions = Local::_real(symmetrizer.second);

                if(rnk == 2 * subspace_dims[ns])
                  t_irrep.pseudo_irrep = true;
                else
                  t_irrep.pseudo_irrep = false;

                // extend adapted_subspace (used to simplify next iteration)
                adapted_subspace.block(0, Nfound, dim, rnk) = t_irrep_subspace;//t_irrep.trans_mat.adjoint();


                auto it = result.find(t_irrep);
                if(it != result.end()) {
                  t_irrep.index++;
                  ++it;
                }

                while(it != result.end() && it->index > 0) {
                  t_irrep.index++;
                  ++it;
                }

                result.emplace_hint(it, std::move(t_irrep));


                Nfound += rnk;
                found_new_irreps = true;
              }
            }
            last_i += subspace_dims[ns];
          }
          if(found_new_irreps) {
            qr.compute(adapted_subspace);
            kernel = Eigen::MatrixXcd(qr.householderQ()).rightCols(dim - Nfound);
            kci = 0;
            kci--;
            break;
          }
        }
      }
    }

    return std::vector<SymRepTools::IrrepInfo>(std::make_move_iterator(result.begin()),
                                               std::make_move_iterator(result.end()));
  }

  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd irrep_trans_mat(SymGroupRep const &_rep, SymGroup const &head_group) {
    return full_trans_mat(irrep_decomposition(_rep,
                                              head_group,
    [&](Eigen::Ref<const Eigen::MatrixXcd> const & _subspace) {
      return irrep_symmetrizer(_rep, head_group, _subspace, TOL); //.transpose();
    },
    false));
  }

  //*******************************************************************************************

  SymGroupRep subset_permutation_rep(SymGroupRep const &permute_rep, const std::vector<std::set<Index>> &subsets) {
    SymGroupRep new_rep(SymGroupRep::NO_HOME, permute_rep.size());
    if(permute_rep.has_valid_master())
      new_rep = SymGroupRep(permute_rep.master_group());

    std::vector<Index> tperm(subsets.size());
    for(Index np = 0; np < permute_rep.size(); ++np) {
      Permutation const *perm_ptr(permute_rep.permutation(np));
      if(perm_ptr == NULL) {
        // This check may not make sense. Some representations may not have all operations.  If it's causing problems, comment it out
        default_err_log() << "CRITICAL ERROR: In subset_permutation_rep(SymGroupRep,std::vector<Index>::X2), permute_rep is missing at least one permutation!\n"
                          << "                Exiting...\n";
        assert(0);
        exit(1);
      }
      Permutation iperm = perm_ptr->inverse();
      for(Index ns = 0; ns < subsets.size(); ++ns) {
        std::set<Index> perm_sub;

        for(Index i : subsets[ns]) {
          perm_sub.insert(iperm[i]);
        }

        Index ns2;
        for(ns2 = 0; ns2 < subsets.size(); ++ns2) {
          if(subsets[ns] == perm_sub) {
            tperm[ns] = ns2;
            break;
          }
        }

        if(ns2 >= subsets.size()) {
          default_err_log() << "CRITICAL ERROR: In subset_permutation_rep(SymGroupRep,std::vector<Index>::X2), subsets are not closed under permute_rep permutations!\n"
                            << "                Exiting...\n";
          assert(0);
          exit(1);
        }
      }
      new_rep.set_rep(np, SymPermutation(tperm.begin(), tperm.end()));
    }
    return new_rep;
  }

  //*******************************************************************************************

  SymGroupRep permuted_direct_sum_rep(SymGroupRep const &permute_rep, const std::vector<SymGroupRep const *> &sum_reps) {
    SymGroupRep new_rep(SymGroupRep::NO_HOME, permute_rep.size());
    if(permute_rep.has_valid_master())
      new_rep = SymGroupRep(permute_rep.master_group());

    std::vector<Index> rep_dims(sum_reps.size(), 0);
    if(permute_rep.permutation(0) == NULL) {
      default_err_log() << "CRITICAL ERROR: In permuted_direct_sum_rep(SymGroupRep,std::vector<SymGroupRep const*>), permute_rep does not contain permutations!\n"
                        << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    for(Index i = 0; i < sum_reps.size(); i++) {
      if(sum_reps[i] == NULL) {
        continue;
      }
      if((permute_rep.has_valid_master() != sum_reps[i]->has_valid_master())
         || (permute_rep.has_valid_master() && &(permute_rep.master_group()) != &(sum_reps[i]->master_group()))
         || sum_reps[i]->MatrixXd(0) == NULL) {
        default_err_log() << "CRITICAL ERROR: In permuted_direct_sum_rep(SymGroupRep,std::vector<SymGroupRep const*>), found incompatible SymGroupReps!\n"
                          << "                Exiting...\n";
        assert(0);
        exit(1);
      }
      rep_dims[i] = sum_reps[i]->dim();
    }

    std::vector<Index> sum_inds(1, 0);
    std::partial_sum(rep_dims.begin(), rep_dims.end(), std::back_inserter(sum_inds));
    auto tot = std::accumulate(rep_dims.begin(), rep_dims.end(), 0);
    Eigen::MatrixXd sum_mat(tot, tot);
    Eigen::MatrixXd const *rep_mat_ptr(NULL);
    Permutation const *perm_ptr;
    for(Index g = 0; g < permute_rep.size(); g++) {
      sum_mat.setZero();
      perm_ptr = permute_rep.permutation(g);
      for(Index r = 0; r < perm_ptr->size(); r++) {
        if(rep_dims[r] == 0) continue;
        Index new_r = (*perm_ptr)[r];
        rep_mat_ptr = sum_reps[new_r]->MatrixXd(g);
        sum_mat.block(sum_inds[r], sum_inds[new_r], rep_dims[new_r], rep_dims[new_r]) = *rep_mat_ptr;
      }
      new_rep.set_rep(g, SymMatrixXd(sum_mat));
    }
    return new_rep;
  }

  //*******************************************************************************************

  SymGroupRep kron_rep(SymGroupRep const &LHS, SymGroupRep const &RHS) {
    if((LHS.has_valid_master() != RHS.has_valid_master()) || (LHS.has_valid_master() && &LHS.master_group() != &RHS.master_group())) {
      default_err_log() << "CRITICAL ERROR: In kron_rep(SymGroupRep,SymGroupRep), taking product of two incompatible SymGroupReps!\n"
                        << "                Exiting...\n";
      exit(1);
    }
    SymGroupRep new_rep(SymGroupRep::NO_HOME, LHS.size());
    if(LHS.has_valid_master())
      new_rep = SymGroupRep(LHS.master_group());

    Eigen::MatrixXd tmat;
    for(Index i = 0; i < LHS.size(); i++) {
      //std::cout << "rep1 matrix:\n";
      //std::cout << *(rep1->MatrixXd(i)) << '\n';
      //std::cout << "rep2 matrix:\n";
      //std::cout << *(rep2->MatrixXd(i)) << '\n';
      assert(LHS.MatrixXd(i) && RHS.MatrixXd(i));
      kroneckerProduct(*(LHS.MatrixXd(i)), *(RHS.MatrixXd(i)), tmat);
      //std::cout << "Total matrix:\n" << tmat << '\n';
      new_rep.set_rep(i, SymMatrixXd(tmat));
    }
    return new_rep;
  }

  namespace SymRepTools {

    IrrepWedge IrrepWedge::make_dummy_irrep_wedge(const Eigen::MatrixXd& axes){
        IrrepWedge irrep_wedge(IrrepInfo::make_dummy(axes), axes);
        irrep_wedge.mult.reserve(axes.cols());
        for (Index i=0; i < axes.cols(); ++i){
            irrep_wedge.mult.push_back(1);
        }
        return irrep_wedge;
    }

    Eigen::MatrixXd SubWedge::_subwedge_to_trans_mat(std::vector<IrrepWedge> const &_iwedges) {

      if(_iwedges.empty())
        return Eigen::MatrixXd();
      Eigen::MatrixXd result(_iwedges[0].axes.rows(), _iwedges[0].axes.rows());
      Index i = 0;
      for(Index w = 0; w < _iwedges.size(); ++w) {
        result.block(0, i, _iwedges[w].axes.rows(), _iwedges[w].axes.cols()) = _iwedges[w].axes;
        i += _iwedges[w].axes.cols();
      }
      if(i < result.cols())
        result.conservativeResize(Eigen::NoChange, i);
      return result;
    }

    //*******************************************************************************************

    std::vector<IrrepWedge> irrep_wedges(SymGroup const &head_group,
                                         SymGroupRepID id,
                                         Eigen::Ref<const Eigen::MatrixXd> const &_subspace) {
      return irrep_wedges(head_group.master_group().representation(id), head_group, _subspace);
    }

    //*******************************************************************************************

    std::vector<IrrepWedge> irrep_wedges(SymGroupRep const &_rep,
                                         SymGroup const &head_group,
                                         Eigen::Ref<const Eigen::MatrixXd> const &_subspace) {

      std::vector<SymRepTools::IrrepInfo> irreps = irrep_decomposition(_rep, head_group, _subspace, false);
      std::vector<IrrepWedge> wedges;
      wedges.reserve(irreps.size());
      double best_proj, tproj;

      //std::cout << "subspace: \n" << _subspace << std::endl;
      //std::cout << "irrep decomposition size: " << irreps.size() << std::endl;

      for(SymRepTools::IrrepInfo const &irrep : irreps) {
        //1D irreps directions can have positive and negative directions, but we only want to include one.
        //only 1D irreps can have singly-degenerate directions and two singly-degenerate directions indicate
        //the same vector duplicated in positive and negative direction (because they are not equivalent by symmetry)
        //If irrep.directions[0] is singly degenerate (orbits size == 1) then irrepdim is 1 and we only need one direction to
        //define wedge
        wedges.emplace_back(irrep, Eigen::MatrixXd::Zero(irrep.vector_dim(), irrep.irrep_dim()));
        //std::cout << "Irrep characters: \n" << irrep.characters << std::endl;
        //std::cout << "Irrep directions: " << irrep.directions.size() << std::endl;
        if(irrep.directions.empty()) {
          wedges.back() = Local::_wedge_from_pseudo_irrep(irrep, _rep, head_group);
          continue;
        }

        //std::cout << "Irrep direction orbit" << 0 << " : " << irrep.directions[0].size() << std::endl;
        //std::cout << "Irrep direction: " << irrep.directions[0][0].transpose() << std::endl;
        wedges.back().axes.col(0) = irrep.directions[0][0];
        wedges.back().mult.push_back(irrep.directions[0].size());
        for(Index i = 1; i < irrep.irrep_dim(); i++) {
          //std::cout << "Irrep direction orbit" << i << " : " << irrep.directions[i].size() << std::endl;
          //std::cout << "Irrep direction: " << irrep.directions[i][0].transpose() << std::endl;
          Index j_best = 0;
          best_proj = (wedges.back().axes.transpose() * irrep.directions[i][0]).sum();
          for(Index j = 1; j < irrep.directions[i].size(); j++) {
            tproj = (wedges.back().axes.transpose() * irrep.directions[i][j]).sum();
            if(tproj > best_proj) {
              best_proj = tproj;
              j_best = j;
            }
          }

          wedges.back().axes.col(i) = irrep.directions[i][j_best];
          wedges.back().mult.push_back(irrep.directions[i].size());
        }
        //std::cout << "New irrep wedge: \n" << wedges.back().axes.transpose() << std::endl;
      }
      return wedges;
    }

    //*******************************************************************************************

    std::vector<SubWedge> symrep_subwedges(SymGroup const &head_group, SymGroupRepID id) {
      return symrep_subwedges(head_group.master_group().representation(id), head_group);
    }

    //*******************************************************************************************

    std::vector<SubWedge> symrep_subwedges(SymGroup const &head_group,
                                           SymGroupRepID id,
                                           Eigen::Ref<const Eigen::MatrixXd> const &_subspace) {
      return symrep_subwedges(head_group.master_group().representation(id), head_group, _subspace);
    }

    //*******************************************************************************************

    std::vector<SubWedge> symrep_subwedges(SymGroupRep const &_rep,
                                           SymGroup const &head_group) {
      return symrep_subwedges(_rep,
                              head_group,
                              Eigen::MatrixXd::Identity(_rep.dim(), _rep.dim()));

    }
    //*******************************************************************************************

    std::vector<SubWedge> symrep_subwedges(SymGroupRep const &_rep,
                                           SymGroup const &head_group,
                                           Eigen::Ref<const Eigen::MatrixXd> const &_subspace) {

      auto irrep_wedge_compare = [](const IrrepWedge & a, const IrrepWedge & b)->bool {
        return Eigen::almost_equal(a.axes, b.axes);
      };

      auto tot_wedge_compare =
      [irrep_wedge_compare](const std::vector<IrrepWedge> &a, const std::vector<IrrepWedge> &b)->bool {
        if(a.size() != b.size())
          return false;
        for(auto ita = a.begin(), itb = b.begin(); ita != a.end(); ++ita, ++itb)
          if(!irrep_wedge_compare(*ita, *itb))
            return false;
        //std::cout << "Equal: \n" << SubWedge(a).trans_mat() << ",\n" << SubWedge(b).trans_mat() << "\n\n";
        return true;
      };

      if(!_rep.MatrixXd(0))
        throw std::runtime_error("In symrep_subwedges, SymGroupRep does not describe matrix representation");

      std::vector<IrrepWedge> init_wedges = irrep_wedges(_rep, head_group, _subspace);


      std::vector<SubWedge> result;

      //irrep_wedge_orbits[w] is orbit of wedges[w]
      multivector<IrrepWedge>::X<2> irrep_wedge_orbits;
      irrep_wedge_orbits.reserve(init_wedges.size());
      //max_equiv[w] is irrep_wedge_orbits[w].size()-1
      std::vector<Index> max_equiv;
      max_equiv.reserve(init_wedges.size());
      //std::cout << "irreducible wedges for group of order " << head_group.size() << std::endl;
      Index imax = 0;
      multivector<Index>::X<2> subgroups;
      for(IrrepWedge const &wedge : init_wedges) {
        //std::cout << "Working wedge with axes: \n" << wedge.axes.transpose() << std::endl;
        irrep_wedge_orbits.push_back({wedge});

        //Start getting orbit of wedges[w]
        subgroups.push_back({});
        for(SymOp const &op : head_group) {
          IrrepWedge test_wedge{wedge};
          test_wedge.axes = (*(_rep[op.index()]->MatrixXd())) * wedge.axes;
          Index o = 0;
          for(; o < irrep_wedge_orbits.back().size(); ++o) {
            if(irrep_wedge_compare(irrep_wedge_orbits.back()[o], test_wedge)) {
              if(o == 0) {
                subgroups.back().push_back(op.index());
              }
              break;
            }
          }
          if(o < irrep_wedge_orbits.back().size())
            continue;
          irrep_wedge_orbits.back().push_back(test_wedge);
        }
        //std::cout << "wedge mult: " << irrep_wedge_orbits.back().size();
        //std::cout << "; subgroups[" << subgroups.size() << "]: " << subgroups.back() << std::endl;
        //std::cout << "N equiv wedges found: " << irrep_wedge_orbits.back().size() << std::endl;
        max_equiv.push_back(irrep_wedge_orbits.back().size() - 1);
        if(max_equiv.back() > max_equiv[imax])
          imax = max_equiv.size() - 1;
      }
      max_equiv[imax] = 0;

      //Counter over combinations of equivalent wedges
      Counter<std::vector<Index> > wcount(std::vector<Index>(init_wedges.size(), 0),
                                          max_equiv,
                                          std::vector<Index>(init_wedges.size(), 1));

      //std::cout << "max_equiv: " << max_equiv << std::endl;
      //std::cout << "init wcount: " << wcount() << std::endl;
      multivector<IrrepWedge>::X<3> tot_wedge_orbits;
      //std::cout << "Starting slow bit!\n";
      for(; wcount; ++wcount) {
        std::vector<IrrepWedge> twedge = init_wedges;
        for(Index i = 0; i < init_wedges.size(); i++)
          twedge[i].axes = irrep_wedge_orbits[i][wcount[i]].axes;


        if(contains_if(tot_wedge_orbits,
        [tot_wedge_compare, &twedge](const multivector<IrrepWedge>::X<2> &wedge_orbit)->bool {
        return contains(wedge_orbit,
                        twedge,
                        tot_wedge_compare);
        }))
        continue;

        tot_wedge_orbits.push_back({twedge});
        result.emplace_back(twedge);
        for(Index p : subgroups[imax]) {
          for(Index i = 0; i < twedge.size(); i++)
            twedge[i].axes = (*(_rep[p]->MatrixXd())) * result.back().irrep_wedges()[i].axes;
          if(!contains(tot_wedge_orbits.back(), twedge, tot_wedge_compare)) {
            tot_wedge_orbits.back().push_back(twedge);
          }
        }
      }
      //std::cout << "Num subwedges: " << result.size() << std::endl;
      return result;
    }

  }
}
