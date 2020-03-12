#include <numeric>
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

namespace Local {

  struct IrrepCompare {
    bool operator()(IrrepInfo const &irrep_a, IrrepInfo const &irrep_b) const {
      Eigen::VectorXcd const &a = irrep_a.characters;
      Eigen::VectorXcd const &b = irrep_b.characters;

      typedef std::complex<double> C;

      // Check if a and b are identity irrep, identity always compares less than other irreps
      bool a_id = almost_equal(a[0], C(1., 0.)) && almost_equal(a.sum(), C(a.size(), 0.));
      bool b_id = almost_equal(b[0], C(1., 0.)) && almost_equal(b.sum(), C(b.size(), 0.));

      if(a_id != b_id)
        return a_id;
      else if(a_id) {
        // Index is used to break ties if they are both identity
        return irrep_a.index < irrep_b.index;
      }

      // Low-dimensional irreps come before higher dimensional
      if(!almost_equal(a[0], b[0]))
        return a[0] < b[0];

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
}

namespace CASM {

  Eigen::MatrixXd block_shape_matrix(SymGroupRep const &_rep, SymGroup const &head_group) {
    if(!_rep.size() || !head_group.size() || !_rep.MatrixXd(head_group[0]))
      return Eigen::MatrixXd();

    Eigen::MatrixXd block_shape(_rep.MatrixXd(head_group[0])->cwiseProduct(*_rep.MatrixXd(head_group[0])));

    for(Index i = 1; i < head_group.size(); i++) {
      block_shape += _rep.MatrixXd(head_group[i])->cwiseProduct(*_rep.MatrixXd(head_group[i]));
    }

    return block_shape;
  }

  //*******************************************************************************************

  Index num_blocks(SymGroupRep const &_rep, SymGroup const &head_group) {
    Eigen::MatrixXd bmat(block_shape_matrix(_rep, head_group));

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
  SymRepTools::Symmetrizer irrep_symmetrizer_and_directions(SymGroupRep const &_rep,
                                                            SymGroup const &head_group,
                                                            double vec_compare_tol) {
    Index dim = _rep.MatrixXd(0)->cols();
    return irrep_symmetrizer_and_directions(_rep, head_group, Eigen::MatrixXcd::Identity(dim, dim), vec_compare_tol);
  }

  //*******************************************************************************************

  //assumes that representation is real-valued and irreducible
  SymRepTools::Symmetrizer irrep_symmetrizer_and_directions(SymGroupRep const &_rep,
                                                            SymGroup const &head_group,
                                                            Eigen::Ref<const Eigen::MatrixXcd> const &_subspace,
                                                            double vec_compare_tol) {

    // Result of sdirs is ordered by symmetry breaking properties of different directions.
    // We won't do any element-based comparisons after this point, in order to ensure that
    // our choice of axes is totally determined by symmetry
    multivector<Eigen::VectorXcd>::X<2>  sdirs = special_irrep_directions(_rep, head_group, _subspace, vec_compare_tol);
    return std::make_pair(irrep_symmetrizer_from_directions(sdirs, _subspace, vec_compare_tol), std::move(sdirs));
  }
  //*******************************************************************************************

  //assumes that representation is real-valued and irreducible
  Eigen::MatrixXcd irrep_symmetrizer_from_directions(multivector<Eigen::VectorXcd>::X<2> const &special_directions,
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
      result = _subspace.colPivHouseholderQr().solve(representation_prepare_impl(obj, vec_compare_tol));
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

  Eigen::MatrixXd full_trans_mat(std::vector<IrrepInfo> const &irreps) {
    Index row = 0;
    Index col = 0;
    for(auto const &irrep : irreps) {
      col = irrep.vector_dim();
      row += irrep.irrep_dim();
    }
    Eigen::MatrixXd trans_mat(row, col);
    row = 0;
    for(auto const &irrep : irreps) {
      trans_mat.block(row, 0, irrep.irrep_dim(), irrep.vector_dim()) = irrep.trans_mat;
      row += irrep.irrep_dim();
    }
    return trans_mat;
  }

  //*******************************************************************************************

  multivector< Eigen::VectorXd >::X<3>
  special_total_directions(SymGroupRep const &_rep,
                           SymGroup const &head_group) {

    std::vector<IrrepInfo> irreps = irrep_decomposition(_rep, head_group);

    multivector< Eigen::VectorXd >::X<3> result;
    result.reserve(irreps.size());

    for(auto const &irrep : irreps) {
      result.push_back(irrep.directions);
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
    Index dim = _rep.dim()

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
    using VectorOrbit = GenericOrbit<SymCompareType>;
    std::set<VectorOrbit> orbit_result;
    make_orbits(tdirs.begin(),
                tdirs.end(),
                head_group,
                SymCompareType(vec_compare_tol, _rep),
                std::inserter(orbit_result, orbit_result.begin()));

    multivector<Eigen::VectorXd>::X<2> result;
    Index o(0);
    for(VectorOrbit const &orbit : orbit_result) {
      //std::cout << "Orbit " << ++o << ": \n";
      //for(auto const &v : orbit) {
      //std::cout << v.transpose() << "\n";
      //}
      //std::cout << "---------\n";
      result.emplace_back(orbit.begin(), orbit.end());
    }

    if(all_subgroups || !result.empty()) {
      //std::cout << "special_irrep_direcions RETURNING RESULT\n";
      return result;
    }
    else {
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
      tsub = QR.householderQ().eval().leftCols(QR.rank());

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
  //assumes that representation is real-valued, and combines complex-valued irreps with their complex conjugate
  std::vector<Index> num_each_real_irrep(SymGroupRep const &_rep, SymGroup const &head_group, bool verbose) {
    const std::vector<bool > &complex_irrep(head_group.get_complex_irrep_list());
    std::vector<Index> tarray(num_each_irrep(_rep, head_group, verbose));

    if(tarray.size() != complex_irrep.size()) {
      default_err_log() << "CRITICAL ERROR: Dimension mismatch in num_each_real_irrep. Exiting..\n";
      assert(0);
      exit(1);
    }

    // Go through and remove double counting of complex irreps
    for(Index i = 0; i < tarray.size(); i++) {
      if(tarray[i] && complex_irrep[i]) {
        if(i + 1 == tarray.size() || tarray[i] != tarray[i + 1]) {
          default_err_log() << "CRITICAL ERROR: Invalid irrep decomposition found in num_each_real_irrep. Exiting..\n";
          assert(0);
          exit(1);
        }
        tarray[i + 1] = 0;
        i++;
      }
    }
    return tarray;
  }

  //*******************************************************************************************

  std::vector<Index> num_each_irrep(SymGroupRep const &_rep, SymGroup const &head_group, bool verbose) {
    const std::vector<std::vector<Index> > &conj_class(head_group.get_conjugacy_classes());
    const std::vector<std::vector<std::complex<double> > > &char_table(head_group.character_table());

    std::vector<Index> tdec(conj_class.size(), 0);
    std::vector<double> repchar(conj_class.size(), 0);

    if(verbose) {
      std::cout << "Decomposing representation:\n ";
      for(Index i = 0; i < head_group.size(); i++) {
        std::cout << *(_rep.MatrixXd(head_group[i])) << "\n\n";
      }
    }

    //std::cout << "Character table is\n";
    for(Index i = 0; i < conj_class.size(); i++) {
      Index rep_index(head_group[conj_class[i][0]].index());
      repchar[i] = _rep[rep_index]->character();
      //std::cout << char_table[i] << "\n\n";
    }

    if(verbose)
      std::cout << " Irrep decomposition: ";
    for(Index i = 0; i < char_table.size(); i++) { // Loop over irreducible representations
      double temp(0);
      for(Index j = 0; j < char_table[i].size(); j++) { // Loop over conjugacy classes
        temp += double(conj_class[j].size()) * repchar[j] * (char_table[i][j]).real();
      }
      tdec[i] = round(temp / double(head_group.size()));
      if(verbose)
        std::cout << "  " << tdec[i];
    }
    if(verbose)
      std::cout << '\n';


    return tdec;
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
                                               bool calc_wedge = false) {

    VectorSpaceSymReport result;


    if(calc_wedge) {
      auto wedges_and_irreps = SymRepTools::symrep_subwedges(_rep, head_group, _subspace);
      result.irreducible_wedge = std::move(wedges_and_irreps.first);
      result.irreps = std::move(wedges_and_irreps.second);
    }
    else {
      result.irreps = irrep_decomposition(_rep, head_group, _subspace);
    }

    std::vector<SymRepTools::SubWedge> wedges;
    std::vector<int> dims;
    SymGroup pg = make_point_group(_config.group(), _config.supercell().sym_info().supercell_lattice());
    DoFKey strain_dof_key;
    std::vector<DoFKey> tdof_types = global_dof_types(_primclex.prim());
    Index istrain = find_index_if(tdof_types,
    [](DoFKey const & other) {
      return other.find("strain") != std::string::npos;
    });
    if(istrain == tdof_types.size())
      throw std::runtime_error("Cannot enumerate strains for project in which strain has not been specified as a degree of freedom.");
    strain_dof_key = tdof_types[istrain];
    if(!sym_axes)
      wedges.push_back(SymRepTools::SubWedge({SymRepTools::IrrepWedge(_axes, std::vector<Index>(_axes.cols(), 1))}));
    else
      wedges = SymRepTools::symrep_subwedges(pg, _primclex.prim().global_dof(strain_dof_key).symrep_ID());

    ///LOCAL

    std::pair<Eigen::MatrixXd, std::vector<Index>> normcoords = collective_dof_normal_coords_and_irrep_dims(config.sites().begin(),
                                                                config.sites().end(),
                                                                config.supercell().sym_info(),
                                                                _dof,
                                                                config.group(),
                                                                _axes);
    axes = normcoords.first.transpose();

    //std::cout << "Axes:\n" << axes.transpose().format(tformat) << "\n";
    if(axes.cols() != _axes.cols()) {
      throw std::runtime_error("In ConfigEnumSiteDoFs, symmetry-adapted axes do not have same dimension as provided axes. "
                               "Please ensure that provided axes completely span one or more of subspaces listed above.");
    }
  }



  //*******************************************************************************************

  std::vector<IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                             SymGroup const &head_group) {
    return irrep_decomposition(_rep,
                               head_group,
                               Eigen::MatrixXd::Identity(_rep.dim(), _rep.dim()));
  }

  //*******************************************************************************************

  std::vector<IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                             SymGroup const &head_group,
                                             Eigen::Ref<const Eigen::MatrixXd> const &_subspace) {
    return irrep_decomposition(_rep,
                               head_group,
    [&](Eigen::Ref<const Eigen::MatrixXd> const & f_subspace) {
      return irrep_symmetrizer_and_directions(_rep, head_group, f_subspace, TOL);
    },
    _subspace);
  }

  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  std::vector<IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                             SymGroup const &head_group,
                                             SymRepTools::SymmetrizerFunction const &symmetrizer_func,
                                             Eigen::Ref<const Eigen::MatrixXd> const &_subspace) {

    Eigen::MatriXd symspace(_subspace.rows(), _subspace.cols()*head_group.size());

    Index l = 0;
    if(_rep.dim() != _subspace.rows()) {
      throw std::runtime_error("In irrep_decomposition, subspace matrix does not have proper number of rows (should have "
                               + std::to_string(_rep.dim()) + ", but has " + std::to_string(_subspace.rows()));
    }

    Eigen::MatrixXd subspace;
    if(!_subspace.isIdentity()) {
      for(auto const &op : head_group) {
        symspace.block(0, l, _subspace.rows(), _subspace.cols()) = (*(_rep[op.index()]->MatrixXd())) * _subspace;
      }
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(symspace);
      qr.setThreshold(TOL);
      subspace = (qr.householderQ())eval().leftcols(qr.rank());
    }
    else {
      subspace = _subspace;
    }

    SymGroupRep sub_rep = _rep.coord_transformed_copy(subspace.transpose());

    auto irreps = irrep_decomposition(sub_rep, head_group, symmetrizer_func);
    for(auto &irrep : irreps) {
      irrep.trans_mat *= subspace.transpose();
    }
    return irreps;
  }
  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  std::vector<Local::IrrepInfo> irrep_decomposition(SymGroupRep const &_rep,
                                                    SymGroup const &head_group,
                                                    SymRepTools::SymmetrizerFunction const &symmetrizer_func,
                                                    bool allow_complex) {
    std::set<Local::IrrepInfo, Local::CharCompare > result(CharCompare());
    if(!_rep.size() || !head_group.size() || !_rep.MatrixXd(head_group[0])) {
      default_err_log() << "WARNING:  In irrep_decomposition(), size of representation is " << _rep.size() << " and MatrixXd address is " << _rep.MatrixXd(head_group[0]) << std::endl
                        << "          No valid irreps will be returned.\n";
      return std::make_pair(Eigen::MatrixXd(), char_table);
    }
    assert(Local::_rep_check(_rep, head_group, true) && "REPRESENTATION IS ILL-DEFINED!!");
    int dim(_rep.MatrixXd(head_group[0])->rows());
    std::vector<Eigen::MatrixXcd> commuters;
    std::cout.precision(8);
    std::cout.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    // Identity always commutes, and by putting it first we prevent some accidental degeneracies
    commuters.push_back(Eigen::MatrixXcd::Identity(dim, dim) / sqrt(double(dim)));
    typedef std::complex<double> cplx;
    std::vector<cplx> phase;
    phase.push_back(cplx(1.0, 0.0)); // 1+0i
    phase.push_back(cplx(0.0, 1.0)); // 0+1i

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esolve;
    Eigen::HouseholderQR<Eigen::MatrixXd> qr;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> colqr;
    colqr.setThreshold(0.001);
    int Nfound(0);
    Eigen::MatrixXcd tmat(Eigen::MatrixXcd::Zero(dim, dim)), tcommute(Eigen::MatrixXcd::Zero(dim, dim));
    Eigen::MatrixXd trans_mat(Eigen::MatrixXd::Zero(dim, dim));

    // Initialize kernel as a random orthogonal matrix
    Eigen::MatrixXd kernel(((1000 * iround(Eigen::MatrixXd::Random(dim, dim))).cast<double>() / 1000.).householderQr().householderQ());

    // Build 'commuters', which span space of real matrices that commute with SymOpReps
    // 'kernel' is the kernel of trans_mat, and as the loop progresses, 'kernel' shrinks and the rank of trans_mat increases
    for(Index kci = 0; kci < kernel.cols(); kci++) {
      bool found_new_irreps(false);
      for(Index kcj = kci; kcj < kernel.cols(); kcj++) {
        for(Index nph = 0; nph < 2; nph++) {
          if(kci == kcj && nph > 0) continue;
          tcommute.setZero();

          // form commuters by taking outer product of i'th column of the kernel with the j'th column
          // commuters are constructed to be self-adjoint, which assures eigenvalues are real
          tmat = phase[nph] * kernel.col(kci) * kernel.col(kcj).transpose() // outer product
                 + std::conj(phase[nph]) * kernel.col(kcj) * kernel.col(kci).transpose(); // adjoint of outer product

          //apply reynolds operator
          for(Index ns = 0; ns < head_group.size(); ns++) {
            tcommute += (*(_rep.MatrixXd(head_group[ns]))) * tmat * (*(_rep.MatrixXd(head_group[ns]))).transpose();
          }

          //Do Gram-Shmidt while building 'commuters'

          for(Index nc = 0; nc < commuters.size(); nc++) {
            cplx tproj((commuters[nc].array().conjugate()*tcommute.array()).sum()); //Frobenius product
            tcommute -= tproj * commuters[nc];
          }

          double tnorm((tcommute.array().conjugate()*tcommute.array()).sum().real());
          if(tnorm > TOL) {
            commuters.push_back(tcommute / tnorm);
          }
          else continue;  // Attempt to construct the next commuter...


          //Finished building commuter now we can try to harvest irreps from it

          // construct trans_mat from the non-degenerate irreps obtained from the commuter

          Index nc = commuters.size() - 1;

          // magnify the range of eigenvalues to be (I think) independent of matrix dimension by multiplying by dim^{3/2}
          esolve.compute(double(dim)*sqrt(double(dim))*kernel.transpose()*commuters[nc]*kernel);

          std::vector<Index> subspace_dims = partition_distinct_values(esolve.eigenvalues());

          // Columns of tmat are orthonormal eigenvectors of commuter in terms of natural basis
          // (they were calculated in terms of kernel as basis)
          tmat = kernel * (esolve.eigenvectors().householderQr().householderQ());


          // make transformed copy of the representation
          std::vector<Eigen::MatrixXcd> trans_rep(head_group.size());
          Eigen::MatrixXd block_shape(Eigen::MatrixXd::Zero(kernel.cols(), kernel.cols()));

          for(Index i = 0; i < head_group.size(); i++) {
            trans_rep[i] = tmat.adjoint() * (*_rep.MatrixXd(head_group[i])) * tmat;
            block_shape += (trans_rep[i].cwiseProduct(trans_rep[i].conjugate())).real();
          }

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

            if(almost_equal(sqnorm, double(head_group.size()))) { // this representation is irreducible
              Eigen::MatrixXd ttrans_mat(dim, 2 * subspace_dims[ns]);

              ttrans_mat.leftCols(subspace_dims[ns]) = sqrt(2.0) * tmat.block(0, last_i, dim, subspace_dims[ns]).real();
              ttrans_mat.rightCols(subspace_dims[ns]) = sqrt(2.0) * tmat.block(0, last_i, dim, subspace_dims[ns]).imag();

              //Only append to trans_mat if the new columns extend the space (i.e., are orthogonal)

              if(almost_zero((ttrans_mat.transpose()*trans_mat).norm(), 0.001)) {
                qr.compute(ttrans_mat);
                //it seems stupid to use two different decompositions that do almost the same thing, but
                // HouseholderQR is not rank-revealing, and colPivHouseholder mixes up the columns of the Q matrix.
                colqr.compute(ttrans_mat);
                Index rnk = colqr.rank();

                ttrans_mat = Eigen::MatrixXd(qr.householderQ()).leftCols(rnk);
                auto symmetrizer = symmetrizer_func(ttrans_mat);
                Local::IrrepInfo t_info(ttrans_mat * symmetrizer.first, char_vec);
                t_info.directions = symmetrizer.second;

                // extend trans_mat (used to simplify next iteration)
                trans_mat.block(0, Nfound, dim, rnk) = t_info.trans_mat;


                auto it = result.find(t_info);
                if(it != result.end()) {
                  t_info.index++;
                  ++it;
                }

                while(it != result.end() && it->index > 0) {
                  t_info.index++;
                  ++it;
                }
                result.emplace_hint(it, std::move(t_info));


                Nfound += rnk;
                found_new_irreps = true;
              }
            }
            last_i += subspace_dims[ns];
          }
          if(found_new_irreps) {
            qr.compute(trans_mat);
            kernel = Eigen::MatrixXd(qr.householderQ()).rightCols(dim - Nfound);
            break;
          }
        }
        if(found_new_irreps) {
          kci = 0;
          break;
        }
      }
    }


    return std::vector<IrrepInfo>(std::make_move_iterator(result.begin()),
                                  std::make_move_iterator(result.end()));

  }

  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd irrep_trans_mat(SymGroupRep const &_rep, SymGroup const &head_group) {
    return full_trans_mat(irrep_decomposition(_rep,
                                              head_group,
    [&](Eigen::Ref<const Eigen::MatrixXd> const & _subspace) {
      return irrep_symmetrizer_and_directions(_rep, head_group, _subspace, TOL); //.transpose();
    }));
  }

  //*******************************************************************************************

  std::vector<Eigen::MatrixXd> irrep_projection_operators(SymGroupRep const &_rep, SymGroup const &head_group) {

    const std::vector<std::vector<Index> > &conj_class(head_group.get_conjugacy_classes());
    const std::vector<std::vector<std::complex<double> > > &char_table(head_group.character_table());

    Eigen::MatrixXd tmat((*(_rep.MatrixXd(0))).rows(), (*(_rep.MatrixXd(0))).cols());
    tmat.setZero();
    std::vector<Eigen::MatrixXd> tarray(conj_class.size());

    for(Index i = 0; i < conj_class.size(); i++) {
      double dimension = char_table[i][0].real();
      for(Index j = 0; j < head_group.size(); j++) {
        double character =  char_table[i][head_group.class_of_op(j)].real();
        Eigen::MatrixXd rep_mat = *(_rep.MatrixXd(j));
        tmat += (character * rep_mat);
      }
      tarray[i] = ((dimension) / (head_group.size())) * tmat;
    }

    return tarray;
  }

  //*******************************************************************************************

  /// \brief Make copy of (*this) that is transformed so that axes are oriented along high-symmetry direction
  /// and confined to subspaces that transform as irreps.
  SymGroupRep symmetry_adapted_copy(SymGroupRep const &_rep, SymGroup const &head_group) {
    return coord_transformed_copy(_rep, irrep_trans_mat(_rep, head_group));
  }

  //*******************************************************************************************

  SymGroupRep subset_permutation_rep(SymGroupRep const &permute_rep, const std::vector<std::vector<Index>> &subsets) {
    SymGroupRep new_rep(SymGroupRep::NO_HOME, permute_rep.size());
    if(permute_rep.has_valid_master())
      new_rep = SymGroupRep(permute_rep.master_group());

    std::vector<Index> perm_sub;
    if(subsets.size())
      perm_sub.resize(subsets[0].size());

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
        for(Index ne = 0; ne < subsets[ns].size(); ++ne) {
          perm_sub[ne] = iperm[subsets[ns][ne]];
        }
        Index ns2;
        for(ns2 = 0; ns2 < subsets.size(); ++ns2) {
          if(contains_all(subsets[ns2], perm_sub)) {
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

    std::pair<std::vector<IrrepWedge>, std::vector<IrrepInfo> > irrep_wedges(SymGroup const &head_group, SymGroupRepID id) {
      return irrep_wedges(head_group.master_group().representation(id), head_group);
    }

    //*******************************************************************************************

    std::pair<std::vector<IrrepWedge>, std::vector<IrrepInfo> > irrep_wedges(SymGroupRep const &_rep, SymGroup const &head_group) {

      std::vector<IrrepInfo> irreps = irrep_decomposition(_rep, head_group);
      std::vector<IrrepWedge> wedges;
      wedges.reserve(sdirs.size());
      double best_proj, tproj;

      for(IrrepInfo const &irrep : IrrepInfo) {
        //1D irreps directions can have positive and negative directions, but we only want to include one.
        //only 1D irreps can have singly-degenerate directions and two singly-degenerate directions indicate
        //the same vector duplicated in positive and negative direction (because they are not equivalent by symmetry)
        //If irrep.directions[0] is singly degenerate (orbits size == 1) then irrepdim is 1 and we only need one direction to
        //define wedge
        wedges.emplace_back(Eigen::MatrixXd::Zero(irrep.vector_dim(), irrep.irrep_dim()), {});

        wedges.back().axes.col(0) = irrep.directions[0][0];
        wedges.back().mult.push_back(irrep.directions[0].size());
        for(Index i = 1; i < irrep.irrep_dim(); i++) {
          Index j_best = 0;
          best_proj = (wedges.back().axes.transpose() * irrep.directions[i][0]).sum();
          for(Index j = 1; j < irrep.directions[i].size(); j++) {
            tproj = (wedges().axes.transpose() * irrep.directions[i][j]).sum();
            if(tproj > best_proj) {
              best_proj = tproj;
              j_best = j;
            }
          }

          wedges.back().axes.col(i) = irrep.directions[i][j_best];
          wedges.back().mult.push_back(irrep.directions[i].size());
        }
      }
      return std::pair<std::vector<IrrepWedge>, std::vector<IrrepInfo> >(std::move(wedges), std::move(irreps));
    }

    //*******************************************************************************************

    std::pair<std::vector<SubWedge>, std::vector<IrrepInfo> > symrep_subwedges(SymGroup const &head_group, SymGroupRepID id) {
      return symrep_subwedges(head_group.master_group().representation(id), head_group);
    }

    //*******************************************************************************************

    std::pair<std::vector<SubWedge>, std::vector<IrrepInfo> > symrep_subwedges(SymGroupRep const &_rep, SymGroup const &head_group) {
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

      auto init_wedges_and_irreps = irrep_wedges(head_group, id);
      std::vector<IrrepWedge> const &init_wedges = init_wedges_and_irreps.first;


      std::vector<SubWedge> result;

      //irrep_wedge_orbits[w] is orbit of wedges[w]
      multivector<IrrepWedge>::X<2> irrep_wedge_orbits;
      irrep_wedge_orbits.reserve(init_wedges.size());
      //max_equiv[w] is irrep_wedge_orbits[w].size()-1
      std::vector<Index> max_equiv;
      max_equiv.reserve(init_wedges.size());

      for(IrrepWedge const &wedge : init_wedges) {
        irrep_wedge_orbits.push_back({wedge});

        //Start getting orbit of wedges[w]
        for(Index p = 0; p < head_group.size(); p++) {
          IrrepWedge test_wedge((*(_rep[head_group[p].index()]->MatrixXd()))*wedge.axes,
                                wedge.mult);

          if(contains(irrep_wedge_orbits.back(), test_wedge, irrep_wedge_compare))
            continue;
          irrep_wedge_orbits.back().push_back(test_wedge);
        }

        max_equiv.push_back(irrep_wedge_orbits.back().size() - 1);
      }
      max_equiv[0] = 0;

      //Counter over combinations of equivalent wedges
      Counter<std::vector<Index> > wcount(std::vector<Index>(init_wedges.size(), 0),
                                          max_equiv,
                                          std::vector<Index>(init_wedges.size(), 1));

      multivector<IrrepWedge>::X<3> tot_wedge_orbits;
      for(; wcount; ++wcount) {
        std::vector<IrrepWedge> twedge = init_wedges;
        for(Index i = 1; i < init_wedges.size(); i++)
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
        for(Index p = 0; p < head_group.size(); p++) {
          for(Index i = 0; i < twedge.size(); i++)
            twedge[i].axes = (*(_rep[head_group[p].index()]->MatrixXd())) * result.back().irrep_wedges()[i].axes;
          if(!contains(tot_wedge_orbits.back(), twedge, tot_wedge_compare)) {
            //std::cout << "Adding subwedge!\n";
            tot_wedge_orbits.back().push_back(twedge);
          }
        }
      }

      return std::pair<std::vector<SubWedge>, std::vector<IrrepInfo> >(std::move(result), std::move(init_wedges_and_irreps.second));
    }

  }
}
