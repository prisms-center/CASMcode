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
#include "casm/casm_io/stream_io/container.hh"
#include "casm/misc/algorithm.hh"
#include "casm/container/Counter.hh"

namespace CASM {

  //*******************************************************************************************
  Eigen::MatrixXd block_shape_matrix(SymGroupRep const &_rep, const SymGroup &head_group) {
    if(!_rep.size() || !head_group.size() || !_rep.MatrixXd(head_group[0]))
      return Eigen::MatrixXd();

    Eigen::MatrixXd block_shape(_rep.MatrixXd(head_group[0])->cwiseProduct(*_rep.MatrixXd(head_group[0])));

    for(Index i = 1; i < head_group.size(); i++) {
      block_shape += _rep.MatrixXd(head_group[i])->cwiseProduct(*_rep.MatrixXd(head_group[i]));
    }

    return block_shape;
  }

  //*******************************************************************************************

  Index num_blocks(SymGroupRep const &_rep, const SymGroup &head_group) {
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
  Eigen::MatrixXd irrep_symmetrizer(SymGroupRep const &_rep,
                                    const SymGroup &head_group,
                                    double vec_compare_tol) {
    Index dim = _rep.MatrixXd(0)->cols();
    return irrep_symmetrizer(_rep, head_group, Eigen::MatrixXd::Identity(dim, dim), vec_compare_tol);
  }
  //*******************************************************************************************

  //assumes that representation is real-valued and irreducible
  Eigen::MatrixXd irrep_symmetrizer(SymGroupRep const &_rep,
                                    const SymGroup &head_group,
                                    Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                                    double vec_compare_tol) {

    // Result of sdirs is ordered by symmetry breaking properties of different directions.
    // We won't do any element-based comparisons after this point, in order to ensure that
    // our choice of axes is totally determined by symmetry
    multivector<Eigen::VectorXd>::X<2>  sdirs = special_irrep_directions(_rep, head_group, _subspace, vec_compare_tol);

    // Four strategies, in order of desparation
    // 1) find a spanning set of orthogonal axes within a single orbit of special directions
    // 2) find a spanning set of orthogonal axes within the total set of special directions
    // 3) perform qr decomposition on lowest-multiplicity orbit to find spanning set of axes
    // 4) find column-echelon form of _subspace matrix to get a sparse/pretty set of axes
    Eigen::MatrixXd result;
    Index dim = _subspace.cols();
    Index min_mult = 10000;
    bool orb_orthog = false;
    bool tot_orthog = false;

    Eigen::MatrixXd axes, orb_axes, tot_axes;
    Index tot_col(0);
    tot_axes.setZero(_subspace.rows(), dim);

    for(auto const &orbit : sdirs) {
      // Strategy 1
      if((orb_orthog && orbit.size() < min_mult) || !orb_orthog) {
        orb_axes.setZero(_subspace.rows(), dim);
        Index col = 0;
        for(auto const &el : orbit) {
          if(almost_zero((el.transpose()*orb_axes).eval(), vec_compare_tol)) {
            if(col < orb_axes.cols())
              orb_axes.col(col++) = el;
            else {
              std::stringstream errstr;
              errstr << "Error in irrep_symmtrizer(). Constructing coordinate axes from special directions of space spanned by row vectors:\n "
                     << _subspace.transpose()
                     << "\nAxes collected thus far are the row vectors:\n" << orb_axes.transpose()
                     << "\nBut an additional orthogonal row vector has been found:\n" << el.transpose()
                     << "\nindicating that subspace matrix is malformed.";
              throw std::runtime_error(errstr.str());
            }
          }
        }
        //std::cout << "PLAN A-- col: " << col << "; dim: " << dim << "; min_mult: " << min_mult << "; orthog: " << orb_orthog << ";\naxes: \n" << axes << "\norb_axes: \n" << orb_axes << "\n\n";
        if(col == dim) {
          orb_orthog = true;
          min_mult = orbit.size();
          axes = orb_axes;
        }
      }

      // Greedy(ish) implementation of strategy 2 -- may not find a solution, even if it exists
      if(!orb_orthog && !tot_orthog) {
        for(auto const &el : orbit) {
          if(almost_zero((el.transpose()*tot_axes).eval(), vec_compare_tol)) {
            if(tot_col < tot_axes.cols())
              tot_axes.col(tot_col++) = el;
            else {
              std::stringstream errstr;
              errstr << "Error in irrep_symmtrizer(). Constructing coordinate axes from special directions of space spanned by row vectors:\n "
                     << _subspace.transpose()
                     << "\nAxes collected thus far are the row vectors:\n" << tot_axes.transpose()
                     << "\nBut an additional orthogonal row vector has been found:\n" << el.transpose()
                     << "\nindicating that subspace matrix is malformed.";
              throw std::runtime_error(errstr.str());
            }
          }
        }
        //std::cout << "PLAN B-- col: " << tot_col << "; dim: " << dim << "; min_mult: " << min_mult << "; orthog: " << tot_orthog << ";\naxes: \n" << axes << "\ntot_axes: \n" << tot_axes << "\n\n";
        if(tot_col == dim) {
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
        axes = Eigen::MatrixXd(orb_axes.colPivHouseholderQr().matrixQ()).leftCols(dim);
        //std::cout << "\norb_axes: \n" << orb_axes << "\n\n";
      }
    }
    //std::cout << "axes: \n" << axes << "\n";
    if(axes.cols() == 0) {
      SubspaceSymCompare<Eigen::MatrixXd, EigenSymRepApply<Eigen::MatrixXd>> tcompare(vec_compare_tol, _rep);
      result = *(tcompare.canonical_transform(_subspace)->MatrixXd());
    }
    else {
      // Strategy 4
      result = _subspace.colPivHouseholderQr().solve(axes);
    }

    //std::cout << "result: \n" << result << "\n"
    //<< "_subspace*result: \n" << _subspace*result << "\n";
    return result;//.transpose();
  }

  //*******************************************************************************************
  multivector< Eigen::VectorXd >::X<3>
  special_total_directions(SymGroupRep const &_rep,
                           SymGroup const &head_group) {

    std::pair<Eigen::MatrixXd, std::vector<Index>> transmat = irrep_trans_mat_and_dims(_rep,
                                                              head_group,
    [&](Eigen::Ref<const Eigen::MatrixXd> const & _subspace) {
      return Eigen::MatrixXd::Identity(_subspace.cols(), _subspace.cols());
    });
    multivector< Eigen::VectorXd >::X<3> result;
    result.reserve(transmat.second.size());

    Index l = 0;
    for(Index sdim : transmat.second) {
      Eigen::MatrixXd subspace = transmat.first.block(l, 0, sdim, transmat.first.cols()).transpose();
      result.push_back(special_irrep_directions(_rep, head_group, subspace, TOL));
      l += sdim;
    }
    return result;
  }

  //*******************************************************************************************
  /// Returns array of orbits of high-symmetry directions in the vector space on
  /// which this representation is defined. This routine is different from special_total_directions
  /// in that it does not rely on the character tables to generate the special directions. The
  /// routine also adopts the faster method to generating the irreducible transformation matrix

  multivector< Eigen::VectorXd >::X<2>
  special_irrep_directions(SymGroupRep const &_rep,
                           SymGroup const &head_group,
                           double vec_compare_tol) {
    Index dim = _rep.MatrixXd(0)->cols();
    return special_irrep_directions(_rep, head_group, Eigen::MatrixXd::Identity(dim, dim), vec_compare_tol);

  }
  //*******************************************************************************************
  /// Returns array of orbits of high-symmetry directions in the vector space on
  /// which this representation is defined. This routine is different from special_total_directions
  /// in that it does not rely on the character tables to generate the special directions. The
  /// routine also adopts the faster method to generating the irreducible transformation matrix

  multivector< Eigen::VectorXd >::X<2>
  special_irrep_directions(SymGroupRep const &_rep,
                           SymGroup const &head_group,
                           Eigen::Ref<const Eigen::MatrixXd> const &_subspace,
                           double vec_compare_tol,
                           bool all_subgroups) {
    auto sgroups = (all_subgroups ? head_group.subgroups() : head_group.small_subgroups());

    std::vector< Eigen::VectorXd > tdirs;
    Eigen::MatrixXd R;
    Index dim = _rep.MatrixXd(0)->cols();

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
      Eigen::MatrixXd Q = QR.matrixQ();
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

    using SymCompareType = DirectionSymCompare<Eigen::VectorXd, EigenSymRepApply<Eigen::VectorXd>>;
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
      //std::cout << "Orbit " << o << ": \n";
      //for(auto const & v : orbit){
      //std::cout << v.transpose() << "\n";
      //}
      //std::cout << "---------\n";
      result.emplace_back(orbit.begin(), orbit.end());
    }

    if(all_subgroups || !result.empty())
      return result;
    else
      return special_irrep_directions(_rep, head_group, _subspace, vec_compare_tol, true);
  }

  //*******************************************************************************************
  std::vector<std::vector< Eigen::VectorXd> > special_irrep_directions_old(SymGroupRep const &_rep, const SymGroup &head_group) {
    if(!_rep.size() || !(_rep.MatrixXd(head_group[0]))) {
      default_err_log() << "CRITICAL ERROR: In special_irrep_directions() called on imcompatible SymGroupRep.\n Exiting...\n";
      exit(1);
    }

    Index i, j;

    std::vector<std::vector<Eigen::VectorXd> > sdirs; //std::vector of orbits of special directions

    Eigen::VectorXd tdir;

    auto const &isubs(head_group.subgroups());
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(*(_rep.MatrixXd(head_group[0])));

    Index dim = _rep.MatrixXd(head_group[0])->cols();

    QR.setThreshold(10 * TOL);

    // i loops over subgroup "orbits". There are 0 to 'dim' special directions associated with each orbit.
    Eigen::MatrixXd Reynolds(*(_rep.MatrixXd(0)));
    for(i = 0; i < isubs.size(); i++) {
      if(isubs[i].size() == 0 || isubs[i].begin()->size() == 0) {
        default_err_log() << "CRITICAL ERROR: In special_irrep_directions(), attempting to use a zero-size subgroup.\n Exiting...\n";
        exit(1);
      }

      Reynolds.setZero();
      // j loops over elements of "prototype" subgroup to get the Reynold's operator for the vector
      // space on which the representation is defined
      for(Index op : * (isubs[i].begin())) {
        Reynolds += *_rep.MatrixXd(head_group[op]);
      }

      Reynolds /= double(isubs[i].begin()->size());

      //Need this because Eigen computes nonzero rank for matrix filled with small numbers close to zero.
      //QR.setThreshold() doesn't help
      if(almost_zero(Reynolds))
        continue;

      // Column space of Reynold's matrix is invariant under the subgroup
      QR.compute(Reynolds);

      // We're only interested in 1-D invariant spaces
      if(QR.rank() != 1)
        continue;

      //QR.matrixQ().col(0) is a special direction

      Eigen::MatrixXd matrixQ(QR.householderQ());

      //std::cout << "QR matrixR:\n" << Eigen::MatrixXd(QR.matrixQR().template triangularView<Eigen::Upper>()) << "\n";
      //std::cout << "QR matrixQ:\n" << matrixQ << "\n";


      auto vector_almost_equal = [&](const Eigen::VectorXd & val1, const Eigen::VectorXd & val2) {
        return almost_equal(val1, val2);
      };

      // See if QR.matrixQ().col(0) has been found
      for(j = 0; j < sdirs.size(); j++) {
        if(contains(sdirs[j], Eigen::VectorXd(matrixQ.col(0)), vector_almost_equal)) {
          break;
        }
      }
      if(j == sdirs.size()) {

        //Get equivalents
        sdirs.push_back(std::vector<Eigen::VectorXd>());
        for(j = 0; j < head_group.size(); j++) {
          tdir = (*_rep.MatrixXd(head_group[j])) * matrixQ.col(0);
          if(!contains(sdirs.back(), tdir, vector_almost_equal))
            sdirs.back().push_back(tdir);
        }
      }

      for(j = 0; j < sdirs.size(); j++) {
        if(contains(sdirs[j], Eigen::VectorXd(-matrixQ.col(0)), vector_almost_equal)) {
          break;
        }
      }
      if(j == sdirs.size() && dim > 1) {

        //Get equivalents
        sdirs.push_back(std::vector<Eigen::VectorXd>());
        for(j = 0; j < head_group.size(); j++) {
          tdir = -(*_rep.MatrixXd(head_group[j])) * matrixQ.col(0);
          if(!contains(sdirs.back(), tdir, vector_almost_equal))
            sdirs.back().push_back(tdir);
        }
      }
    }
    return sdirs;
  }

  //*******************************************************************************************

  // Similar to algorithm for finding special directions, except we also find special planes,
  // n-planes, etc. The matrices that are returned have column vectors that span the special
  // subspaces.  The subspaces are arranged in orbits of subspaces that are equivalent by symmetry

  std::vector<std::vector< Eigen::MatrixXd> > special_subspaces(SymGroupRep const &_rep, const SymGroup &head_group) {
    if(!_rep.size() || !_rep.MatrixXd(0)) {
      default_err_log() << "CRITICAL ERROR: In special_subspaces() called on imcompatible SymGroupRep.\n Exiting...\n";
      exit(1);
    }

    Index i, j, k;

    std::vector<std::vector<Eigen::MatrixXd> > ssubs; //std::vector of orbits of special directions

    Eigen::MatrixXd tsub, ttrans;

    // Operation indices of subgroups
    auto const &sg_op_inds(head_group.subgroups());

    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> QR(*(_rep.MatrixXd(0)));

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
      tsub = QR.matrixQ().leftCols(QR.rank());

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
  std::vector<Index> num_each_real_irrep(SymGroupRep const &_rep, const SymGroup &head_group, bool verbose) {
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

  std::vector<Index> num_each_irrep(SymGroupRep const &_rep, const SymGroup &head_group, bool verbose) {
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

  std::vector<SymGroupRepID> irrep_IDs(SymGroupRep const &_rep, const SymGroup &head_group) {
    std::vector<Index> irrep_decomp(num_each_real_irrep(_rep, head_group));
    std::vector<SymGroupRepID> result;
    for(Index i = 0; i < irrep_decomp.size(); i++) {
      if(irrep_decomp[i] == 0)
        continue;

      if(head_group.get_irrep_ID(i).empty())
        calc_new_irreps(_rep, head_group);

      if(head_group.get_irrep_ID(i).empty()) {
        default_err_log() << "CRITICAL ERROR: In irrep_IDs(), cannot resolve irrep " << i << " which has multiplicity " << irrep_decomp[i] << ". Exiting...\n";
        assert(0);
        exit(1);
      }
      auto tail = std::vector<SymGroupRepID>(irrep_decomp[i], head_group.get_irrep_ID(i));
      result.insert(result.end(), tail.begin(), tail.end());
    }

    return result;

  }


  //*******************************************************************************************
  //
  // assumes that representation is real-valued. This means that complex-valued
  // irreps in the character table get combined into single 2-D irreps

  void calc_new_irreps(SymGroupRep const &_rep, const SymGroup &head_group, int max_iter) {

    if(!_rep.size() || !_rep.MatrixXd(0)) {
      default_err_log() << "WARNING:  In calc_new_irreps, size of representation is " << _rep.size() << " and MatrixXd address is " << _rep.MatrixXd(0) << std::endl
                        << "          No valid irreps will be returned.\n";
      return;
    }
    if(!head_group.size()) {
      assert(0 && "In calc_new_irreps(), passed an empty SymGroup!");
      return;
    }

    MasterSymGroup const *master_ptr(NULL);
    if(_rep.has_valid_master())
      master_ptr = &_rep.master_group();
    else if(head_group[0].has_valid_master())
      master_ptr = &head_group[0].master_group();
    else {
      default_err_log() << "CRITICAL ERROR: In calc_new_irreps(), there is no MasterSymGroup available!\n"
                        << "                Exiting...\n";
    }
    int dim((_rep.MatrixXd(0))->rows());

    std::vector<Index> i_decomp(num_each_real_irrep(_rep, head_group));
    int Nirrep(std::accumulate(i_decomp.begin(), i_decomp.end(), 0));
    if(Nirrep == 0) {
      default_err_log() << "In calc_new_irreps(), there are no valid irreps, so no valid irreps will be returned.\n";
      return;
    }
    if(Nirrep == 1) {
      //std::cout << "There is only 1 irrep, so I'm returning a symmetrized copy of *this\n";
      head_group.set_irrep_ID(find_index(i_decomp, 1), master_ptr->add_representation(coord_transformed_copy(_rep, irrep_symmetrizer(_rep, head_group, TOL))));
      return;
    }

    //std::cout << "This representation contains " <<  Nirrep << " irreps... About to build commuters!\n";
    std::vector<Eigen::MatrixXd> commuters;

    // Build 'commuters', which span space of matrices that commute with SymOpReps
    for(int i = 0; i < dim; i++) {
      for(int j = i; j < dim; j++) {
        Eigen::MatrixXd tmat(Eigen::MatrixXd::Zero(dim, dim)), tcommute(Eigen::MatrixXd::Zero(dim, dim));

        //We only build symmetric commuters, which assures we only select out real-valued representations
        tmat(i, j) = 1;
        tmat(j, i) = 1;
        //apply reynolds operator
        for(Index ns = 0; ns < head_group.size(); ns++) {
          tcommute += (*(_rep.MatrixXd(head_group[ns]))) * tmat * (*(_rep.MatrixXd(head_group[ns]))).transpose();
        }
        //Do Gram-Shmidt while building 'commuters'
        for(Index nc = 0; nc < commuters.size(); nc++) {
          double tproj((commuters[nc].array()*tcommute.array()).sum()); //Frobenius product
          tcommute -= tproj * commuters[nc];
        }
        double tnorm(tcommute.norm());
        if(tnorm > TOL) {
          commuters.push_back(tcommute / tnorm);
          //std::cout << commuters.back() << "\n\n";
        }
      }
    }
    for(Index nc = 0; nc < commuters.size(); nc++) {
      commuters[nc] *= (dim * sqrt(dim)); //magnify the range of eigenvalues
    }
    //Finished building commuters

    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    double min_diff(-1);
    int Nstep(0), Nfound(0), max_found = 0;
    Eigen::MatrixXd tmat(dim, dim), trans_mat(dim, dim);
    int rseed(5284163); //use constant seed for reproducibility
    Eigen::VectorXd mags(commuters.size());
    while((max_found < Nirrep || min_diff < 0.05 || Nstep < 10) && Nstep < max_iter) {
      Nstep++;
      Nfound = 1;
      double tmindiff(1e+8), tdiff(0);

      tmat.setZero();
      for(EigenIndex i = 0; i < mags.size(); i++)
        mags[i] = ran0(rseed);
      mags.normalize();

      for(Index nc = 0; nc < commuters.size(); nc++) {
        tmat += mags[nc] * commuters[nc];
      }
      svd.compute(tmat, Eigen::ComputeFullU);

      for(EigenIndex ne = 1; ne < svd.singularValues().size(); ne++) {
        tdiff = std::abs(svd.singularValues()[ne - 1] - svd.singularValues()[ne]);
        if(tdiff > TOL) {
          Nfound++;
          tmindiff = std::min(tdiff, tmindiff);
        }
      }

      if(Nfound > max_found || (Nfound == max_found && tmindiff > min_diff)) {
        max_found = Nfound;
        //trans_mat defined to go from old basis to new basis
        trans_mat = svd.matrixU().transpose();
        min_diff = tmindiff;
      }
    }
    //std::cout << "After " << Nstep << " iterations, found suitable decomposition \n ***U:" << svd.matrixU() << "\n\n ***S" << svd.singularValues() << "\n\n";
    //std::cout << "trans_mat is:\n" << trans_mat << "\n\n";
    if(Nstep >= max_iter) {
      default_err_log() << "CRITICAL ERROR: In calc_new_irreps, reached maximum iterations (" << max_iter << ") without finding all irreps.\n"
                        << "                I found " << max_found << " of " << Nirrep << " distinguishable within " << min_diff << " (0.1 required). Exiting...\n";
      assert(0);
      exit(1);

    }

    // Almost there.  We have a transformation matrix that block-diagonalizes the symmetry matrices. Now we apply it.

    //get the block_shape_matrix of the transformed representation
    //std::cout << "Transformed representation will be: \n";
    Eigen::MatrixXd block_shape(Eigen::MatrixXd::Zero(dim, dim));
    for(Index i = 0; i < head_group.size(); i++) {
      //std::cout << "Operation " << i << " of " << head_group.size() << ":\n";
      tmat = trans_mat * (*(_rep.MatrixXd(head_group[i]))) * trans_mat.transpose();
      //std::cout << tmat << "\n";
      block_shape.array() += tmat.array() * tmat.array();
    }
    //std::cout << "Block shape matrix of transformed representation is: \n" << block_shape << '\n';
    //now loop over the 1st off-diagonal of the block-shape matrix and look for breaks
    int i1(0);
    int nblocks = 0;
    for(int i2 = 1; i2 < dim; i2++) {
      if(block_shape(i2 - 1, i2) > 0.5)
        continue;
      nblocks++;
      //std::cout <<  "Found block number " << nblocks << " of " << Nirrep << '\n';
      //Found end of block
      tmat.resize(i2 - i1, dim);
      tmat = trans_mat.block(i1, 0, i2 - i1, dim);
      SymGroupRep new_irrep(coord_transformed_copy(_rep, tmat));

      std::vector<Index> tdecomp(num_each_real_irrep(new_irrep, head_group));
      if(std::accumulate(tdecomp.begin(), tdecomp.end(), 0) != 1) {
        default_err_log() << "CRITICAL ERROR: In calc_new_irreps(), representation identified as irreducible did not pass final tests. Exiting...\n";
        assert(0);
        exit(1);
      }
      Index i_irrep = find_index(tdecomp, 1);
      if(!head_group.get_irrep_ID(i_irrep).empty()) {
        i1 = i2;
        continue;
      }

      SymGroupRepID irrep_ID = master_ptr->add_representation(coord_transformed_copy(new_irrep, irrep_symmetrizer(new_irrep, head_group, TOL)));

      head_group.set_irrep_ID(i_irrep, irrep_ID);
      //std::cout << "Setting the rep_ID of irrep " << i_irrep << " to " << irrep_ID;
      //std::cout << " and validating -- " << head_group.get_irrep_ID(i_irrep) << '\n';

      i1 = i2;

    }

    // Process the final block
    nblocks++;
    //std::cout <<  "Found block number " << nblocks << " of " << Nirrep << '\n';

    tmat.resize(dim - i1, dim);
    tmat = trans_mat.block(i1, 0, dim - i1, dim);
    SymGroupRep new_irrep(coord_transformed_copy(_rep, tmat));

    //std::cout << "\nIts representation looks like this:\n";
    //new_irrep->print_MatrixXd(std::cout, head_group);
    //std::cout << '\n';
    std::vector<Index> tdecomp(num_each_real_irrep(new_irrep, head_group));
    if(std::accumulate(tdecomp.begin(), tdecomp.end(), 0) != 1) {
      default_err_log() << "CRITICAL ERROR: In calc_new_irreps(), representation identified as irreducible did not pass final tests. Exiting...\n";
      assert(0);
      exit(1);
    }
    Index i_irrep = find_index(tdecomp, 1);
    if(!head_group.get_irrep_ID(i_irrep).empty()) {
      return;
    }

    SymGroupRepID irrep_ID = master_ptr->add_representation(coord_transformed_copy(new_irrep, irrep_symmetrizer(new_irrep, head_group, TOL)));

    //std::cout << "Setting the rep_ID of irrep " << i_irrep << " to " << irrep_ID;
    head_group.set_irrep_ID(i_irrep, irrep_ID);
    //std::cout << " and validating -- " << head_group.get_irrep_ID(i_irrep) << '\n';

    return;

  }

  //*******************************************************************************************

  bool is_irrep(SymGroupRep const &_rep, const SymGroup &head_group) {
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

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd irrep_trans_mat_old(SymGroupRep const &_rep, const SymGroup &head_group) {
    return irrep_trans_mat_and_dims_old(_rep, head_group).first;
  }

  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  // Also populate 'subspaces' with lists of columns that form irreps
  std::pair<Eigen::MatrixXd, std::vector<Index>> irrep_trans_mat_and_dims_old(SymGroupRep const &_rep, const SymGroup &head_group) {

    //std::cout << "INSIDE IRREP_TRANS_MAT\n";
    const std::vector<bool > &complex_irrep(head_group.get_complex_irrep_list());
    //const std::vector<std::vector<std::complex<double> > > &char_table(_rep.master_group().character_table());
    const std::vector<std::vector<std::complex<double> > > &char_table(head_group.character_table());

    std::vector<Index> irrep_decomp(num_each_real_irrep(_rep, head_group));
    std::vector<Index> irrep_dims;
    //irrep_IDs() harvests any new irreps that appear in this representation
    std::vector<SymGroupRepID> IDs(irrep_IDs(_rep, head_group));

    int rep_dim = _rep.MatrixXd(head_group[0])->rows();


    // We will fill the columns of trans_mat and then return its transpose at the end
    Eigen::MatrixXd tproj(Eigen::MatrixXd::Zero(rep_dim, rep_dim)),
          tmat(Eigen::MatrixXd::Zero(rep_dim, rep_dim)),
          trans_mat(Eigen::MatrixXd::Zero(rep_dim, rep_dim));

    // set a relatively high threshold for determining matrix rank, since we know values should be close to 1
    Eigen::FullPivLU<Eigen::MatrixXd> LU;
    LU.setThreshold(0.001);
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> QR;
    QR.setThreshold(0.001);

    int col_count = 0;
    for(Index i = 0; i < irrep_decomp.size(); i++) {
      if(irrep_decomp[i] == 0)
        continue;

      int c_mult(1);
      if(complex_irrep[i])
        c_mult = 2;

      const SymGroupRep &t_irrep(head_group.get_irrep(i));

      int irrep_dim = c_mult * round(char_table[i][0].real());

      // make sure t_irrep is valid and has the correct dimension.
      assert(irrep_dim == (t_irrep.MatrixXd(head_group[0])->rows()));

      //get the (0,0) projection operator -> projects out vectors that transform as the first basis vector of irrep 'i'
      tproj.setZero();
      //std::cout << "PROJECTING OUT THE FOLLOWING IRREP\n";
      for(Index j = 0; j < head_group.size(); j++) {
        //std::cout << "OP " << j << ":\n";

        tproj += (*(t_irrep.MatrixXd(head_group[j])))(0, 0) * (*_rep.MatrixXd(head_group[j]));
      }
      tproj *= double(irrep_dim) / double(head_group.size());

      // tproj is now the (0,0) projection operator for t_irrep. That means its rows span
      // a subspace spanned by (0,0) basis functions of t_irrep, and we can do Gram Shmidt to find a minimal set

      //std::cout << "Irrep decomp is " << irrep_decomp << "\n";
      //std::cout << "For irrep " << i << " tproj is:\n" << tproj << "\n\n";

      // We want to get a set of spanning vectors that are sparse and orthonormal
      // first we do LU (i.e., Gaussian elimination) to clean things up a bit
      LU.compute(tproj);
      tmat = LU.matrixLU().triangularView<Eigen::Upper>();

      //std::cout << "Gaussian elimination yields:\n" << tmat << "\n\n";

      //find reduced-row-echelon form -- why is this not an Eigen routine?
      double tval;
      for(int j = LU.rank() - 1; j > 0; j--) {
        tval = tmat(j, j); // <- in case aliasing is an issue
        tmat.row(j) /= tval;
        for(int k = 0; k < j; k++) {
          tval = tmat(k, j);
          tmat.row(k) -= tval * tmat.row(j);
        }
      }
      // rearrange columns so that they match initial definition
      tmat = tmat * LU.permutationQ().inverse();

      //std::cout << "And finally:\n" << tmat << "\n\n";

      // rows of tmat span the space, but they aren't orthonormal
      // do QR (aka Gram-Shmidt), and first QR.rank() columns of QR will span the space
      QR.compute(tmat.transpose());
      if(QR.rank() != irrep_decomp[i]) {
        default_err_log() << "CRITICAL ERROR: In irrep_trans_mat(), did not find correct number of basis vectors to account for\n"
                          << "                multiplicity of irrep #" << i << ". There should be " << irrep_decomp[i] << " but I found " << QR.rank() << "\n"
                          << "                Exiting...\n";
        exit(1);
      }

      //std::cout << "QR yields:\n " << QR.matrixQ() << "\n\n";
      for(EigenIndex n = 0; n < QR.rank(); n++) {
        irrep_dims.push_back(irrep_dim);
        trans_mat.col(col_count + n * irrep_dim) = QR.matrixQ().col(n);
      }
      //std::cout << "Initial trans_mat is:\n" << trans_mat << "\n\n";

      // NOW, let's make the (0,j) projection matrices to project out the other basis functions
      for(int j = 1; j < irrep_dim; j++) {

        tproj.setZero();
        for(Index k = 0; k < head_group.size(); k++) {
          //tproj += (*(t_irrep->MatrixXd(head_group[k])))(0,j) * (*MatrixXd(head_group[k])); // <-- this is wrong, not totally sure why
          tproj += (*(t_irrep.MatrixXd(head_group[k])))(j, 0) * (*_rep.MatrixXd(head_group[k]));
        }
        tproj *= double(irrep_dim) / double(head_group.size());

        //std::cout << "(0, " << j << ") projection operator is:\n" << tproj << "\n\n";
        //apply (0,j) projection operator to trans_mat. It will project the 0'th basis vectors to the j'th
        tmat = tproj * QR.matrixQ();

        //store the result in trans_mat
        for(EigenIndex n = 0; n < QR.rank(); n++) {
          trans_mat.col(col_count + n * irrep_dim + j) = tmat.col(n);
        }
        //std::cout << "NOW trans_mat is:\n" << trans_mat << "\n\n";
      }

      //update col_count:
      col_count += irrep_decomp[i] * irrep_dim;
    }

    return std::pair<Eigen::MatrixXd, std::vector<Index>>(trans_mat.transpose(), irrep_dims);
  }

  bool rep_check(SymGroupRep const &_rep, SymGroup const &head_group, bool verbose) {
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
  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  std::pair<Eigen::MatrixXd, std::vector<Index>> irrep_trans_mat_and_dims(SymGroupRep const &_rep,
                                                                          SymGroup const &head_group,
                                                                          std::function<Eigen::MatrixXd(Eigen::Ref<const Eigen::MatrixXd> const &f_subspace)> symmetrizer_func,
  Eigen::Ref<const Eigen::MatrixXd> const &_subspace) {

    auto result = irrep_trans_mat_and_dims(_rep, head_group, symmetrizer_func);
    if(_subspace.rows() != result.first.cols()) {
      throw std::runtime_error("In irrep_trans_mat_and_dims, subspace matrix does not have proper number of rows (should have "
                               + std::to_string(result.first.cols()) + ", but has " + std::to_string(_subspace.rows()));
    }
    Eigen::VectorXd trans_mags = (result.first * _subspace).array().square().rowwise().sum();
    std::pair<Eigen::MatrixXd, std::vector<Index>> new_result;

    Index l = 0;
    Index keep = 0;
    for(Index i : result.second) {
      double sqmag(0.);
      Index firstl = l;
      for(; l < firstl + i; ++l)
        sqmag += trans_mags[l];
      if(sqmag > 0.5) {
        new_result.first.conservativeResize(new_result.first.rows() + i, result.first.cols());
        new_result.first.block(keep, 0, i, result.first.cols()) = result.first.block(firstl, 0, i, result.first.cols());
        new_result.second.push_back(i);
        keep += i;
      }
    }

    return new_result;
  }
  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  std::pair<Eigen::MatrixXd, std::vector<Index>> irrep_trans_mat_and_dims(SymGroupRep const &_rep,
                                                                          const SymGroup &head_group,
  std::function<Eigen::MatrixXd(Eigen::Ref<const Eigen::MatrixXd> const &f_subspace)> symmetrizer_func) {
    std::vector<Eigen::VectorXcd> char_table;
    std::vector<Index> irrep_dims;
    if(!_rep.size() || !head_group.size() || !_rep.MatrixXd(head_group[0])) {
      default_err_log() << "WARNING:  In calc_new_irreps, size of representation is " << _rep.size() << " and MatrixXd address is " << _rep.MatrixXd(head_group[0]) << std::endl
                        << "          No valid irreps will be returned.\n";
      return std::make_pair(Eigen::MatrixXd(), irrep_dims);
    }
    assert(rep_check(_rep, head_group, true) && "REPRESENTATION IS ILL-DEFINED!!");
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



    //std::cout << "Commuter are:\n";
    // Build 'commuters', which span space of real matrices that commute with SymOpReps
    // 'kernel' is the kernel of trans_mat, and as the loop progresses, 'kernel' shrinks and the rank for trans_mat increases
    //std::vector<Index> irrep_dims;
    //std::cout << "~~~~~~~Dimension is " << dim << "\n";
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
          //std::cout << "Counting indices: ";
          for(Index ns = 0; ns < head_group.size(); ns++) {
            //std::cout << head_group[ns].index() << "  ";
            tcommute += (*(_rep.MatrixXd(head_group[ns]))) * tmat * (*(_rep.MatrixXd(head_group[ns]))).transpose();
            //tcommute += phase[nph]*(_rep.MatrixXd(ns)->col(i)) * (_rep.MatrixXd(ns)->row(j))+std::conj(phase[nph])*(_rep.MatrixXd(ns)->col(j)) * (_rep.MatrixXd(ns)->row(i));
          }

          //Do Gram-Shmidt while building 'commuters'

          for(Index nc = 0; nc < commuters.size(); nc++) {
            cplx tproj((commuters[nc].array().conjugate()*tcommute.array()).sum()); //Frobenius product
            tcommute -= tproj * commuters[nc];
          }

          double tnorm((tcommute.array().conjugate()*tcommute.array()).sum().real());
          if(tnorm > TOL) {
            commuters.push_back(tcommute / tnorm);
            //std::cout << "\nCommuter " << commuters.size() << "\n";
            //std::cout << commuters.back().real() << "\n\n";
          }
          else continue;  // Attempt to construct the next commuter...


          //Finished building commuter now we can try to harvest irreps from it

          // construct trans_mat from the non-degenerate irreps obtained from the commuter

          Index nc = commuters.size() - 1;

          // magnify the range of eigenvalues to be (I think) independent of matrix dimension by multiplying by dim^{3/2}
          esolve.compute(double(dim)*sqrt(double(dim))*kernel.transpose()*commuters[nc]*kernel);
          //esolve.compute(double(dim)*sqrt(double(dim))*commuters[nc]);
          //std::cout << "KERNEL IS:\n" << kernel << "\n\n";
          //std::cout << esolve.eigenvalues().size() << " EIGENVALUES of commuter " << nc << " are \n" << esolve.eigenvalues().transpose() << "\n";
          std::vector<Index> subspace_dims = partition_distinct_values(esolve.eigenvalues());

          // Columns of tmat are orthonormal eigenvectors of commuter in terms of natural basis (they were calculated in terms of kernel as basis)
          tmat = kernel * (esolve.eigenvectors().householderQr().householderQ());
          //tmat = esolve.eigenvectors().householderQr().householderQ();


          // make transformed copy of the representation
          std::vector<Eigen::MatrixXcd> trans_rep(head_group.size());
          Eigen::MatrixXd block_shape(Eigen::MatrixXd::Zero(kernel.cols(), kernel.cols()));
          //std::cout << "Transformed representation is:\n";
          for(Index i = 0; i < head_group.size(); i++) {
            trans_rep[i] = tmat.adjoint() * (*_rep.MatrixXd(head_group[i])) * tmat;
            block_shape += (trans_rep[i].cwiseProduct(trans_rep[i].conjugate())).real();
            //std::cout << trans_rep[i] << "\n\n";
          }
          //std::cout << "Mini block shape:\n" << block_shape << "\n";
          Index last_i = 0;
          //std::cout << "subspace_dims is : " << subspace_dims << "\n";
          for(Index ns = 0; ns < subspace_dims.size(); ns++) {
            //std::cout << "ns is " << ns << "\n";
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

            //std::cout << "REPRESENTATION norm = "<< sqnorm << "\n";
            if(almost_equal(sqnorm, double(head_group.size()))) { // this representation is irreducible
              //std::cout << "THE REPRESENTATION IS IRREDUCIBLE with norm = "<< sqnorm << "\n";
              //std::cout << "Its characters are: " << char_vec << "\n\n";
              Eigen::MatrixXd ttrans_mat(dim, 2 * subspace_dims[ns]);
              //std::cout << "Vectors in this subspace are:\n" << tmat.block(0, last_i, dim, subspace_dims[ns]) << "\n\n";
              ttrans_mat.leftCols(subspace_dims[ns]) = sqrt(2.0) * tmat.block(0, last_i, dim, subspace_dims[ns]).real();
              ttrans_mat.rightCols(subspace_dims[ns]) = sqrt(2.0) * tmat.block(0, last_i, dim, subspace_dims[ns]).imag();
              //std::cout << "Separated into components:\n" << ttrans_mat << "\n\n";
              //Only append to trans_mat if the new columns extend the space (i.e., are orthogonal)

              if(almost_zero((ttrans_mat.transpose()*trans_mat).norm(), 0.001)) {
                qr.compute(ttrans_mat);
                //it seems stupid to use two different decompositions that do almost the same thing, but
                // HouseholderQR is not rank-revealing, and colPivHouseholder mixes up the columns of the Q matrix.
                colqr.compute(ttrans_mat);
                Index rnk = colqr.rank();
                //std::cout << "rank is " << rnk << " and R matrix is \n" << colqr.matrixR() << "\n";
                irrep_dims.push_back(subspace_dims[ns]);
                if(rnk == 2 * subspace_dims[ns])
                  irrep_dims.push_back(subspace_dims[ns]);
                ttrans_mat = Eigen::MatrixXd(qr.householderQ()).leftCols(rnk);
                //SymGroupRep t_rep(coord_transformed_copy(_rep, ttrans_mat.transpose()));
                //std::cout << "***Adding columns!\n" << Eigen::MatrixXd(qr.householderQ()).leftCols(rnk) << "\n\n";
                //trans_mat.block(0, Nfound, dim, rnk) = Eigen::MatrixXd(qr.householderQ()).leftCols(rnk);
                //trans_mat.block(0, Nfound, dim, rnk) = ttrans_mat * (irrep_symmetrizer(t_rep, head_group, tol)).transpose();
                //std::cout << "Nfound, dim, rnk: " << Nfound << ", " << dim << ", " << rnk << "\n\n ---------------\n";
                trans_mat.block(0, Nfound, dim, rnk) = ttrans_mat * symmetrizer_func(ttrans_mat);
                //std::cout << "trans_mat is \n" << trans_mat << "\n\n";
                //std::cout << "ttrans_mat is \n" << ttrans_mat << "\n\n";
                //std::cout << "symmetrizer_func is \n" << symmetrizer_func(ttrans_mat) << "\n\n";
                Nfound += rnk;
                found_new_irreps = true;
                char_table.push_back(char_vec);
                //std::cout << "trans_mat is now: \n" << trans_mat << "\n\n";
              }
            }
            else {
              //std::cout << "THE REPRESENTATION IS **REDUCIBLE with norm = "<< sqnorm << "\n";
              //std::cout << "It's characters are: " << char_array << "\n\n";
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

    // make transformed copy of the representation
    Eigen::MatrixXd block_shape(Eigen::MatrixXd::Zero(dim, dim));
    //std::cout << "Transformed representation is:\n";
    for(Index i = 0; i < head_group.size(); i++) {
      block_shape += (trans_mat.transpose() * (*_rep.MatrixXd(head_group[i])) * trans_mat).cwiseAbs2();
    }


    return std::make_pair(trans_mat.transpose(), irrep_dims);

  }

  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd irrep_trans_mat(SymGroupRep const &_rep, const SymGroup &head_group) {
    return irrep_trans_mat_and_dims(_rep,
                                    head_group,
    [&](Eigen::Ref<const Eigen::MatrixXd> const & _subspace) {
      return irrep_symmetrizer(_rep, head_group, _subspace, TOL); //.transpose();
    }).first;
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

  /// \brief Make copy of (*this) that is transformed so that axes are oriented along high-symmetry direction
  /// and confined to subspaces that transform as irreps.
  SymGroupRep symmetry_adapted_copy(SymGroupRep const &_rep, const SymGroup &head_group) {
    return coord_transformed_copy(_rep, irrep_trans_mat(_rep, head_group));
  }
  //*******************************************************************************************

  SymGroupRep subset_permutation_rep(const SymGroupRep &permute_rep, const std::vector<std::vector<Index>> &subsets) {
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

  SymGroupRep permuted_direct_sum_rep(const SymGroupRep &permute_rep, const std::vector<SymGroupRep const *> &sum_reps) {
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
      rep_dims[i] = (*sum_reps[i]).MatrixXd(0)->cols();
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

  SymGroupRep kron_rep(const SymGroupRep &LHS, const SymGroupRep &RHS) {
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

    std::vector<IrrepWedge> irreducible_wedges(const SymGroup &head_group, SymGroupRepID id) {
      SymGroupRep const &srep(head_group.master_group().representation(id));
      Index dim = srep.MatrixXd(head_group[0])->cols();

      //Handle for strain symrep
      SymGroupRep::RemoteHandle trep(head_group, id);

      multivector<Eigen::VectorXd>::X<3> sdirs = special_total_directions(srep, head_group);
      std::vector<IrrepWedge> result(sdirs.size());
      double best_proj, tproj;

      for(Index s = 0; s < sdirs.size(); s++) {
        //1D irreps directions can have positive and negative directions, but we only want to include one.
        //only 1D irreps can have singly-degenerate directions and two singly-degenerate directions indicate
        //the same vector duplicated in positive and negative direction (because they are not equivalent by symmetry)
        //If sdirs[s][0] is singly degenerate (orbits size == 1) then irrepdim is 1 and we only need one direction to
        //define wedge
        Index irrepdim = sdirs[s].size();
        if(irrepdim == 2 && sdirs[s][0].size() == 1)
          irrepdim = 1;

        result[s].axes = Eigen::MatrixXd::Zero(dim, irrepdim);
        result[s].axes.col(0) = sdirs[s][0][0];
        result[s].mult.push_back(sdirs[s][0].size());
        for(Index i = 1; i < irrepdim; i++) {
          Index j_best = 0;
          best_proj = (result[s].axes.transpose() * sdirs[s][i][0]).sum();
          for(Index j = 1; j < sdirs[s][i].size(); j++) {
            tproj = (result[s].axes.transpose() * sdirs[s][i][j]).sum();
            if(tproj > best_proj) {
              best_proj = tproj;
              j_best = j;
            }
          }

          result[s].axes.col(i) = sdirs[s][i][j_best];
        }
      }
      return result;
    }




    //*******************************************************************************************

    std::vector<SubWedge > symrep_subwedges(SymGroup const &head_group, SymGroupRepID id) {
      auto irrep_wedge_compare = [](const IrrepWedge & a, const IrrepWedge & b)->bool {
        return Eigen::almost_equal(a.axes, b.axes);
      };

      auto tot_wedge_compare = [irrep_wedge_compare](const std::vector<IrrepWedge> &a, const std::vector<IrrepWedge> &b)->bool {
        for(auto ita = a.begin(), itb = b.begin(); ita != a.end(); ++ita, ++itb)
          if(!irrep_wedge_compare(*ita, *itb))
            return false;
        //std::cout << "Equal: \n" << SubWedge(a).trans_mat() << ",\n" << SubWedge(b).trans_mat() << "\n\n";
        return true;
      };


      std::vector<SubWedge> result;
      SymGroupRep const &srep(head_group.master_group().representation(id));
      if(!srep.MatrixXd(0))
        throw std::runtime_error("In symrep_subwedges, SymGroupRep does not describe matrix representation");
      Index dim = srep.MatrixXd(0)->cols();

      //for(SymOp const &op : head_group) {
      //std::cout << "OP " << op.index() << ":\n" << *(srep[op.index()]->MatrixXd()) << "\n\n";
      //}

      std::vector<IrrepWedge> init_wedges = irreducible_wedges(head_group, id);
      //for(Index i=0; i<init_wedges.size(); ++i){
      //std::cout << "IrrepWedge #" << i+1 << ":\n" << init_wedges[i].axes.transpose() << "\n\n";
      //}


      //Handle for strain symrep
      SymGroupRep::RemoteHandle trep(head_group, id);
      //irrep_wedge_orbits[w] is orbit of wedges[w]
      multivector<IrrepWedge>::X<2> irrep_wedge_orbits;
      irrep_wedge_orbits.reserve(init_wedges.size());
      //max_equiv[w] is irrep_wedge_orbits[w].size()-1
      std::vector<Index> max_equiv;
      max_equiv.reserve(init_wedges.size());

      for(IrrepWedge const &wedge : init_wedges) {
        irrep_wedge_orbits.push_back({wedge});

        //Start getting orbit of wedges[w]
        for(Index p = 0; p < trep.size(); p++) {
          IrrepWedge test_wedge((*(trep[p]->MatrixXd()))*wedge.axes,
                                wedge.mult);

          if(contains(irrep_wedge_orbits.back(), test_wedge, irrep_wedge_compare))
            continue;
          irrep_wedge_orbits.back().push_back(test_wedge);
        }

        max_equiv.push_back(irrep_wedge_orbits.back().size() - 1);
        max_equiv[0] = 0;
      }



      //Counter over combinations of equivalent wedges

      Counter<std::vector<Index> > wcount(std::vector<Index>(init_wedges.size(), 0), max_equiv, std::vector<Index>(init_wedges.size(), 1));
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
        result.push_back({twedge});
        for(Index p = 0; p < trep.size(); p++) {
          for(Index i = 0; i < twedge.size(); i++)
            twedge[i].axes = (*(trep[p]->MatrixXd())) * result.back().irrep_wedges()[i].axes;
          if(!contains(tot_wedge_orbits.back(), twedge, tot_wedge_compare)) {
            //std::cout << "Adding subwedge!\n";
            tot_wedge_orbits.back().push_back(twedge);
          }
        }
      }

      return result;
    }

  }
}
