#include "casm/symmetry/SymGroupRep.hh"

#include <numeric>
#include "casm/external/Eigen/CASM_AddOns"
#include "casm/misc/CASM_math.hh"

#include "casm/container/Permutation.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"


namespace CASM {

  SymGroupRep::SymGroupRep(const SymGroupRep &RHS) :
    Array<SymOpRepresentation * > (RHS.size(), NULL) {
    (*this) = RHS;
  }

  //*******************************************************************************************

  SymGroupRep::~SymGroupRep() {
    clear();
  }

  //*******************************************************************************************

  SymGroupRep &SymGroupRep::operator=(const SymGroupRep &RHS) {
    m_master_group = RHS.m_master_group;
    clear();
    if(RHS.size() > 0 && has_valid_master() && master_group().size() != RHS.size()) {
      throw std::runtime_error("Invalid assignment of SymGroupRep.  Sizes are incompatible.\n");
    }
    resize(RHS.size(), NULL);
    for(Index i = 0; i < RHS.size(); i++) {
      if(RHS[i])
        set_rep(i, *RHS[i]);
    }
    return *this;
  }

  //*******************************************************************************************

  void SymGroupRep::set_master_group(const MasterSymGroup &master, const SymGroupRepID &_rep_ID) {
    m_master_group = &master;
    if(_rep_ID.empty() || &(master.representation(_rep_ID)) != this) {
      throw std::runtime_error(std::string("SymGroupRep::set_master_group() attempted to assign SymGroupRepID that does not match the current SymGroupRep!\n"));
    }
    m_rep_ID = _rep_ID;
    if(size() == 0)
      Array<SymOpRepresentation *>::resize(master.size());
    else if(size() == master.size()) {
      for(Index i = 0; i < size(); i++) {
        if(at(i))
          at(i)->set_identifiers(master, m_rep_ID, i);
      }
    }
    else {
      throw std::runtime_error(std::string("SymGroupRep::set_master_group() passed new master whose size is incompatible with the size\n") +
                               "of the current representation.\n");
    }
  }

  //*******************************************************************************************

  void SymGroupRep::set_rep(const SymOp &base_op, const SymOpRepresentation &new_rep) const {
    if(!has_valid_master()) {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::set_rep(), you are trying to assign the representation of a SymOp whose factor_group is not specified!\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }
    if(valid_index(base_op.index()))
      return set_rep(base_op.index(), new_rep);

    //else:
    return set_rep(master_group().find_periodic(base_op), new_rep);

  }

  //*******************************************************************************************

  void SymGroupRep::set_rep(const SymOpRepresentation &base_op, const SymOpRepresentation &new_rep) const {
    if(!has_valid_master()) {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::set_rep(), you are trying to assign the representation of a SymOpRepresentation whose factor_group is not specified!\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }
    if(valid_index(base_op.index()))
      set_rep(base_op.index(), new_rep);
    else {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::set_rep(), you are trying to assign the representation of a SymOpRepresentation whose index is not specified!\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

  }

  //*******************************************************************************************

  void SymGroupRep::set_rep(Index op_index, const SymOpRepresentation &new_rep) const {

    assert(valid_index(op_index) && op_index < size() && "In SymGroupRep::set_rep(), reference representation is improperly initialized.");
    if(at(op_index)) {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::set_rep(), representation already exists for operation " << op_index << ".\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    SymOpRepresentation *tcopy = new_rep.copy();
    if(has_valid_master())
      tcopy->set_identifiers(master_group(), m_rep_ID, op_index);

    // avoid doing this elsewhere in CASM
    const_cast<SymOpRepresentation *&>(at(op_index)) = tcopy;

  }

  //*******************************************************************************************

  SymGroupRepID SymGroupRep::add_copy_to_master() const {
    return master_group().add_representation(*this);
  }
  //*******************************************************************************************

  void SymGroupRep::clear() {
    for(Index i = 0; i < size(); i++)
      delete at(i);
    Array<SymOpRepresentation *>::clear();
  }

  //John G 050513

  //*******************************************************************************************

  void SymGroupRep::print_permutation(std::ostream &stream) const {
    for(Index i = 0; i < this->size(); i++) {
      if(at(i)->get_permutation())
        stream << *(at(i)->get_permutation()) << '\n';
      else
        stream << "nullptr\n";
    }
    return;
  }

  //\John G

  //*******************************************************************************************

  void SymGroupRep::print_MatrixXd(std::ostream &stream) const {
    for(Index i = 0; i < size(); i++) {
      stream << "SymRep Matrix " << i + 1 << " of " << size() << ":\n";
      if(at(i)->get_MatrixXd())
        stream << *(at(i)->get_MatrixXd()) << '\n';
      else
        stream << "nullptr\n";
      stream << "\n";
    }
    std::cout << "\n";
    return;
  }

  //*******************************************************************************************

  void SymGroupRep::print_MatrixXd(std::ostream &stream, const SymGroup &head_group) const {
    for(Index i = 0; i < head_group.size(); i++) {
      stream << "SymRep Matrix " << i + 1 << " of " << head_group.size() << ":\n";
      if(at(head_group[i].index())->get_MatrixXd())
        stream << *(at(head_group[i].index())->get_MatrixXd()) << '\n';
      else
        stream << "nullptr\n";
      stream << "\n";
    }
    return;
  }

  //*******************************************************************************************
  Eigen::MatrixXd SymGroupRep::block_shape_matrix() const {
    if(!size() || !get_MatrixXd(0))
      return Eigen::MatrixXd();

    Eigen::MatrixXd block_shape(get_MatrixXd(0)->cwiseProduct(*get_MatrixXd(0)));

    for(Index i = 1; i < size(); i++) {
      block_shape += get_MatrixXd(i)->cwiseProduct(*get_MatrixXd(i));

    }
    return block_shape;
  }
  //*******************************************************************************************
  Eigen::MatrixXd SymGroupRep::block_shape_matrix(const SymGroup &head_group) const {
    if(!size() || !head_group.size() || !get_MatrixXd(head_group[0].index()))
      return Eigen::MatrixXd();

    Eigen::MatrixXd block_shape(get_MatrixXd(head_group[0].index())->cwiseProduct(*get_MatrixXd(head_group[0].index())));

    for(Index i = 1; i < head_group.size(); i++) {
      block_shape += get_MatrixXd(head_group[i].index())->cwiseProduct(*get_MatrixXd(head_group[i].index()));
    }

    return block_shape;
  }

  //*******************************************************************************************
  Index SymGroupRep::num_blocks() const {
    Eigen::MatrixXd bmat(block_shape_matrix());

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

  Index SymGroupRep::num_blocks(const SymGroup &head_group) const {
    Eigen::MatrixXd bmat(block_shape_matrix(head_group));

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
  void SymGroupRep::push_back_copy(const SymOpRepresentation &_pushed) {
    if(has_valid_master()) {
      throw std::runtime_error("SymGroupRep::push_back_copy() only callable on SymGroupRep that has no master\n");
    }
    Array<SymOpRepresentation *>::push_back(_pushed.copy());

  }
  //*******************************************************************************************

  Eigen::MatrixXd const *SymGroupRep::get_MatrixXd(Index i) const {
    return at(i)->get_MatrixXd();
  }

  //*******************************************************************************************

  Eigen::MatrixXd const *SymGroupRep::get_MatrixXd(const SymOpRepresentation &op) const {
    return at(op.index())->get_MatrixXd();
  }

  //*******************************************************************************************

  Permutation const *SymGroupRep::get_permutation(Index i) const {
    return at(i)->get_permutation();
  }

  //*******************************************************************************************

  Permutation const *SymGroupRep::get_permutation(const SymOpRepresentation &op) const {
    return at(op.index())->get_permutation();
  }


  //*******************************************************************************************
  // Calculates new SymGroupRep that is the results of performing coordinate transformation specified by trans_mat
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  SymGroupRep SymGroupRep::coord_transformed_copy(const Eigen::MatrixXd &trans_mat) const {
    SymGroupRep new_rep(master_group());
    if(!size())
      return new_rep;
    if(at(0) && !(at(0)->get_MatrixXd())) {
      std::cerr << "CRITICAL ERROR: Trying to perform matrix transformation on a non-matrix SymRep. Exiting...\n";
      assert(0);
      exit(1);
    }
    Eigen::MatrixXd rightmat;
    rightmat = Eigen::JacobiSVD<Eigen::MatrixXd>(trans_mat.transpose()).solve(Eigen::MatrixXd::Identity(trans_mat.rows(), trans_mat.rows()));

    for(Index i = 0; i < size(); i++) {
      if(!at(i))
        continue;

      new_rep.set_rep(i, SymMatrixXd(trans_mat * (*(at(i)->get_MatrixXd())) * rightmat));
    }
    return new_rep;
  }

  //*******************************************************************************************

  //assumes that representation is real-valued and irreducible
  Eigen::MatrixXd SymGroupRep::_symmetrized_irrep_trans_mat(const SymGroup &head_group) const {

    // High symmetry subspace matrices (columns span subspaces)
    Array<Array< Eigen::MatrixXd> > ssubs;
    ssubs = calc_special_subspaces(head_group);
    /*for(Index i = 0; i < ssubs.size(); i++) {
      std::cout << "In SymGroupRep::coord_symmetrized_copy -- subspace orbit " << i << " of " << ssubs.size() << " (dimensionality " << ssubs[i][0].cols() << "):\n";
      for(Index j = 0; j < ssubs[i].size(); j++) {
    std::cout << ssubs[i][j].transpose() << "\n";
      }
      std::cout << "\n";
    }
    */

    //std::cout << "\n\nSymGroupRep before reorientation is:\n";
    //print_MatrixXd(std::cout, head_group);
    //std::cout << '\n';

    Array<bool> sub_valid(ssubs.size(), true);
    Array<Index> sub_ranks(ssubs.size());


    Index dim(get_MatrixXd(head_group[0].index())->rows());
    Eigen::MatrixXd trans_mat(Eigen::MatrixXd::Zero(dim, dim)), tmat;
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> QR;
    QR.setThreshold(TOL);

    // calculate the rank of each subspace orbit
    for(Index i = 0; i < ssubs.size(); i++) {
      //Fill tmat by concatenating columns of subspace matrices
      tmat.resize(dim, ssubs[i][0].cols()*ssubs[i].size());
      for(Index j = 0; j < ssubs[i].size(); j++) {
        tmat.block(0, j * ssubs[i][j].cols(), dim, ssubs[i][j].cols()) = ssubs[i][j];
      }
      sub_ranks[i] = QR.compute(tmat).rank();
    }

    Index i_best, n_cols(0);
    while(n_cols < dim) {

      // find the best subspace orbit to add
      i_best = -1;
      for(Index i = 0; i < ssubs.size(); i++) {
        // The subspace orbit must be valid (as we add subspaces, we invalidate them)
        if(!sub_valid[i])
          continue;

        if(!valid_index(i_best)) {
          i_best = i;
          continue;
        }

        //compare cols() also
        if(ssubs[i][0].cols() < ssubs[i_best][0].cols() // prefer 1D subspaces (high symmetry directions are better than high-symmetry planes)
           || (ssubs[i][0].cols() == ssubs[i_best][0].cols() && sub_ranks[i] < sub_ranks[i_best])  // prefer low-rank orbits (likely irreps)
           // prefer lower multiplicity (higher symmetry)
           || (ssubs[i][0].cols() == ssubs[i_best][0].cols() && sub_ranks[i] == sub_ranks[i_best] && ssubs[i].size() < ssubs[i_best].size())) {
          i_best = i;
        }
      }

      // make big matrix that combines trans_mat with best subspace
      // we fill it in this particular way so that the resulting basis is most likely to transform as an irreducible representation
      tmat.resize(dim, n_cols + ssubs[i_best].size()*ssubs[i_best][0].cols());
      tmat.leftCols(n_cols) = trans_mat.leftCols(n_cols);
      for(EigenIndex k = 0; k < ssubs[i_best][0].cols(); k++) {
        for(Index j = 0; j < ssubs[i_best].size(); j++) {
          tmat.col(n_cols + k * ssubs[i_best].size() + j) = ssubs[i_best][j].col(k);
        }
      }

      // number of basis vectors is now QR.rank()
      n_cols = QR.compute(tmat).rank();
      trans_mat.leftCols(n_cols) = Eigen::MatrixXd(tmat.householderQr().householderQ()).leftCols(n_cols);
      sub_valid[i_best] = false; // take i_best out of consideration

      if(n_cols == dim)
        break;

      // invalidate any subspaces that are redundant
      for(Index i = 0; i < ssubs.size(); i++) {
        if(!sub_valid[i])
          continue;
        tmat.resize(dim, n_cols + ssubs[i].size()*ssubs[i][0].cols());
        tmat.leftCols(n_cols) = trans_mat.leftCols(n_cols);
        for(Index j = 0; j < ssubs[i].size(); j++) {
          tmat.block(0, n_cols + j * ssubs[i][j].cols(), dim, ssubs[i][j].cols()) = ssubs[i][j];
        }

        if(QR.compute(tmat).rank() < n_cols + sub_ranks[i]) {
          sub_valid[i] = false;
        }
      }

    }
    //std::cout << "TRANSFORMATION MATRIX IS:\n" << trans_mat << "\n\n";

    return trans_mat.transpose();

  }

  //*******************************************************************************************
  /// Returns array of orbits of high-symmetry directions in the vector space on which this representation is defined.
  ///     defining: result = calc_special_directions(my_group)
  /// result[i][j] is a direction that is invariant to a subgroup of 'my_group'
  /// result[i][k] is a direction that is equivalent to result[i][j]
  /// For any result[i][j], -result[i][j] must also be a special direction, but the properties of the group determine whether
  /// it is in result[i] or in some other orbit result[l].
  /// If, for a result[i][j], -result[i][j] is not among the special directions, it means that result[i][j] is a 1d irrep
  /// that is invariant to 'my_group' (the irreducible 'wedge' spans that entire axis).  For a 2d irrep, the irreducible
  /// wedge must span <= one quadrant.  For a 3d irrep, it must span <= one octant
  multivector< Eigen::VectorXd>::X<3> SymGroupRep::calc_special_total_directions(const SymGroup &head_group)const {
    Array<SymGroupRepID> irrep_IDs = get_irrep_IDs(head_group);

    multivector<Eigen::VectorXd>::X<3> sdirs(irrep_IDs.size());
    Array<Array<Eigen::VectorXd> > irrep_dirs;
    Index irow(0);
    Eigen::MatrixXd trans_mat = get_irrep_trans_mat(head_group);
    for(Index i = 0; i < irrep_IDs.size(); i++) {
      RemoteHandle irrep(head_group, irrep_IDs[i]);
      Eigen::MatrixXd subtrans_mat = trans_mat.block(irow, 0, irrep.dim(), trans_mat.cols()).transpose();
      if(i == 0 || irrep_IDs[i] != irrep_IDs[i - 1]) {
        irrep_dirs = irrep->_calc_special_irrep_directions(head_group);
      }
      for(Index s = 0; s < irrep_dirs.size(); s++) {
        //std::cout << "Irrep dir " << i << ", " << s << ": \n";
        sdirs[i].push_back(std::vector<Eigen::VectorXd>());
        for(Index d = 0; d < irrep_dirs[s].size(); d++) {
          //std::cout << irrep_dirs[s][d].transpose() << std::endl;
          sdirs[i].back().push_back(subtrans_mat * irrep_dirs[s][d]);
        }
      }
      irow += irrep.dim();
    }
    return sdirs;
  }

  //*******************************************************************************************
  std::vector<Eigen::MatrixXd> SymGroupRep::irreducible_wedges(const SymGroup &head_group, std::vector<Index> &multiplicities)const {
    Index dim = get_MatrixXd(head_group[0])->cols();

    multivector<Eigen::VectorXd>::X<3> sdirs = calc_special_total_directions(head_group);
    std::vector<Eigen::MatrixXd> wedges(sdirs.size());
    double best_proj, tproj;
    multiplicities.clear();

    for(Index s = 0; s < sdirs.size(); s++) {
      wedges[s] = Eigen::MatrixXd::Zero(dim, sdirs[s].size());
      wedges[s].col(0) = sdirs[s][0][0];
      multiplicities.push_back(sdirs[s][0].size());
      for(Index i = 1; i < sdirs[s].size(); i++) {
        Index j_best = 0;
        best_proj = (wedges[s].transpose() * sdirs[s][i][0]).sum();
        for(Index j = 1; j < sdirs[s][i].size(); j++) {
          tproj = (wedges[s].transpose() * sdirs[s][i][j]).sum();
          if(tproj > best_proj) {
            best_proj = tproj;
            j_best = j;
          }
        }
        multiplicities.push_back(sdirs[s][i].size());
        wedges[s].col(i) = sdirs[s][i][j_best];
      }
    }
    return wedges;
  }
  //*******************************************************************************************
  ReturnArray<Array< Eigen::VectorXd> > SymGroupRep::_calc_special_irrep_directions(const SymGroup &head_group)const {
    if(!size() || !(get_MatrixXd(head_group[0]))) {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::_calc_special_irrep_directions() called on imcompatible SymGroupRep.\n Exiting...\n";
      exit(1);
    }

    Index i, j;

    Array<Array<Eigen::VectorXd> > sdirs; // Array of orbits of special directions

    Eigen::VectorXd tdir;

    const Array<Array<Array<Index> > > &isubs(head_group.subgroups());
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(*(get_MatrixXd(head_group[0])));

    Index dim = get_MatrixXd(head_group[0])->cols();

    QR.setThreshold(10 * TOL);

    // i loops over subgroup "orbits". There are 0 to 'dim' special directions associated with each orbit.
    Eigen::MatrixXd Reynolds(*(at(0)->get_MatrixXd()));
    for(i = 0; i < isubs.size(); i++) {
      if(isubs[i].size() == 0 || isubs[i][0].size() == 0) {
        std::cerr << "CRITICAL ERROR: In SymGroupRep::_calc_special_irrep_directions(), attempting to use a zero-size subgroup.\n Exiting...\n";
        exit(1);
      }
      Index s = 0;
      Reynolds.setZero();
      // j loops over elements of "prototype" subgroup to get the Reynold's operator for the vector
      // space on which the representation is defined
      for(j = 0; j < isubs[i][s].size(); j++) {
        Reynolds += *get_MatrixXd(head_group[isubs[i][s][j]]);
      }

      Reynolds /= double(isubs[i][s].size());

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

      // See if QR.matrixQ().col(0) has been found
      for(j = 0; j < sdirs.size(); j++) {
        if(sdirs[j].almost_contains(matrixQ.col(0))) {
          break;
        }
      }
      if(j == sdirs.size()) {

        //Get equivalents
        sdirs.push_back(Array<Eigen::VectorXd>());
        for(j = 0; j < head_group.size(); j++) {
          tdir = (*get_MatrixXd(head_group[j])) * matrixQ.col(0);
          if(!sdirs.back().almost_contains(tdir))
            sdirs.back().push_back(tdir);
        }
      }

      for(j = 0; j < sdirs.size(); j++) {
        if(sdirs[j].almost_contains(-matrixQ.col(0))) {
          break;
        }
      }
      if(j == sdirs.size() && dim > 1) {

        //Get equivalents
        sdirs.push_back(Array<Eigen::VectorXd>());
        for(j = 0; j < head_group.size(); j++) {
          tdir = -(*get_MatrixXd(head_group[j])) * matrixQ.col(0);
          if(!sdirs.back().almost_contains(tdir))
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

  ReturnArray<Array< Eigen::MatrixXd> > SymGroupRep::calc_special_subspaces(const SymGroup &head_group)const {
    if(!size() || !at(0)->get_MatrixXd()) {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::calc_special_subspaces() called on imcompatible SymGroupRep.\n Exiting...\n";
      exit(1);
    }

    Index i, j, k;

    Array<Array<Eigen::MatrixXd> > ssubs; // Array of orbits of special directions

    Eigen::MatrixXd tsub, ttrans;

    // Operation indices of subgroups
    const Array<Array<Array<Index> > > &sg_op_inds(head_group.subgroups());

    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> QR(*(at(0)->get_MatrixXd()));

    QR.setThreshold(TOL);

    Eigen::MatrixXd Reynolds(*(at(0)->get_MatrixXd()));

    // 'i' loops over subgroup "orbits". There are 0 to 'dim' special directions associated with each orbit.
    for(i = 0; i < sg_op_inds.size(); i++) {
      if(sg_op_inds[i].size() == 0 || sg_op_inds[i][0].size() == 0) {
        std::cerr << "CRITICAL ERROR: In SymGroupRep::calc_special_subspaces(), attempting to use a zero-size subgroup.\n Exiting...\n";
        exit(1);
      }

      Reynolds.setZero();
      // j loops over elements of "prototype" subgroup to get the Reynold's operator for the vector
      // space on which the representation is defined
      for(j = 0; j < sg_op_inds[i][0].size(); j++) {
        Reynolds += *(at(head_group[sg_op_inds[i][0][j]].index())->get_MatrixXd());
      }

      Reynolds /= double(sg_op_inds[i][0].size());

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
      ssubs.push_back(Array<Eigen::MatrixXd>());
      for(j = 0; j < head_group.size(); j++) {
        ttrans = (*(at(head_group[j].index())->get_MatrixXd())) * tsub;
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
          ssubs.swap_elem(i, j);
        if(ssubs[i][0].cols() == ssubs[j][0].cols() && ssubs[i].size() > ssubs[j].size())
          ssubs.swap_elem(i, j);
      }
    }

    return ssubs;
  }

  //*******************************************************************************************
  //assumes that representation is real-valued, and combines complex-valued irreps with their complex conjugate
  ReturnArray<Index> SymGroupRep::num_each_real_irrep(const SymGroup &head_group, bool verbose) const {
    const Array<bool > &complex_irrep(head_group.get_complex_irrep_list());
    Array<Index> tarray(num_each_irrep(head_group, verbose));

    if(tarray.size() != complex_irrep.size()) {
      std::cerr << "CRITICAL ERROR: Dimension mismatch in SymGroupRep::num_each_real_irrep. Exiting..\n";
      assert(0);
      exit(1);
    }

    // Go through and remove double counting of complex irreps
    for(Index i = 0; i < tarray.size(); i++) {
      if(tarray[i] && complex_irrep[i]) {
        if(i + 1 == tarray.size() || tarray[i] != tarray[i + 1]) {
          std::cerr << "CRITICAL ERROR: Invalid irrep decomposition found in SymGroupRep::num_each_real_irrep. Exiting..\n";
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
  //assumes that representation is real-valued
  ReturnArray<Index> SymGroupRep::num_each_irrep() const {
    if(!has_valid_master()) {
      std::cerr << "In SymGroupRep::num_each_irrep() attempting to find irrep decomposition of a SymGroupRep that has no valid master_group.\n";
      assert(0);
      exit(1);
    }
    return num_each_irrep(master_group());
  }

  //*******************************************************************************************

  ReturnArray<Index> SymGroupRep::num_each_irrep(const SymGroup &head_group, bool verbose) const {
    const Array<Array<Index> > &conj_class(head_group.get_conjugacy_classes());
    const Array<Array<std::complex<double> > > &char_table(head_group.character_table());

    Array<Index> tdec(conj_class.size(), 0);
    Array<double> repchar(conj_class.size(), 0);

    if(verbose) {
      std::cout << "Decomposing representation:\n ";
      for(Index i = 0; i < head_group.size(); i++) {
        std::cout << *(get_MatrixXd(head_group[i])) << "\n\n";
      }
    }

    //std::cout << "Character table is\n";
    for(Index i = 0; i < conj_class.size(); i++) {
      Index rep_index(head_group[conj_class[i][0]].index());
      repchar[i] = at(rep_index)->character();
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

  ReturnArray<SymGroupRepID> SymGroupRep::get_irrep_IDs(const SymGroup &head_group) const {
    Array<Index> irrep_decomp(num_each_real_irrep(head_group));
    Array<SymGroupRepID> irrep_IDs;
    for(Index i = 0; i < irrep_decomp.size(); i++) {
      if(irrep_decomp[i] == 0)
        continue;

      if(head_group.get_irrep_ID(i).empty())
        calc_new_irreps(head_group);

      if(head_group.get_irrep_ID(i).empty()) {
        std::cerr << "CRITICAL ERROR: In SymGroupRep::get_irrep_IDs(), cannot resolve irrep " << i << " which has multiplicity " << irrep_decomp[i] << ". Exiting...\n";
        assert(0);
        exit(1);
      }

      irrep_IDs.append(Array<SymGroupRepID>(irrep_decomp[i], head_group.get_irrep_ID(i)));
    }

    return irrep_IDs;

  }

  //*******************************************************************************************

  void SymGroupRep::calc_new_irreps(int max_iter) const {
    if(!has_valid_master()) {
      std::cerr << "In SymGroupRep::calc_new_irreps() attempting to calculate new irreps of a SymGroupRep that has no valid master_group.\n";
      assert(0);
      exit(1);
    }
    calc_new_irreps(master_group(), max_iter);
    return;
  }

  //*******************************************************************************************
  //
  // assumes that representation is real-valued. This means that complex-valued
  // irreps in the character table get combined into single 2-D irreps

  void SymGroupRep::calc_new_irreps(const SymGroup &head_group, int max_iter) const {

    if(!size() || !at(0)->get_MatrixXd()) {
      std::cerr << "WARNING:  In SymGroupRep::calc_new_irreps, size of representation is " << size() << " and MatrixXd address is " << at(0)->get_MatrixXd() << std::endl
                << "          No valid irreps will be returned.\n";
      return;
    }
    if(!head_group.size()) {
      assert(0 && "In SymGroupRep::calc_new_irreps(), passed an empty SymGroup!");
      return;
    }

    MasterSymGroup const *master_ptr(NULL);
    if(has_valid_master())
      master_ptr = &master_group();
    else if(head_group[0].has_valid_master())
      master_ptr = &head_group[0].master_group();
    else {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::calc_new_irreps(), there is no MasterSymGroup available!\n"
                << "                Exiting...\n";
    }
    int dim((at(0)->get_MatrixXd())->rows());

    Array<Index> i_decomp(num_each_real_irrep(head_group));
    int Nirrep(i_decomp.sum());
    if(Nirrep == 0) {
      std::cerr << "In SymGroupRep::calc_new_irreps(), there are no valid irreps, so no valid irreps will be returned.\n";
      return;
    }
    if(Nirrep == 1) {
      //std::cout << "There is only 1 irrep, so I'm returning a symmetrized copy of *this\n";
      head_group.set_irrep_ID(i_decomp.find(1), master_ptr->add_representation(coord_transformed_copy(_symmetrized_irrep_trans_mat(head_group))));
      return;
    }

    //std::cout << "This representation contains " <<  Nirrep << " irreps... About to build commuters!\n";
    Array<Eigen::MatrixXd> commuters;

    // Build 'commuters', which span space of matrices that commute with SymOpReps
    for(int i = 0; i < dim; i++) {
      for(int j = i; j < dim; j++) {
        Eigen::MatrixXd tmat(Eigen::MatrixXd::Zero(dim, dim)), tcommute(Eigen::MatrixXd::Zero(dim, dim));

        //We only build symmetric commuters, which assures we only select out real-valued representations
        tmat(i, j) = 1;
        tmat(j, i) = 1;
        //apply reynolds operator
        for(Index ns = 0; ns < head_group.size(); ns++) {
          tcommute += (*(get_MatrixXd(head_group[ns]))) * tmat * (*(get_MatrixXd(head_group[ns]))).transpose();
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
      std::cerr << "CRITICAL ERROR: In SymGroupRep::calc_new_irreps, reached maximum iterations (" << max_iter << ") without finding all irreps.\n"
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
      tmat = trans_mat * (*(at(head_group[i].index())->get_MatrixXd())) * trans_mat.transpose();
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
      SymGroupRep new_irrep(coord_transformed_copy(tmat));

      Array<Index> tdecomp(new_irrep.num_each_real_irrep(head_group));
      if(tdecomp.sum() != 1) {
        std::cerr << "CRITICAL ERROR: In SymGroupRep::calc_new_irreps(), representation identified as irreducible did not pass final tests. Exiting...\n";
        assert(0);
        exit(1);
      }
      Index i_irrep = tdecomp.find(1);
      if(!head_group.get_irrep_ID(i_irrep).empty()) {
        i1 = i2;
        continue;
      }

      SymGroupRepID irrep_ID = master_ptr->add_representation(new_irrep.coord_transformed_copy(new_irrep._symmetrized_irrep_trans_mat(head_group)));

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
    SymGroupRep new_irrep(coord_transformed_copy(tmat));

    //std::cout << "\nIts representation looks like this:\n";
    //new_irrep->print_MatrixXd(std::cout, head_group);
    //std::cout << '\n';
    Array<Index> tdecomp(new_irrep.num_each_real_irrep(head_group));
    if(tdecomp.sum() != 1) {
      std::cerr << "CRITICAL ERROR: In SymGroupRep::calc_new_irreps(), representation identified as irreducible did not pass final tests. Exiting...\n";
      assert(0);
      exit(1);
    }
    Index i_irrep = tdecomp.find(1);
    if(!head_group.get_irrep_ID(i_irrep).empty()) {
      return;
    }

    SymGroupRepID irrep_ID = master_ptr->add_representation(new_irrep.coord_transformed_copy(new_irrep._symmetrized_irrep_trans_mat(head_group)));

    //std::cout << "Setting the rep_ID of irrep " << i_irrep << " to " << irrep_ID;
    head_group.set_irrep_ID(i_irrep, irrep_ID);
    //std::cout << " and validating -- " << head_group.get_irrep_ID(i_irrep) << '\n';

    return;

  }

  //*******************************************************************************************

  bool SymGroupRep::is_irrep() const {
    double tvalue = 0;

    for(Index i = 0; i < size(); i++) {
      tvalue += (at(i)->character()) * (at(i)->character());
    }

    if(Index(round(tvalue)) == master_group().size()) {
      return true;
    }
    else {
      return false;
    }
  }

  //*******************************************************************************************

  bool SymGroupRep::is_irrep(const SymGroup &head_group) const {
    double tvalue = 0;

    for(Index i = 0; i < head_group.size(); i++) {
      tvalue += (at(head_group[i].index())->character()) * (at(head_group[i].index())->character());
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
  Eigen::MatrixXd SymGroupRep::get_irrep_trans_mat(const SymGroup &head_group) const {
    std::vector< std::vector<Index> > subspaces;
    return get_irrep_trans_mat(head_group, subspaces);
  }

  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  // Also populate 'subspaces' with lists of columns that form irreps
  Eigen::MatrixXd SymGroupRep::get_irrep_trans_mat(const SymGroup &head_group, std::vector< std::vector<Index> > &subspaces) const {
    //std::cout << "INSIDE GET_IRREP_TRANS_MAT\n";
    const Array<bool > &complex_irrep(head_group.get_complex_irrep_list());
    //const Array<Array<std::complex<double> > > &char_table(master_group().character_table());
    const Array<Array<std::complex<double> > > &char_table(head_group.character_table());

    Array<Index> irrep_decomp(num_each_real_irrep(head_group));

    //get_irrep_IDs() harvests any new irreps that appear in this representation
    Array<SymGroupRepID> irrep_IDs(get_irrep_IDs(head_group));

    int rep_dim = get_MatrixXd(head_group[0])->rows();

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
      assert(irrep_dim == (t_irrep.get_MatrixXd(head_group[0])->rows()));

      //get the (0,0) projection operator -> projects out vectors that transform as the first basis vector of irrep 'i'
      tproj.setZero();
      //std::cout << "PROJECTING OUT THE FOLLOWING IRREP\n";
      for(Index j = 0; j < head_group.size(); j++) {
        //std::cout << "OP " << j << ":\n";

        tproj += (*(t_irrep.get_MatrixXd(head_group[j])))(0, 0) * (*get_MatrixXd(head_group[j]));
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
        std::cerr << "CRITICAL ERROR: In SymGroupRep::get_irrep_trans_mat(), did not find correct number of basis vectors to account for\n"
                  << "                multiplicity of irrep #" << i << ". There should be " << irrep_decomp[i] << " but I found " << QR.rank() << "\n"
                  << "                Exiting...\n";
        exit(1);
      }

      //std::cout << "QR yields:\n " << QR.matrixQ() << "\n\n";
      for(EigenIndex n = 0; n < QR.rank(); n++) {
        subspaces.push_back(std::vector<Index>(irrep_dim));
        std::iota(subspaces.back().begin(), subspaces.back().end(), col_count + n * irrep_dim);
        trans_mat.col(col_count + n * irrep_dim) = QR.matrixQ().col(n);
      }
      //std::cout << "Initial trans_mat is:\n" << trans_mat << "\n\n";

      // NOW, let's make the (0,j) projection matrices to project out the other basis functions
      for(int j = 1; j < irrep_dim; j++) {

        tproj.setZero();
        for(Index k = 0; k < head_group.size(); k++) {
          //tproj += (*(t_irrep->get_MatrixXd(head_group[k])))(0,j) * (*get_MatrixXd(head_group[k])); // <-- this is wrong, not totally sure why
          tproj += (*(t_irrep.get_MatrixXd(head_group[k])))(j, 0) * (*get_MatrixXd(head_group[k]));
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

    return trans_mat.transpose();
  }
  //*******************************************************************************************

  // Finds the transformation matrix that block-diagonalizes this representation into irrep blocks
  // The ROWS of trans_mat are the new basis vectors in terms of the old such that
  // new_symrep_matrix = trans_mat * old_symrep_matrix * trans_mat.transpose();
  Eigen::MatrixXd SymGroupRep::get_irrep_trans_mat_blind(const SymGroup &head_group) const {

    if(!size() || !head_group.size() || !get_MatrixXd(head_group[0].index())) {
      std::cerr << "WARNING:  In SymGroupRep::calc_new_irreps, size of representation is " << size() << " and MatrixXd address is " << get_MatrixXd(head_group[0].index()) << std::endl
                << "          No valid irreps will be returned.\n";
      return Eigen::MatrixXd();
    }

    int dim(get_MatrixXd(head_group[0].index())->rows());
    Array<Eigen::MatrixXcd> commuters;
    std::cout.precision(8);
    std::cout.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    // Identity always commutes, and by putting it first we prevent some accidental degeneracies
    commuters.push_back(Eigen::MatrixXcd::Identity(dim, dim) / sqrt(double(dim)));
    typedef std::complex<double> cplx;
    Array<cplx> phase;
    phase.push_back(cplx(1.0, 0.0)); // 1+0i
    phase.push_back(cplx(0.0, 1.0)); // 0+1i

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esolve;
    Eigen::HouseholderQR<Eigen::MatrixXd> qr;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> colqr;
    colqr.setThreshold(0.001);
    int Nfound(0);
    Eigen::MatrixXcd tmat(Eigen::MatrixXcd::Zero(dim, dim)), tcommute(Eigen::MatrixXcd::Zero(dim, dim));
    Eigen::MatrixXd trans_mat(Eigen::MatrixXd::Zero(dim, dim)), kernel(Eigen::MatrixXd::Random(dim, dim).householderQr().householderQ());



    //std::cout << "Commuter are:\n";
    // Build 'commuters', which span space of real matrices that commute with SymOpReps
    // 'kernel' is the kernel of trans_mat, and as the loop progresses, 'kernel' shrinks and the rank for trans_mat increases
    Array<Index> irrep_dims;
    std::cout << "~~~~~~~Dimension is " << dim << "\n";
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
            tcommute += (*(get_MatrixXd(head_group[ns].index()))) * tmat * (*(get_MatrixXd(head_group[ns].index()))).transpose();
            //tcommute += phase[nph]*(get_MatrixXd(ns)->col(i)) * (get_MatrixXd(ns)->row(j))+std::conj(phase[nph])*(get_MatrixXd(ns)->col(j)) * (get_MatrixXd(ns)->row(i));
          }

          //Do Gram-Shmidt while building 'commuters'
          for(Index nc = 0; nc < commuters.size(); nc++) {
            double tproj((commuters[nc].array().conjugate()*tcommute.array()).sum().real()); //Frobenius product
            tcommute -= tproj * commuters[nc];
          }
          double tnorm(tcommute.norm());
          if(tnorm > TOL) {
            commuters.push_back(tcommute / tnorm);
            //std::cout << commuters.back() << "\n\n";
          }
          else continue;  // Attempt to construct the next commuter...

          //Finished building commuter now we can try to harvest irreps from it

          // construct trans_mat from the non-degenerate irreps obtained from the commuter

          Index nc = commuters.size() - 1;

          // magnify the range of eigenvalues to be (I think) independent of matrix dimension by multiplying by dim^{3/2}
          esolve.compute(double(dim)*sqrt(double(dim))*kernel.transpose()*commuters[nc]*kernel);
          //std::cout << "KERNEL IS:\n" << kernel << "\n\n";
          //std::cout << esolve.eigenvalues().size() << " EIGENVALUES of commuter " << nc << " are \n" << esolve.eigenvalues().transpose() << "\n";
          Array<Index> subspace_dims = partition_distinct_values(esolve.eigenvalues());
          //std::cout << "Eigenvectors are:\n" << esolve.eigenvectors().real() << "\n\n" << esolve.eigenvectors().imag() << "\n\n";
          //std::cout << "QR Eigenvectors are:\n"
          //<< Eigen::MatrixXcd(esolve.eigenvectors().householderQr().householderQ()).real() << "\n\n"
          //<< Eigen::MatrixXcd(esolve.eigenvectors().householderQr().householderQ()).imag() << "\n\n";

          tmat = kernel * (esolve.eigenvectors().householderQr().householderQ());
          //std::cout << "Considering vectors:\n" << tmat.real() << "\n\n" << tmat.imag() << "\n\n";
          //std::cout << "orthonormality check 1:\n" << (tmat.adjoint()*tmat).real() << "\n\n" << (tmat.adjoint()*tmat).imag() << "\n\n";
          //std::cout << "orthonormality check 2:\n" << (tmat * tmat.adjoint()).real() << "\n\n" << (tmat * tmat.adjoint()).imag() << "\n\n";

          // make transformed copy of the representation
          Array<Eigen::MatrixXcd> trans_rep(head_group.size());
          Eigen::MatrixXd block_shape(Eigen::MatrixXd::Zero(kernel.cols(), kernel.cols()));
          //std::cout << "Transformed representation is:\n";
          for(Index i = 0; i < head_group.size(); i++) {
            trans_rep[head_group[i].index()] = tmat.adjoint() * (*get_MatrixXd(head_group[i].index())) * tmat;
            block_shape += (trans_rep[head_group[i].index()].cwiseProduct(trans_rep[head_group[i].index()].conjugate())).real();
            //std::cout << trans_rep[i] << "\n\n";
          }

          Index last_i = 0;
          for(Index ns = 0; ns < subspace_dims.size(); ns++) {
            double tnorm(0);
            Array<std::complex<double> > char_array;

            for(Index ng = 0; ng < trans_rep.size(); ng++) {
              cplx tchar(0, 0);
              for(Index i = last_i; i < last_i + subspace_dims[ns]; i++)
                tchar += trans_rep[ng](i, i);
              char_array.push_back(tchar);
              tnorm += std::norm(tchar);
            }

            if(almost_equal(tnorm, double(head_group.size()))) { // this representation is irreducible
              std::cout << "THE REPRESENTATION IS IRREDUCIBLE!!\n";
              std::cout << "It's characters are: " << char_array << "\n\n";
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
                irrep_dims.push_back(subspace_dims[ns]);
                if(rnk == 2 * subspace_dims[ns])
                  irrep_dims.push_back(subspace_dims[ns]);
                ttrans_mat = Eigen::MatrixXd(qr.householderQ()).leftCols(rnk);
                SymGroupRep t_rep(coord_transformed_copy(ttrans_mat.transpose()));
                //std::cout << "***Adding columns!\n" << Eigen::MatrixXd(qr.householderQ()).leftCols(rnk) << "\n\n";
                //trans_mat.block(0, Nfound, dim, rnk) = Eigen::MatrixXd(qr.householderQ()).leftCols(rnk);
                trans_mat.block(0, Nfound, dim, rnk) = ttrans_mat * (t_rep._symmetrized_irrep_trans_mat(head_group)).transpose();
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

    // make transformed copy of the representation
    Eigen::MatrixXd block_shape(Eigen::MatrixXd::Zero(dim, dim));
    //std::cout << "Transformed representation is:\n";
    for(Index i = 0; i < head_group.size(); i++) {
      block_shape += (trans_mat.transpose() * (*get_MatrixXd(head_group[i].index())) * trans_mat).cwiseAbs2();
    }
    std::cout << "BLOCK MATRIX IS:\n"
              << block_shape << "\n\n";
    std::cout << "IRREP DIMENSIONS ARE: " << irrep_dims << "\n\n";
    return trans_mat.transpose();

  }

  //*******************************************************************************************

  ReturnArray<Eigen::MatrixXd> SymGroupRep::get_projection_operators() const {

    const Array<Array<Index> > &conj_class(master_group().get_conjugacy_classes());
    const Array<Array<std::complex<double> > > &char_table(master_group().character_table());

    Eigen::MatrixXd tmat((*(at(0)->get_MatrixXd())).rows(), (*(at(0)->get_MatrixXd())).cols());
    tmat.setZero();
    Array<Eigen::MatrixXd> tarray(conj_class.size());

    for(Index i = 0; i < conj_class.size(); i++) {
      double dimension = char_table[i][0].real();
      for(Index j = 0; j < master_group().size(); j++) {
        double character =  char_table[i][master_group().class_of_op(j)].real();
        Eigen::MatrixXd rep_mat = *(at(j)->get_MatrixXd());
        tmat += (character * rep_mat);
      }
      tarray[i] = ((dimension) / (master_group().size())) * tmat;
    }

    return tarray;
  }

  //*******************************************************************************************

  jsonParser &SymGroupRep::to_json(jsonParser &json) const {
    json.put_obj();

    // Member not included in json:
    //
    //   Pointer to the m_master_group that generated this SymGroupRep
    //   MasterSymGroup const *m_master_group;

    // class SymGroupRep : public Array<SymOpRepresentation *>
    json["symop_representations"].put_array();
    for(Index i = 0; i < size(); i++) {
      at(i)->to_json(json["symop_representations"]);
    }

    // int m_rep_ID;
    json["m_rep_ID"] = m_rep_ID;

    return json;
  }

  //*******************************************************************************************

  SymGroupRep subset_permutation_rep(const SymGroupRep &permute_rep, const Array<Index>::X2 &subsets) {
    SymGroupRep new_rep(SymGroupRep::NO_HOME, permute_rep.size());
    if(permute_rep.has_valid_master())
      new_rep = SymGroupRep(permute_rep.master_group());

    Array<Index> perm_sub;
    if(subsets.size())
      perm_sub.resize(subsets[0].size());

    Array<Index> tperm(subsets.size());
    for(Index np = 0; np < permute_rep.size(); ++np) {
      Permutation const *perm_ptr(permute_rep.get_permutation(np));
      if(perm_ptr == NULL) {
        // This check may not make sense. Some representations may not have all operations.  If it's causing problems, comment it out
        std::cerr << "CRITICAL ERROR: In subset_permutation_rep(SymGroupRep,Array<Index>::X2), permute_rep is missing at least one permutation!\n"
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
          if(perm_sub.all_in(subsets[ns2])) {
            tperm[ns] = ns2;
            break;
          }
        }
        if(ns2 >= subsets.size()) {
          std::cerr << "CRITICAL ERROR: In subset_permutation_rep(SymGroupRep,Array<Index>::X2), subsets are not closed under permute_rep permutations!\n"
                    << "                Exiting...\n";
          assert(0);
          exit(1);
        }
      }
      new_rep.set_rep(np, SymPermutation(tperm));
    }
    return new_rep;
  }

  //*******************************************************************************************

  SymGroupRep permuted_direct_sum_rep(const SymGroupRep &permute_rep, const Array<SymGroupRep const *> &sum_reps) {
    SymGroupRep new_rep(SymGroupRep::NO_HOME, permute_rep.size());
    if(permute_rep.has_valid_master())
      new_rep = SymGroupRep(permute_rep.master_group());

    Array<Index> rep_dims(sum_reps.size(), 0);
    if(permute_rep.get_permutation(0) == NULL) {
      std::cerr << "CRITICAL ERROR: In permuted_direct_sum_rep(SymGroupRep,Array<SymGroupRep const*>), permute_rep does not contain permutations!\n"
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
         || sum_reps[i]->get_MatrixXd(0) == NULL) {
        std::cerr << "CRITICAL ERROR: In permuted_direct_sum_rep(SymGroupRep,Array<SymGroupRep const*>), found incompatible SymGroupReps!\n"
                  << "                Exiting...\n";
        assert(0);
        exit(1);
      }
      rep_dims[i] = (*sum_reps[i])[0]->get_MatrixXd()->cols();
    }
    //std::cout << "SUM_REPS is " << sum_reps << "\n";
    //std::cout << "REP_DIMS is " << rep_dims << "\n";
    Array<Index> sum_inds(cum_sum(rep_dims));
    Eigen::MatrixXd sum_mat(rep_dims.sum(), rep_dims.sum());
    Eigen::MatrixXd const *rep_mat_ptr(NULL);
    Permutation const *perm_ptr;
    for(Index g = 0; g < permute_rep.size(); g++) {
      sum_mat.setZero();
      perm_ptr = permute_rep.get_permutation(g);
      for(Index r = 0; r < perm_ptr->size(); r++) {
        if(rep_dims[r] == 0) continue;
        Index new_r = (*perm_ptr)[r];
        rep_mat_ptr = sum_reps[new_r]->get_MatrixXd(g);
        sum_mat.block(sum_inds[r], sum_inds[new_r], rep_dims[new_r], rep_dims[new_r]) = *rep_mat_ptr;
      }
      new_rep.set_rep(g, SymMatrixXd(sum_mat));
    }
    return new_rep;
  }

  //*******************************************************************************************
  SymGroupRep kron_rep(const SymGroupRep &LHS, const SymGroupRep &RHS) {
    if((LHS.has_valid_master() != RHS.has_valid_master()) || (LHS.has_valid_master() && &LHS.master_group() != &RHS.master_group())) {
      std::cerr << "CRITICAL ERROR: In kron_rep(SymGroupRep,SymGroupRep), taking product of two incompatible SymGroupReps!\n"
                << "                Exiting...\n";
      exit(1);
    }
    SymGroupRep new_rep(SymGroupRep::NO_HOME, LHS.size());
    if(LHS.has_valid_master())
      new_rep = SymGroupRep(LHS.master_group());

    Eigen::MatrixXd tmat;
    for(Index i = 0; i < LHS.size(); i++) {
      //std::cout << "rep1 matrix:\n";
      //std::cout << *(rep1->get_MatrixXd(i)) << '\n';
      //std::cout << "rep2 matrix:\n";
      //std::cout << *(rep2->get_MatrixXd(i)) << '\n';
      assert(LHS.get_MatrixXd(i) && RHS.get_MatrixXd(i));
      kroneckerProduct(*(LHS.get_MatrixXd(i)), *(RHS.get_MatrixXd(i)), tmat);
      //std::cout << "Total matrix:\n" << tmat << '\n';
      new_rep.set_rep(i, SymMatrixXd(tmat));
    }
    return new_rep;
  }


  //*******************************************************************************************

  // If 'm_home_group' is not nullptr, should be initialized accordingly
  void SymGroupRep::from_json(const jsonParser &json) {
    // Member not included in json:
    //
    //   Pointer to the m_master_group that generated this SymGroupRep
    //   MasterSymGroup const *m_master_group;

    // class SymGroupRep : public Array<SymOpRepresentation *>

    for(Index i = 0; i < size(); i++) {
      delete at(i);
    }
    clear();
    //std::cout << "Resizing SymGroupRep to " << json["symop_representations"].size() << std::endl;
    resize(json["symop_representations"].size());
    //std::cout << "Reading in the symmetry operations" << std::endl;
    for(int i = 0; i < json["symop_representations"].size(); i++) {
      /// This allocates a new object to 'at(i)'.
      CASM::from_json(at(i), json["symop_representations"][i]);
    }

    //std::cout << "Reading in m_rep_id" << std::endl;
    // int m_rep_ID;
    CASM::from_json(m_rep_ID, json["m_rep_ID"]);

    //std::cout << "Done reading in the permute_group" << std::endl;
  }

  //*******************************************************************************************

  jsonParser &to_json(const SymGroupRep &rep, jsonParser &json) {
    return rep.to_json(json);
  }


  // If 'm_master_group' is not NULL, should be initialized accordingly
  void from_json(SymGroupRep &rep, const jsonParser &json) {
    rep.from_json(json);
  }


}

