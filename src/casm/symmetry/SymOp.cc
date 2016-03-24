#include "casm/symmetry/SymOp.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"

namespace CASM {

  //SymOpRepresentation::~SymOpRepresentation() {};

  const Eigen::Matrix3d &SymOp::matrix() const {
    return m_mat;
  }

  //**********************************************************

  const Vector3<double> &SymOp::eigenvec() const {
    return m_eigenvec;
  }


  //**********************************************************

  const double &SymOp::map_error() const {
    return m_map_error;
  }

  //**********************************************************

  void SymOp::set_map_error(const double &value) {
    m_map_error = value;
    return;
  }


  //**********************************************************

  const Coordinate &SymOp::tau() const {
    return m_tau;
  }

  //**********************************************************

  double SymOp::tau(int i) const {
    return m_tau[i];
  }


  //**********************************************************

  void SymOp::set_index(const MasterSymGroup &new_group, Index new_index) {
    if((valid_index(new_index) && new_index < new_group.size())
       && (this == &(new_group[new_index]) ||
           matrix().is_equal(new_group[new_index].matrix()))) {
      m_master_group = &new_group;
      op_index = new_index;
    }
    else {
      m_master_group = NULL;
      op_index = -1;
    }
  }

  //**********************************************************
  void SymOp::set_rep(Index rep_ID, const SymOpRepresentation &op_rep) const {
    SymGroupRep const *tRep(master_group().representation(rep_ID));
    if(!tRep) {
      std::cerr << "CRITICAL ERROR: In SymOp::set_matrix_rep(" << rep_ID << "), representation was not found.\n"
                << "                Exiting...\n";
      exit(1);
    }

    tRep->set_rep(op_index, op_rep);

    return;
  }

  //**********************************************************
  Eigen::MatrixXd const *SymOp::get_matrix_rep(Index rep_ID) const {
    SymGroupRep const *tRep(master_group().representation(rep_ID));
    if(!tRep) return NULL;

    return (tRep->at(op_index))->get_MatrixXd();
  }

  //**********************************************************
  SymBasisPermute const *SymOp::get_basis_permute_rep(Index rep_ID) const {

    SymGroupRep const *tRep(master_group().representation(rep_ID));
    if(!tRep) {
      std::cerr << "Warning: You have requested information from a nonexistent representation!\n"
                << "m_master_group pointer is " << m_master_group << '\n';
      return NULL;
    }

    return (tRep->at(op_index))->get_ucc_permutation();
  }

  //**********************************************************
  Permutation const *SymOp::get_permutation_rep(Index rep_ID) const {
    SymGroupRep const *tRep(master_group().representation(rep_ID));
    if(!tRep) return NULL;

    return (tRep->at(op_index))->get_permutation();
  }

  //**********************************************************
  Array<Eigen::MatrixXd const * > SymOp::get_matrix_reps(Array<Index> rep_IDs) const {
    Array<Eigen::MatrixXd const * > tmat;
    for(Index i = 0; i < rep_IDs.size(); i++) {
      tmat.push_back(get_matrix_rep(rep_IDs[i]));

    }
    return tmat;
  }

  //**********************************************************

  SymOp SymOp::operator*(const SymOp &RHS) const {
    SymOp t_op(matrix() * RHS.matrix(),
               tau() + matrix() * RHS.tau(),
               *home);

    if(m_master_group && (m_master_group == RHS.m_master_group)) {
      t_op.set_index(master_group(), master_group().ind_prod(index(), RHS.index()));
    }
    else if(RHS.m_master_group && !m_master_group && symmetry == identity_op) {
      t_op.set_index(*RHS.m_master_group, RHS.index());
    }
    else if(m_master_group && !RHS.m_master_group && RHS.symmetry == identity_op) {
      t_op.set_index(master_group(), index());
    }
    else if(symmetry == identity_op && RHS.symmetry == identity_op) {
      t_op.symmetry = identity_op;
      //t_op.op_index = 0;
    }
    //else{
    //std::cout << "This symmetry is " << symmetry << " with head " << m_master_group << " and RHS symmetry is " << RHS.symmetry << " with head " << RHS.m_master_group << "\n";
    //}

    return t_op;

  }

  //**********************************************************

  SymOp &SymOp::operator+=(const Coordinate &RHS) {
    tau_vec(CART) += RHS(CART) - matrix() * RHS(CART);
    location(CART) += RHS(CART);
    return (*this);
  }

  //**********************************************************

  SymOp &SymOp::operator-=(const Coordinate &RHS) {
    tau_vec(CART) -= RHS(CART) - matrix() * RHS(CART);
    location(CART) -= RHS(CART);
    return (*this);
  }

  //**********************************************************
  // SymOp matrix is unitary, so inverse is equivalent to transpose.
  // To do inverse of translation, you must perform
  // inverse matrix operaton on translation and subtract
  SymOp SymOp::inverse() const {
    SymOp t_op(matrix().transpose(),
               -(matrix().transpose() * tau(CART)),
               *home,
               CART);
    if(m_master_group) {
      t_op.set_index(master_group(), master_group().ind_inverse(index()));
    }
    else if(is_identity()) {
      t_op.op_index = 0;
    }

    return t_op;
  }

  //**********************************************************

  SymOp SymOp::no_trans() const {

    SymOp t_op(*this);

    t_op.tau_vec(CART) = Vector3<double>(0, 0, 0);

    return t_op;
  };

  //**********************************************************

  void SymOp::within() {
    tau_vec.within();
  }

  //**********************************************************

  bool SymOp::compare(const SymOp &RHS, double eq_tol) const {
    return calc(CART) &&
           RHS.calc(CART) &&
           matrix().is_equal(RHS.matrix(), eq_tol) &&
           tau_vec.min_dist(RHS.tau_vec) < eq_tol;

  }

  //**********************************************************

  bool SymOp::operator==(const SymOp &RHS) const {
    return calc(CART) &&
           RHS.calc(CART) &&
           matrix().is_equal(RHS.matrix()) &&
           tau(CART).is_equal(RHS.tau(CART));
  };

  //******************************************************
  //
  // Functions to check symmetry type:
  //
  //******************************************************

  bool SymOp::is_identity() const {
    if(symmetry == invalid_op)
      get_sym_type();
    return symmetry == identity_op;
  };

  //******************************************************
  /**
   * Find the invariant point after applying symmetry

   *	- General equation

   * P is the invariant point
   * P(CART) = Sym_op*P(CART) + tau(CART)
   * P(CART) = (I-Sym_op).inverse()*tau(CART)

   *	- The General equation is only valid for point sysmetry

   *	- Plane symmetry ( Mirror or Glide)
   *  SP = -P (always sture at the origin)
   *  P(CART) = -P(CART) + tau_per (perpendicular to eigenvec)
   *  P = tau_per/2

   *	- aixes sysmmetry ( Rotation or Screw)
   *   change the coordinate system having eigen vector as a z coordinate,
   *	 so, it makes 3D rotation to 2D roation
   *   find M coordinate transfer matrix
   *   M = ( tau_pp ; eigen x tau_pp ; eigen)
   *	 define Inew = (1 0 0; 0 1 0; 0 0 0)
   *   Pnew = (Inew-Snew)^-1*tau_pp(new)
   *   P  = (MInewM.inverse() - S).inverse*tau_pp
   */

  void SymOp::find_location() const {

    location(CART) = Vector3<double>(0, 0, 0);
    Eigen::Matrix3d tMat, inv_tMat;
    Vector3<double> tau_pp, tau_ll;

    if(type() == invalid_op) {
      get_sym_type();
    }


    if(type() == identity_op) {
      //std::cout<<"all points are high symmetry  points \n";
      return;
    }


    if((type() == rotoinversion_op) || (type() == inversion_op)) {
      //  std::cout << "eigenvec norm is " << eigenvec(CART).norm() << '\n';

      tMat = Eigen::Matrix3d::identity() - m_mat;
      inv_tMat = tMat.inverse();
      location(CART) = inv_tMat * tau(CART);

      return;
    }
    if((type() == mirror_op) || (type() == glide_op)) {
      //component of tau parallel to eigenvector:
      tau_ll = (eigenvec(CART).dot(tau(CART)) / eigenvec(CART).dot(eigenvec(CART))) * eigenvec(CART);

      location(CART) = tau_ll / 2.0;
      return;
    }
    if(((type() == rotation_op) || (type() == screw_op))) {


      // std::cout << "eigenvec norm is " << eigenvec(CART).norm() << '\n';

      //component of tau parallel to rotation axis
      tau_ll = (eigenvec(CART).dot(tau(CART)) / eigenvec(CART).dot(eigenvec(CART))) * eigenvec(CART);

      //component in the plane of rotation
      tau_pp = tau(CART) - tau_ll;

      if(tau_pp.is_zero()) {
        //rotation axis passes through origin
        //std::cout << location() << "\n";
        return;
      }


      Vector3<double> X, Y, Z, tY;
      Eigen::Matrix3d M, I_new(Eigen::Matrix3d::identity());

      tY = eigenvec(CART).cross(tau_pp);
      Y = tY / tY.norm();
      X = tau_pp / tau_pp.norm();
      Z = eigenvec(CART) / eigenvec(CART).norm();
      M(0, 0) = X.at(0);
      M(1, 0) = X.at(1);
      M(2, 0) = X.at(2);
      M(0, 1) = Y.at(0);
      M(1, 1) = Y.at(1);
      M(2, 1) = Y.at(2);
      M(0, 2) = Z.at(0);
      M(1, 2) = Z.at(1);
      M(2, 2) = Z.at(2);

      I_new(2, 2) = 0;

      inv_tMat = M * I_new * M.inverse() - m_mat;

      location(CART) = inv_tMat.inverse() * tau_pp;

      return;
    }

    std::cerr << "DISASTER in SymOp::find_location!!\n Attempted to find symmetr location, but symmetry type is invalid!\n";
  }

  //******************************************************

  bool SymOp::is_mirror() const {
    if(symmetry == invalid_op)
      get_sym_type();
    return symmetry == mirror_op;
  }

  //******************************************************

  bool SymOp::is_glide() const {
    if(symmetry == invalid_op)
      get_sym_type();
    return symmetry == glide_op;
  }

  //******************************************************

  bool SymOp::is_rotation() const {
    if(symmetry == invalid_op)
      get_sym_type();
    return symmetry == rotation_op;
  }

  //******************************************************

  bool SymOp::is_screw() const {
    if(symmetry == invalid_op)
      get_sym_type();
    return symmetry == screw_op;
  }

  //******************************************************

  bool SymOp::is_inversion() const {
    if(symmetry == invalid_op)
      get_sym_type();
    return symmetry == inversion_op;
  }

  //******************************************************

  bool SymOp::is_rotoinversion() const {
    if(symmetry == invalid_op)
      get_sym_type();
    return symmetry == rotoinversion_op;
  }

  //******************************************************

  bool SymOp::is_invalid() const {
    if(symmetry == invalid_op)
      get_sym_type();
    return symmetry == invalid_op;
  }



  //******************************************************

  SymOp &SymOp::apply_sym(const SymOp &op) {
    (*this) = op * (*this) * (op.inverse());
    //get_sym_type(); // why would we do this?
    return *this;
  }

  //******************************************************
  /**
   * Function tests a symmetry matrix and assigns symmetry type
   * to variable symmetry.
   *
   * The function checks for (in this order) identity, inversion
   * mirror, glide, rotoinversion, rotation, and screw.
   */
  //******************************************************
  void SymOp::get_sym_type() const {

    double det = 0.0;
    double trace = 0.0;

    // temporary matrix
    Eigen::Matrix3d tmat;

    //Check to see if symmetry_mat exists
    // if does not exist, update
    if(!calc(CART)) {
      std::cerr << "WARNING: Inside SymOp::get_sym_type(), can't get symmetry type because Cartesian transformation does not exist.\n";
      return;
    }


    det = matrix().determinant();
    trace = matrix().trace();

    //Check for IDENTITY
    //Trace = 3
    if(std::abs(trace - 3.0) < TOL) {
      symmetry = identity_op;
      return;
    }


    //Check for INVERSION
    //Trace = -3
    if(std::abs(trace + 3.0) < TOL) {
      symmetry = inversion_op;
      return;
    }

    //Check for MIRROR
    //Trace = 1; det = -1
    if((det < (-TOL)) && (std::abs(trace - 1) < TOL)) {
      mirror_check();
      return;
    }

    //Check for ROTATION
    //Check for ROTOINVERSION
    //Rotation symmetry matrix det = 1 (does not change the orientation)
    //Rotation has trace 1+2cos(theta)
    //Rotoinversion det = -1
    //Extract just rotation matrix
    else {
      tmat = matrix();

      if(det < (-TOL)) {  //then rotoinversion
        symmetry = rotoinversion_op;

        //"undo" the inversion part of the trace and det
        //to obtain trace, determinant of rotation matrix
        // and rotation matrix itself.
        trace = -trace;
        det = -det;
        tmat *= -1;
        calc_rotation_angle(tmat, det, trace);

      }
      else { //rotation
        symmetry = rotation_op;
        calc_rotation_angle(tmat, det, trace);

        screw_check();
        return;
      }


    }//end of rotation/rotoinversion check

    //If the symmetry matrix passes none of the tests
    //symmetry remains as invalid_op
  }

  //******************************************************
  /**
   * Checks for glide symmetry operation
   *
   * Function requires that a mirror symmetry operation
   * was found first before calculating the translational
   * component (screw_glide_shift) of the glide. It takes
   * tau_vec, subtracts off the component perpendicular
   * to the mirror plane (parallel to eigenvector), leaving
   * just the component parallel to the mirror plane.
   */
  //******************************************************


  void SymOp::glide_check() const {

    Coordinate tcoord(get_home());

    // Check for GLIDE
    // Mirror + Translation over half lattice vector // to mirror plane
    if(symmetry != mirror_op) return;

    //Find magnitude of translation parallel to sym_op eigenvector (normal to mirror plane)
    double tnorm = tau_vec(CART).dot(eigenvec(CART));

    if(std::abs(tnorm - tau_vec(CART).length()) < TOL) {
      // vector tau and eigenvector parallel -- therefore not glide
      return;
    }

    //subtract off component of translation that is parallel
    //to eigenvector to obtain purely translational component
    tcoord(CART) = tau_vec(CART) - eigenvec(CART) * tnorm;

    // If tcoord is not a lattice translation
    if(!tcoord.is_lattice_shift()) {
      screw_glide_shift(CART) = tcoord(CART);
      symmetry = glide_op;
    } // But if it is a lattice shift and the coord is not zero
    //     else if(!tcoord(CART).is_zero(TOL)) { //Added by Ivy 01/24/13
    //       //std::cout << "LATTICE_GLIDE \n";
    //       screw_glide_shift(CART) = tcoord(CART);
    //       symmetry = lattice_glide_op;
    //     }

    // otherwise, symmetry remains a mirror_op
    // when exiting this routine
    return;

  }



  //******************************************************
  /**
   * Checks for screw symmetry operation
   *
   * Function requires that a rotation symmetry operation
   * was found first before checking for the eigenvector
   * of the screw symmetry operation.
   */
  //******************************************************
  void SymOp::screw_check() const {

    Coordinate tcoord(get_home());

    //Rotation + Translation parallel to rotation axis
    if(symmetry != rotation_op) return;

    //Find magnitude of translation parallel to sym_op eigenvector (rotation axis)
    //i.e. dot product of the two vectors, normalized
    double tnorm = tau_vec(CART).dot(eigenvec(CART));

    if(std::abs(tnorm) < TOL) {
      // vector tau and eigenvector perpendicular -- therefore not screw
      return;
    }

    //Get component of translation parallel to eigenvector
    //i.e. dotproduct(tau, eigenvec)*eigenvec
    tcoord(CART) = eigenvec(CART) * tnorm;

    // If tcoord is not a lattice translation
    if(!tcoord.is_lattice_shift()) {
      screw_glide_shift(CART) = tcoord(CART);
      symmetry = screw_op;
    }// But if it is a lattice shift and the coord is not zero
    //     else if(!tcoord(CART).is_zero(TOL)) { //Added by Ivy 01/24/13
    //       //std::cout << "LATTICE_SCREW \n";
    //       screw_glide_shift(CART) = tcoord(CART);
    //       symmetry = lattice_screw_op;
    //     }

    // otherwise symmetry remains rotation_op
    // when exiting this routine
    return;

  }

  //******************************************************
  /**
   * Checks for mirror symmetry operation.
   */
  //******************************************************
  void SymOp::mirror_check() const {
    int i;

    //Mirror = 180 degree rotation + inversion
    //Check for 1) det < 0 (inversion) and 2) trace = 1 (180 rotation)
    //Get normal vector of mirror plane

    //Turn CART form of symmetry matrix into Eigen class Matrix
    //Eigen::Matrix<double, 3, 3> tmat = (matrix());
    Vector3< Vector3<std::complex<double> > > all_eigenvectors;
    Vector3< std::complex<double> > all_eigenvalues;
    all_eigenvalues = matrix().eigen(all_eigenvectors);


    //Want the eigenvector that is perpendicular to the mirror plane.
    //i.e. the eigenvector corresponding to the -1 eigenvalue.
    //A mirror symmetry operation has 3 eigenvalues: 1, 1, -1
    for(i = 0; i < 3; i++) {
      if(all_eigenvalues.at(i).real() < TOL) {
        //check that the imaginary parts of eigenvector are 0
        if((all_eigenvectors.at(i)[0].imag() < TOL) &&
           (all_eigenvectors.at(i)[1].imag() < TOL) &&
           (all_eigenvectors.at(i)[2].imag() < TOL)) {
          //eigenvec(CART) is a Vector3<double> whereas all_eigenvalues.at(i)
          //is a Vector3<std::complex<double>>
          eigenvec(CART).at(0) = all_eigenvectors.at(i)[0].real();
          eigenvec(CART).at(1) = all_eigenvectors.at(i)[1].real();
          eigenvec(CART).at(2) = all_eigenvectors.at(i)[2].real();

          symmetry = mirror_op;

          glide_check();

          return;
        }
      }
    }
  }//end of MIRROR check

  //******************************************************
  /**
   * Calculates the rotation angle.
   *
   *
   */
  //******************************************************

  double SymOp::get_rotation_angle() const {
    return rotation_angle;
  }

  //*****************************************************

  const Coordinate &SymOp::get_location() const {
    return location;
  }

  //*****************************************************

  const Vector3<double> &SymOp::get_location(COORD_TYPE mode) const {
    return location(mode);
  }
  //*****************************************************
  const Coordinate &SymOp::get_screw_glide_shift() const {
    return screw_glide_shift;
  }
  //*****************************************************
  const Vector3<double> &SymOp::get_screw_glide_shift(COORD_TYPE mode) const {
    return screw_glide_shift(mode);
  }

  //******************************************************
  /**
   * Calculates the rotation angle.
   *
   *
   */
  //******************************************************

  void SymOp::calc_rotation_angle(Eigen::Matrix3d mat, double det, double trace) const {
    int i, j;

    double vec_sum, vec_mag;
    //If rotation is 180 degrees
    //Rotation matrix becomes symmetric; 180 rotation can be
    //decomposed into 2 orthogonal mirror planes
    //handled the same way as in mirror_check
    if(std::abs(trace + 1) < TOL) {
      Vector3<double> test_vec(0, 0, 0), w_vec;

      rotation_angle = 180;

      Vector3< Vector3<std::complex<double> > > all_eigenvectors;
      Vector3< std::complex<double> > all_eigenvalues;
      all_eigenvalues = matrix().eigen(all_eigenvectors);

      //Eigenvalues of 180 rotation are 1, -1, -1
      for(j = 0; j < 3; j++) {
        if(std::abs(all_eigenvalues.at(j).real() - 1) < TOL) {
          if((all_eigenvectors.at(j)[0].imag() < TOL) &&
             (all_eigenvectors.at(j)[1].imag() < TOL) &&
             (all_eigenvectors.at(j)[2].imag() < TOL)) {
            eigenvec(CART).at(0) = all_eigenvectors.at(j)[0].real();
            eigenvec(CART).at(1) = all_eigenvectors.at(j)[1].real();
            eigenvec(CART).at(2) = all_eigenvectors.at(j)[2].real();

            return;
          }
        }
      } //End loop over eigenvectors
    }

    //Extract rotation angle

    // Use this method because arccos is only defined for 0 to pi.
    // Using this axis-angle method produces unique angle of rotation.

    // Following only evaluates if we have non-180 proper rotation
    // Method uses inversion of axis-angle interpretation of a rotation matrix R
    // With axis v=(x,y,z) and angle TH, with ||v||=1
    //  c = cos(TH); s = sin(TH); C = 1-c
    //      [ x*xC+c   xyC-zs   zxC+ys ]
    //  R = [ xyC+zs   y*yC+c   yzC-xs ]
    //      [ zxC-ys   yzC+xs   z*zC+c ]

    vec_sum = 0.0;
    vec_mag = 0.0;

    for(i = 0; i < 3; i++) {
      eigenvec.at(i, CART) = mat((i + 2) % 3, (i + 1) % 3) - mat((i + 1) % 3, (i + 2) % 3);

      //right side of above expression evaluates to
      //2xsin(TH) for i = 0, 2ysin(TH) for i = 1, 2zsin(TH) for i = 2

      //if TH = 0 or 180, second condition won't ever be fulfilled
      if((std::abs(vec_sum) < TOL) && (std::abs(eigenvec.at(i, CART)) > TOL)) {
        //extract sign of eigenvec(i,CART)
        vec_sum = eigenvec.at(i, CART);
      }
    }

    vec_mag = eigenvec(CART).length();

    if(vec_sum < -TOL) {
      //vec_mag = 2sin(TH) and trace - 1 = 2cos(TH)
      //negative vec_mag used to preserve the sign of sin since above,
      //the sign was lost when vec_mag = eigenvec(CART).length()
      rotation_angle = int(round((180.0 / M_PI) * atan2(-vec_mag, (trace - 1)))) + 360;
      eigenvec(CART).normalize();
      eigenvec(CART) *= -1;
    }
    else {
      //vec_mag = 2sin(TH) and trace - 1 = 2cos(TH)
      rotation_angle = int(round((180.0 / M_PI) * atan2(vec_mag, (trace - 1))));
      eigenvec(CART).normalize();
    }

    /*
      eigenvec(CART).normalize();

      if (vec_sum < -TOL)
      {
      rotation_angle = 360 - rotation_angle;
      eigenvec(CART) *= -1;
      }
    */
    return;
  }

  //*****************************************************
  void SymOp::print_short(std::ostream &stream, COORD_TYPE mode) const {
    Matrix3 <double> tsym_mat;

    get_sym_type();

    if(mode == COORD_DEFAULT)
      mode = COORD_MODE::CHECK();


    if(!calc(mode))
      std::cerr << "Trying to print out symmetry matrix, but it is not available in " << COORD_MODE::NAME(mode) << " mode\n";

    stream.precision(3);

    switch(symmetry) {
    case identity_op:
      stream << "Identity Operation \n";
      break;

    case mirror_op:
      stream.setf(std::ios::showpoint);
      stream << "Mirror Operation with plane Normal = " << std::setw(7) << eigenvec(mode) << '\n';
      break;

    case glide_op:
      stream.setf(std::ios::showpoint);
      stream << "Glide Operation with plane Normal = " << std::setw(7) << eigenvec(mode) << '\n'
             << "Glide Vector:" << std::setw(7) << screw_glide_shift(mode) << '\n';
      break;

    case rotation_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << rotation_angle << "-degree Rotation Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << eigenvec(mode)  << '\n';
      break;

    case screw_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << rotation_angle << "-degree Screw Operation along axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << eigenvec(mode) << "\n Screw Vector:" << std::setw(7) << screw_glide_shift(mode) << '\n';
      break;

    case inversion_op:
      stream << "Inversion Operation\n";
      break;

    case rotoinversion_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << rotation_angle << "-degree Rotoinversion Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << eigenvec(mode) << "\n";
      break;

    case invalid_op:
      stream << "Invalid Operation !!! \n";
      break;

    }

  }

  //*****************************************************

  void SymOp::print(std::ostream &stream) const {
    Matrix3 <double> tsym_mat;

    get_sym_type();

    int tprec = stream.precision();
    std::ios::fmtflags tflags = stream.flags();

    stream.precision(3);

    switch(symmetry) {
    case identity_op:
      stream << "Identity Operation \n";
      break;

    case mirror_op:
      stream.setf(std::ios::showpoint);
      stream << "Mirror Operation with plane Normal = " << std::setw(7) << eigenvec() << '\n';
      break;

    case glide_op:
      stream.setf(std::ios::showpoint);
      stream << "Glide Operation with plane Normal = " << std::setw(7) << eigenvec() << '\n'
             << "Glide Vector:" << std::setw(7) << screw_glide_shift() << '\n';
      break;

    case rotation_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << rotation_angle << "-degree Rotation Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << eigenvec()  << '\n';
      break;

    case screw_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << rotation_angle << "-degree Screw Operation along axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << eigenvec() << "\n Screw Vector:" << std::setw(7) << screw_glide_shift() << '\n';
      break;

    case inversion_op:
      stream << "Inversion Operation\n";
      break;

    case rotoinversion_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << rotation_angle << "-degree Rotoinversion Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << eigenvec() << "\n";
      break;

    case invalid_op:
      stream << "Invalid Operation !!! \n";
      break;

    }

    stream.flags(std::ios::left);
    stream << std::setw(53) << "Symmetry Operation Matrix" << "Shift \n"; //SOM has 25 chars, width of those 3 columns are 14 each, so 42 width.  Shift width is 22, so spacing of 9-11 extra characters, so add 5 more to get in the middle

    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream.precision(9);

    for(int i = 0; i < 3; i++) {
      //Print each row of the symmetry matrix separately
      for(int j = 0; j < 3; j++) {
        stream << std::setw(14) << matrix()(i, j);
      }

      stream << std::setw(22) << tau(i) << "\n";
    }

    stream.precision(tprec);
    stream.flags(tflags);
    return;
  }

  //**********************************************************

  jsonParser &SymOp::to_json(jsonParser &json) const {
    json.put_obj();

    // Members not included:
    //
    // From SymOpRepresentation:
    //   MasterSymGroup const master_group();
    //
    // From SymOp:
    //   Lattice const *home;
    //   Array<SymOpRepresentation *> representation;

    json["SymOpRep_type"] = "SymOp";

    ///type of symmetry, given by one of the allowed values of symmetry_type
    json["symmetry"] = symmetry;
    json["op_index"] = op_index;
    json["rep_ID"] = rep_ID;

    // mutable Eigen::Matrix3d symmetry_mat[2];
    json["symmetry_mat"] = matrix();

    // mutable Coordinate tau_vec;
    json["tau"] = tau();

    // mutable Coordinate location;
    find_location();
    location.within();
    if(location.calc())
      json["location"] = location;

    // mutable Coordinate eigenvec;
    if(eigenvec.calc())
      json["eigenvec"] = eigenvec;

    // mutable double rotation_angle;
    json["rotation_angle"] = rotation_angle;

    // mutable Coordinate screw_glide_shift;
    if(screw_glide_shift.calc())
      json["screw_glide_shift"] = screw_glide_shift;

    // double map_error;
    json["map_error"] = map_error;

    return json;
  }

  //**********************************************************

  void SymOp::from_json(const jsonParser &json) {
    try {
      //std::cout<<"Inside of SymOp::from_json"<<std::endl;
      //std::cout<<"Reading in symmetry"<<std::endl;
      CASM::from_json(symmetry, json["symmetry"]);
      //std::cout<<"Reading in op_index"<<std::endl;
      CASM::from_json(op_index, json["op_index"]);
      //std::cout<<"Reading in rep_id"<<std::endl;
      CASM::from_json(rep_ID, json["rep_ID"]);

      // mutable Eigen::Matrix3d symmetry_mat[2];
      //std::cout<<"Reading in symmetry_mat"<<std::endl;
      CASM::from_json(m_mat, json["symmetry_mat"]);

      // mutable Coordinate tau_vec;
      //std::cout<<"Reading in tau_vec"<<std::endl;
      if(json.contains("tau"))
        CASM::from_json(m_tau, json["tau"]);

      //std::cout<<"Reading in location"<<std::endl;
      // mutable Coordinate location;
      if(json.contains("location"))
        CASM::from_json(location, json["location"]);

      //std::cout<<"Reading in eigenvec"<<std::endl;
      // mutable Coordinate eigenvec;
      if(json.contains("eigenvec"))
        CASM::from_json(eigenvec, json["eigenvec"]);

      //std::cout<<"Reading in rotation_angle"<<std::endl;
      // mutable double rotation_angle;
      CASM::from_json(rotation_angle, json["rotation_angle"]);

      //std::cout<<"Reading in screw_glide_shift"<<std::endl;
      // mutable Coordinate screw_glide_shift;
      if(json.contains("screw_glide_shift"))
        CASM::from_json(screw_glide_shift, json["screw_glide_shift"]);

      //std::cout<<"Reading in map_error"<<std::endl;
      // double map_error;
      CASM::from_json(map_error, json["map_error"]);
      //std::cout<<"Done Reading in the SymOp"<<std::endl;
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //**********************************************************

  jsonParser &to_json(const SymOp &sym, jsonParser &json) {
    return sym.to_json(json);
  }

  //**********************************************************
  void from_json(SymOp &sym, const jsonParser &json) {
    try {
      sym.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


}

