#include "casm/symmetry/SymOp.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/casm_io/json_io/container.hh"
namespace CASM {

  //SymOpRepresentation::~SymOpRepresentation() {};

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

  void SymOp::set_index(const MasterSymGroup &new_group, Index new_index) {
    if((valid_index(new_index) && new_index < new_group.size())
       && (this == &(new_group[new_index]) ||
           almost_equal(matrix(), new_group[new_index].matrix()))) {
      m_master_group = &new_group;
      op_index = new_index;
    }
    else {
      m_master_group = NULL;
      op_index = -1;
    }
  }

  //**********************************************************

  SymOp SymOp::operator*(const SymOp &RHS) const {
    SymOp t_op(matrix() * RHS.matrix(),
               tau() + matrix() * RHS.tau());

    if(m_master_group && (m_master_group == RHS.m_master_group)) {
      t_op.set_index(master_group(), master_group().ind_prod(index(), RHS.index()));
    }
    else if(RHS.m_master_group && !m_master_group && is_identity()) {
      t_op.set_index(*RHS.m_master_group, RHS.index());
    }
    else if(m_master_group && !RHS.m_master_group && RHS.is_identity()) {
      t_op.set_index(master_group(), index());
    }
    //The following blocks caused problems at some point (mainly for non-primitive structures)
    //else if(is_identity() && RHS.is_identity()) {
    //t_op.op_index = 0;
    //}
    //else{
    //std::cout << "This symmetry is " << symmetry << " with head " << m_master_group << " and RHS symmetry is " << RHS.symmetry << " with head " << RHS.m_master_group << "\n";
    //}

    return t_op;

  }

  //**********************************************************

  SymOp &SymOp::operator+=(const Eigen::Ref<const SymOp::vector_type> &RHS) {
    m_tau += RHS - matrix() * RHS;
    return (*this);
  }

  //**********************************************************

  SymOp &SymOp::operator-=(const Eigen::Ref<const SymOp::vector_type> &RHS) {
    m_tau -= RHS - matrix() * RHS;
    return (*this);
  }

  //**********************************************************
  // SymOp matrix is unitary, so inverse is equivalent to transpose.
  // To do inverse of translation, you must perform
  // inverse matrix operaton on translation and subtract
  SymOp SymOp::inverse() const {
    SymOp t_op(matrix().transpose(),
               -(matrix().transpose() * tau()));
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

    return SymOp(matrix(), map_error());
  }

  //**********************************************************
  /*
  bool SymOp::compare(const SymOp &RHS, double eq_tol) const {
    return
      almost_equal(matrix(),RHS.matrix(), eq_tol) &&
      tau_vec.min_dist(RHS.tau_vec) < eq_tol;

      }
  */
  //**********************************************************

  bool SymOp::operator==(const SymOp &RHS) const {
    return
      almost_equal(matrix(), RHS.matrix()) &&
      almost_equal(tau(), RHS.tau());
  };

  //******************************************************
  /**
   * Find the invariant point after applying symmetry

   *	- General equation

   * P is the invariant point
   * P(CART) = Sym_op*P(CART) + tau()
   * P(CART) = (I-Sym_op).inverse()*tau()

   *	- The General equation is only valid for point sysmetry

   *	- Plane symmetry ( Mirror or Glide)
   *  SP = -P (always sture at the origin)
   *  P(CART) = -P(CART) + tau_per (perpendicular to eigenvec)
   *  P = tau_per/2

   *	- axial sysmmetry ( Rotation or Screw)
   *   change the coordinate system having eigen vector as a z coordinate,
   *	 so, it makes 3D rotation to 2D roation
   *   find M coordinate transfer matrix
   *   M = ( tau_pp ; eigen x tau_pp ; eigen)
   *	 define Inew = (1 0 0; 0 1 0; 0 0 0)
   *   Pnew = (Inew-Snew)^-1*tau_pp(new)
   *   P  = (MInewM.inverse() - S).inverse*tau_pp
   */
  /*
  void SymOp::find_location() const {


    Eigen::Matrix3d tMat, inv_tMat;
    Eigen::Vector3d tau_pp, tau_ll;

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
      location(CART) = inv_tMat * tau();

      return;
    }
    if((type() == mirror_op) || (type() == glide_op)) {
      //component of tau parallel to eigenvector:
      tau_ll = (eigenvec(CART).dot(tau()) / eigenvec(CART).dot(eigenvec(CART))) * eigenvec(CART);

      location(CART) = tau_ll / 2.0;
      return;
    }
    if(((type() == rotation_op) || (type() == screw_op))) {


      // std::cout << "eigenvec norm is " << eigenvec(CART).norm() << '\n';

      //component of tau parallel to rotation axis
      tau_ll = (eigenvec(CART).dot(tau()) / eigenvec(CART).dot(eigenvec(CART))) * eigenvec(CART);

      //component in the plane of rotation
      tau_pp = tau() - tau_ll;

      if(tau_pp.is_zero()) {
        //rotation axis passes through origin
        //std::cout << location() << "\n";
        return;
      }


      Eigen::Vector3d X, Y, Z, tY;
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
  */
  //******************************************************

  SymOp &SymOp::apply_sym(const SymOp &op) {
    (*this) = op * (*this) * (op.inverse());
    return *this;
  }

  //******************************************************
  /**
   * Calculates the rotation angle.
   *
   *
   */
  //******************************************************

  SymOp::SymInfo SymOp::info() const {
    SymInfo result;
    if(almost_equal(matrix().trace(), -3.0)) {
      result.op_type = identity_op;
      return result;
    }

    if(almost_equal(matrix().trace(), -3.0)) {
      result.op_type = inversion_op;
      return result;
    }

    int det = round(matrix().determinant());

    int i, j;

    //If rotation is 180 degrees
    //Rotation matrix becomes symmetric; 180 rotation can be
    //decomposed into 2 orthogonal mirror planes
    //handled the same way as in mirror_check
    Eigen::EigenSolver<matrix_type> t_eig(det * matrix());

    //Eigenvalues of 180 rotation are 1, -1, -1
    for(i = 0; i < 3; i++) {
      if(almost_equal(t_eig.eigenvalues()(i), std::complex<double>(1, 0))) {
        result.axis = t_eig.eigenvectors().col(i).real();
        break;
      }
    }
    for(i = 0; i < 3; i++) {
      if(!almost_zero(result.axis[i])) {
        result.axis *= sgn(result.axis[i]);
        break;
      }
    }

    vector_type ortho = result.axis.unitOrthogonal();
    vector_type rot = matrix() * ortho;
    if(almost_equal(ortho.dot(rot), 1.0))
      result.angle = 0;

    result.angle = int(round((180.0 / M_PI) * atan2(result.axis.dot(ortho.cross(rot)), ortho.dot(rot)))) + 180;
    if(det < 0) {
      if(almost_equal(result.angle, 180.0)) {
        result.op_type = mirror_op;
        result.screw_glide_shift = tau() - tau().dot(result.axis) * result.axis;
      }
      else {
        result.screw_glide_shift = tau().dot(result.axis) * result.axis;
        result.op_type = rotoinversion_op;
      }
    }
    else {
      result.screw_glide_shift = tau().dot(result.axis) * result.axis;
      result.op_type = rotation_op;
    }
    return result;
  }

  //*****************************************************
  void SymOp::print_short(std::ostream &stream, const Eigen::Ref<const Eigen::Matrix3d> &c2f_mat) const {

    stream.precision(3);
    SymInfo t_info = info();

    switch(t_info.op_type) {
    case identity_op:
      stream << "Identity Operation \n";
      break;

    case mirror_op:
      stream.setf(std::ios::showpoint);
      stream << "Mirror Operation with plane Normal = " << std::setw(7) << c2f_mat *t_info.axis << '\n';
      break;

    case glide_op:
      stream.setf(std::ios::showpoint);
      stream << "Glide Operation with plane Normal = " << std::setw(7) << c2f_mat *t_info.axis << '\n'
             << "Glide Vector:" << std::setw(7) << c2f_mat *t_info.screw_glide_shift << '\n';
      break;

    case rotation_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Rotation Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << c2f_mat *t_info.axis  << '\n';
      break;

    case screw_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Screw Operation along axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << c2f_mat *t_info.axis << "\n Screw Vector:" << std::setw(7) << c2f_mat *t_info.screw_glide_shift << '\n';
      break;

    case inversion_op:
      stream << "Inversion Operation\n";
      break;

    case rotoinversion_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Rotoinversion Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << c2f_mat *t_info.axis << "\n";
      break;

    case invalid_op:
      stream << "Invalid Operation !!! \n";
      break;

    }

  }

  //*****************************************************

  void SymOp::print(std::ostream &stream, const Eigen::Ref<const Eigen::Matrix3d> &c2f_mat) const {
    SymInfo t_info = info();

    int tprec = stream.precision();
    std::ios::fmtflags tflags = stream.flags();

    stream.precision(3);

    switch(info().op_type) {
    case identity_op:
      stream << "Identity Operation \n";
      break;

    case mirror_op:
      stream.setf(std::ios::showpoint);
      stream << "Mirror Operation with plane Normal = " << std::setw(7) << c2f_mat *t_info.axis << '\n';
      break;

    case glide_op:
      stream.setf(std::ios::showpoint);
      stream << "Glide Operation with plane Normal = " << std::setw(7) << c2f_mat *t_info.axis << '\n'
             << "Glide Vector:" << std::setw(7) << c2f_mat *t_info.screw_glide_shift << '\n';
      break;

    case rotation_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Rotation Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << c2f_mat *t_info.axis  << '\n';
      break;

    case screw_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Screw Operation along axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << c2f_mat *t_info.axis << "\n Screw Vector:" << std::setw(7) << c2f_mat *t_info.screw_glide_shift << '\n';
      break;

    case inversion_op:
      stream << "Inversion Operation\n";
      break;

    case rotoinversion_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Rotoinversion Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << c2f_mat *t_info.axis << "\n";
      break;

    case invalid_op:
      stream << "Invalid Operation !!! \n";
      break;

    }

    stream.flags(std::ios::left);
    stream << std::setw(53) << "Symmetry Operation Matrix" << "Shift \n"; //SOM has 25 chars, width of those 3 columns are 14 each, so 42 width.  Shift width is 22, so spacing of 9-11 extra characters, so add 5 more to get in the middle

    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream.precision(9);
    Eigen::Matrix3d tmat(c2f_mat * matrix()*c2f_mat.inverse());
    Eigen::Vector3d ttau(c2f_mat * tau());
    for(int i = 0; i < 3; i++) {
      //Print each row of the symmetry matrix separately
      for(int j = 0; j < 3; j++) {
        stream << std::setw(14) << tmat(i, j);
      }

      stream << std::setw(22) << ttau(i) << "\n";
    }

    stream.precision(tprec);
    stream.flags(tflags);
    return;
  }

  //**********************************************************

  jsonParser &SymOp::to_json(jsonParser &json) const {
    json.put_obj();

    auto t_info = info();

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
    json["symmetry"] = t_info.op_type;
    json["op_index"] = op_index;
    json["rep_ID"] = rep_ID;

    // mutable Eigen::Matrix3d symmetry_mat[2];
    json["symmetry_mat"] = matrix();

    // mutable Coordinate tau_vec;
    json["tau"] = tau();

    // mutable Coordinate location;
    //json["location"] = location;

    // mutable Coordinate eigenvec;
    json["eigenvec"] = t_info.axis;

    // mutable double rotation_angle;
    json["rotation_angle"] = t_info.angle;

    // mutable Coordinate screw_glide_shift;
    json["screw_glide_shift"] = t_info.screw_glide_shift;

    // double map_error;
    json["map_error"] = map_error();

    return json;
  }

  //**********************************************************

  void SymOp::from_json(const jsonParser &json) {
    try {
      //std::cout<<"Inside of SymOp::from_json"<<std::endl;
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

      //std::cout<<"Reading in map_error"<<std::endl;
      // double map_error;
      CASM::from_json(m_map_error, json["map_error"]);
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

