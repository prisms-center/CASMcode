#include "casm/symmetry/SymOp.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/casm_io/json_io/container.hh"
namespace CASM {

  const double &SymOp::map_error() const {
    return m_map_error;
  }

  //*******************************************************************************************

  void SymOp::set_map_error(const double &value) {
    m_map_error = value;
    return;
  }


  //*******************************************************************************************

  void SymOp::set_index(const MasterSymGroup &new_group, Index new_index) {
    if((valid_index(new_index) && new_index < new_group.size())
       && (this == &(new_group[new_index]) ||
           almost_equal(matrix(), new_group[new_index].matrix()))) {
      m_master_group = &new_group;
      m_op_index = new_index;
    }
    else {
      m_master_group = &new_group;
      //m_master_group = NULL;
      m_op_index = -1;
    }
  }

  //*******************************************************************************************

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
    //t_op.m_op_index = 0;
    //}
    //else{
    //std::cout << "This symmetry is " << symmetry << " with head " << m_master_group << " and RHS symmetry is " << RHS.symmetry << " with head " << RHS.m_master_group << "\n";
    //}

    return t_op;

  }

  //*******************************************************************************************

  SymOp &SymOp::operator+=(const Eigen::Ref<const SymOp::vector_type> &RHS) {
    m_tau += RHS - matrix() * RHS;
    return (*this);
  }

  //*******************************************************************************************

  SymOp &SymOp::operator-=(const Eigen::Ref<const SymOp::vector_type> &RHS) {
    m_tau -= RHS - matrix() * RHS;
    return (*this);
  }

  //*******************************************************************************************
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
      t_op.m_op_index = 0;
    }

    return t_op;
  }

  //*******************************************************************************************

  SymOp SymOp::no_trans() const {
    return SymOp(matrix(), vector_type::Zero(), map_error(), index(), m_master_group);
  }

  //*******************************************************************************************

  bool SymOp::operator==(const SymOp &RHS) const {
    return
      almost_equal(matrix(), RHS.matrix()) &&
      almost_equal(tau(), RHS.tau());
  };

  //*******************************************************************************************

  SymOp &SymOp::apply_sym(const SymOp &op) {
    (*this) = op * (*this) * (op.inverse());
    return *this;
  }

  //*******************************************************************************************
  /**
   *
   *
   */
  //*******************************************************************************************

  SymOp::SymInfo SymOp::info() const {
    SymInfo result;

    // Simplest case is identity: has no axis and no location
    if(almost_equal(matrix().trace(), 3.)) {
      result.angle = 0;
      result.op_type = identity_op;
      result.axis = vector_type::Zero();
      result.location = vector_type::Zero();
      return result;
    }

    // second simplest case is inversion: has no axis and location is tau()/2
    if(almost_equal(matrix().trace(), -3.)) {
      result.angle = 0;
      result.op_type = inversion_op;
      result.axis = vector_type::Zero();
      result.location = tau() / 2.;
      return result;
    }

    // det is -1 if improper and +1 if proper
    int det = round(matrix().determinant());

    // Find eigen decomposition of proper operation (by multiplying by determinant)
    Eigen::EigenSolver<matrix_type> t_eig(det * matrix());

    // 'axis' is eigenvector whose eigenvalue is +1
    for(Index i = 0; i < 3; i++) {
      if(almost_equal(t_eig.eigenvalues()(i), std::complex<double>(1, 0))) {
        result.axis = t_eig.eigenvectors().col(i).real();
        break;
      }
    }

    // Sign convention for 'axis': first non-zero element is positive
    for(Index i = 0; i < 3; i++) {
      if(!almost_zero(result.axis[i])) {
        result.axis *= sgn(result.axis[i]);
        break;
      }
    }

    vector_type ortho = result.axis.unitOrthogonal();
    vector_type rot = matrix() * ortho;

    result.angle = int(round((180. / M_PI) * atan2(result.axis.dot(ortho.cross(rot)), ortho.dot(rot)))) + 180;
    if(det < 0) {
      if(almost_equal(result.angle, 180.)) {
        result.op_type = mirror_op;
        // shift is component of tau perpendicular to axis
        result.screw_glide_shift = tau() - tau().dot(result.axis) * result.axis;
        // location is 1/2 of component of tau parallel to axis: matrix()*location+tau() = -location+tau() = location
        result.location = tau().dot(result.axis) * result.axis / 2.;
      }
      else {
        result.op_type = rotoinversion_op;
        // shift is component of tau parallel to axis
        result.screw_glide_shift = tau().dot(result.axis) * result.axis;
        // rotoinversion is point symmetry, so we can solve matrix()*p+tau()=p for invariant point p
        result.location = (matrix_type::Identity() - matrix()).inverse() * tau();
      }
    }
    else {
      result.op_type = rotation_op;
      // shift is component of tau parallel to axis
      result.screw_glide_shift = tau().dot(result.axis) * result.axis;
      // Can only solve 2d location problem
      Eigen::MatrixXd tmat(3, 2);
      tmat << ortho, ortho.cross(result.axis);
      // if A = tmat.transpose()*matrix()*tmat and s=tmat.transpose()*tau()
      // then 2d invariant point 'v' is solution to A*v+s=v
      // implies 3d invariant point 'p' is p=tmat*(eye(2)-A).inverse()*s
      result.location = tmat * (Eigen::MatrixXd::Identity(2, 2) - tmat.transpose() * matrix() * tmat).inverse() * tmat.transpose() * tau();
    }
    return result;
  }


  //*******************************************************************************************
  void SymOp::print_short(std::ostream &stream, const Eigen::Ref<const SymOp::matrix_type> &c2f_mat) const {

    stream.precision(3);
    SymInfo t_info = info();

    switch(t_info.op_type) {
    case identity_op:
      stream << "Identity Operation \n";
      break;

    case mirror_op:
      stream.setf(std::ios::showpoint);
      stream << "Mirror Operation with plane Normal = " << std::setw(7) << (c2f_mat * t_info.axis).transpose() << '\n';
      break;

    case glide_op:
      stream.setf(std::ios::showpoint);
      stream << "Glide Operation with plane Normal = " << std::setw(7) << (c2f_mat * t_info.axis).transpose() << '\n'
             << "Glide Vector:" << std::setw(7) << (c2f_mat * t_info.screw_glide_shift).transpose() << '\n';
      break;

    case rotation_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Rotation Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << (c2f_mat * t_info.axis).transpose()  << '\n';
      break;

    case screw_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Screw Operation along axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << (c2f_mat * t_info.axis).transpose() << "\n Screw Vector:" << std::setw(7) << (c2f_mat * t_info.screw_glide_shift).transpose() << '\n';
      break;

    case inversion_op:
      stream << "Inversion Operation\n";
      break;

    case rotoinversion_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Rotoinversion Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << (c2f_mat * t_info.axis).transpose() << "\n";
      break;

    case invalid_op:
      stream << "Invalid Operation !!! \n";
      break;

    }

  }

  //*******************************************************************************************

  void SymOp::print(std::ostream &stream, const Eigen::Ref<const SymOp::matrix_type> &c2f_mat) const {
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
      stream << "Mirror Operation with plane Normal = " << std::setw(7) << (c2f_mat * t_info.axis).transpose() << '\n';
      break;

    case glide_op:
      stream.setf(std::ios::showpoint);
      stream << "Glide Operation with plane Normal = " << std::setw(7) << (c2f_mat * t_info.axis).transpose() << '\n'
             << "Glide Vector:" << std::setw(7) << (c2f_mat * t_info.screw_glide_shift).transpose() << '\n';
      break;

    case rotation_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Rotation Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << (c2f_mat * t_info.axis).transpose()  << '\n';
      break;

    case screw_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Screw Operation along axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << (c2f_mat * t_info.axis).transpose() << "\n Screw Vector:" << std::setw(7) << (c2f_mat * t_info.screw_glide_shift).transpose() << '\n';
      break;

    case inversion_op:
      stream << "Inversion Operation\n";
      break;

    case rotoinversion_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << t_info.angle << "-degree Rotoinversion Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7) << (c2f_mat * t_info.axis).transpose() << "\n";
      break;

    case invalid_op:
      stream << "Invalid Operation !!! \n";
      break;

    }

    stream.flags(std::ios::left);
    stream << std::setw(53) << "Symmetry Operation Matrix" << "Shift \n"; //SOM has 25 chars, width of those 3 columns are 14 each, so 42 width.  Shift width is 22, so spacing of 9-11 extra characters, so add 5 more to get in the middle

    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream.precision(9);
    matrix_type tmat(c2f_mat * matrix()*c2f_mat.inverse());
    vector_type ttau(c2f_mat * tau());
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

  //*******************************************************************************************

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
    json["op_index"] = m_op_index;
    json["rep_ID"] = m_rep_ID;

    // mutable SymOp::matrix_type symmetry_mat[2];
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

  //*******************************************************************************************

  void SymOp::from_json(const jsonParser &json) {
    try {
      //std::cout<<"Inside of SymOp::from_json"<<std::endl;
      //std::cout<<"Reading in op_index"<<std::endl;
      CASM::from_json(m_op_index, json["op_index"]);
      //std::cout<<"Reading in rep_id"<<std::endl;
      CASM::from_json(m_rep_ID, json["rep_ID"]);

      // mutable SymOp::matrix_type symmetry_mat[2];
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

  //*******************************************************************************************

  jsonParser &to_json(const SymOp &sym, jsonParser &json) {
    return sym.to_json(json);
  }

  //*******************************************************************************************
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
