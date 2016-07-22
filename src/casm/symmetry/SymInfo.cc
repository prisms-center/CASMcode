#include "casm/symmetry/SymInfo.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  const std::string traits<symmetry_type>::name = "symmetry_type";

  const std::multimap<symmetry_type, std::vector<std::string> > traits<symmetry_type>::strval = {
    {symmetry_type::identity_op, {"identity"} },
    {symmetry_type::mirror_op, {"mirror"} },
    {symmetry_type::glide_op, {"glide"} },
    {symmetry_type::rotation_op, {"rotation"} },
    {symmetry_type::screw_op, {"screw"} },
    {symmetry_type::inversion_op, {"inversion"} },
    {symmetry_type::rotoinversion_op, {"rotoinversion"} },
    {symmetry_type::invalid_op, {"invalid"} }
  };


  SymInfo::SymInfo(const SymOp &op, const Lattice &lat) :
    axis(lat),
    screw_glide_shift(lat),
    location(lat) {

    auto matrix = op.matrix();
    auto tau = op.tau();

    vector_type _axis;
    vector_type _screw_glide_shift;
    vector_type _location;

    // Simplest case is identity: has no axis and no location
    if(almost_equal(matrix.trace(), 3.)) {
      angle = 0;
      op_type = symmetry_type::identity_op;
      _axis = vector_type::Zero();
      _location = vector_type::Zero();
      _set(_axis, _screw_glide_shift, _location, lat);
      return;
    }

    // second simplest case is inversion: has no axis and location is tau()/2
    if(almost_equal(matrix.trace(), -3.)) {
      angle = 0;
      op_type = symmetry_type::inversion_op;
      _axis = vector_type::Zero();
      _location = tau / 2.;
      _set(_axis, _screw_glide_shift, _location, lat);
      return;
    }

    // det is -1 if improper and +1 if proper
    int det = round(matrix.determinant());

    // Find eigen decomposition of proper operation (by multiplying by determinant)
    Eigen::EigenSolver<matrix_type> t_eig(det * matrix);

    // 'axis' is eigenvector whose eigenvalue is +1
    for(Index i = 0; i < 3; i++) {
      if(almost_equal(t_eig.eigenvalues()(i), std::complex<double>(1, 0))) {
        _axis = t_eig.eigenvectors().col(i).real();
        break;
      }
    }

    // Sign convention for 'axis': first non-zero element is positive
    for(Index i = 0; i < 3; i++) {
      if(!almost_zero(_axis[i])) {
        _axis *= float_sgn(_axis[i]);
        break;
      }
    }

    // get vector orthogonal to axis: ortho,
    // apply matrix: rot
    // and check angle between ortho and det*rot,
    // using determinant to get the correct angle for improper
    // (i.e. want angle before inversion for rotoinversion)
    vector_type ortho = _axis.unitOrthogonal();
    vector_type rot = det * (matrix * ortho);
    angle = fmod((180. / M_PI) * atan2(_axis.dot(ortho.cross(rot)), ortho.dot(rot)) + 360., 360.);

    /*
    std::cout << "det: " << det << "\n";
    std::cout << "y: " << _axis.dot(ortho.cross(rot)) << "\n";
    std::cout << "x: " << ortho.dot(rot) << "\n";
    std::cout << "angle: " << angle << std::endl;
    */

    if(det < 0) {
      if(almost_equal(angle, 180.)) {

        // shift is component of tau perpendicular to axis
        Coordinate coord(tau - tau.dot(_axis) * _axis, lat, CART);
        _screw_glide_shift = coord.cart();

        // location is 1/2 of component of tau parallel to axis:
        //   matrix*location+tau = -location+tau = location
        _location = tau.dot(_axis) * _axis / 2.;

        op_type = coord.is_lattice_shift() ? symmetry_type::mirror_op : symmetry_type::glide_op;
      }
      else {
        // shift is component of tau parallel to axis
        _screw_glide_shift = tau.dot(_axis) * _axis;

        // rotoinversion is point symmetry, so we can solve matrix*p+tau=p for invariant point p
        _location = (matrix_type::Identity() - matrix).inverse() * tau;

        op_type = symmetry_type::rotoinversion_op;
      }
    }
    else {

      // shift is component of tau parallel to axis
      Coordinate coord(tau.dot(_axis) * _axis, lat, CART);
      _screw_glide_shift = coord.cart();

      // Can only solve 2d location problem
      Eigen::MatrixXd tmat(3, 2);
      tmat << ortho, ortho.cross(_axis);

      // if A = tmat.transpose()*matrix()*tmat and s=tmat.transpose()*tau()
      // then 2d invariant point 'v' is solution to A*v+s=v
      // implies 3d invariant point 'p' is p=tmat*(eye(2)-A).inverse()*s
      _location = tmat * (Eigen::MatrixXd::Identity(2, 2) - tmat.transpose() * matrix * tmat).inverse() * tmat.transpose() * tau;

      op_type = coord.is_lattice_shift() ? symmetry_type::rotation_op : symmetry_type::screw_op;

    }
    _set(_axis, _screw_glide_shift, _location, lat);
    return;
  }

  void SymInfo::_set(const vector_type &_axis,
                     const vector_type &_screw_glide_shift,
                     const vector_type &_location,
                     const Lattice &lat) {
    axis = Coordinate(_axis, lat, CART);
    screw_glide_shift = Coordinate(_screw_glide_shift, lat, CART);
    location = Coordinate(_location, lat, CART);
  }


  /// \brief Print SymInfo to string
  ///
  /// Of the form:
  /// \code
  /// Mirror Operation with plane Normal = 0.25 0.25 0.0
  std::string to_string(const SymInfo &info, COORD_TYPE mode) {

    std::stringstream stream;
    char term(0);
    int prec(3);
    int pad(0);

    auto print_coord = [&](const Coordinate & coord) {
      coord.print(stream, mode, term, prec, pad);
    };

    auto print_axis = [&](const Coordinate & coord) {
      coord.print_axis(stream, mode, term, prec, pad);
    };

    switch(info.op_type) {
    case symmetry_type::identity_op:
      stream << "Identity Operation";
      break;

    case symmetry_type::mirror_op:
      stream.setf(std::ios::showpoint);
      stream << "Mirror Operation with plane normal" << std::setw(7);
      print_axis(info.axis);
      break;

    case symmetry_type::glide_op:
      stream.setf(std::ios::showpoint);
      stream << "Glide Operation with plane normal" << std::setw(7);
      print_axis(info.axis);
      stream << '\n';
      stream << "Glide Vector:" << std::setw(7);
      print_coord(info.screw_glide_shift);
      break;

    case symmetry_type::rotation_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << info.angle << "-degree Rotation Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7);
      print_axis(info.axis);
      break;

    case symmetry_type::screw_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << info.angle << "-degree Screw Operation along axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7);
      print_axis(info.axis);
      stream << "\n Screw Vector:" << std::setw(7);
      print_coord(info.screw_glide_shift);
      break;

    case symmetry_type::inversion_op:
      stream << "Inversion Operation";
      break;

    case symmetry_type::rotoinversion_op:
      stream.unsetf(std::ios::showpoint);
      stream << std::setprecision(3) << info.angle << "-degree Rotoinversion Operation about axis";
      stream.setf(std::ios::showpoint);
      stream << std::setw(7);
      print_axis(info.axis);
      break;

    case symmetry_type::invalid_op:
      stream << "Invalid Operation !!!";
      break;

    }

    return stream.str();
  }

  /// \brief Print SymInfo to string
  std::string description(const SymOp &op, const Lattice &lat, COORD_TYPE mode) {
    return to_string(SymInfo(op, lat), mode);
  }

  /// \brief Add to existing JSON object
  void add_sym_info(const SymInfo &info, jsonParser &j) {
    ///type of symmetry, given by one of the allowed values of symmetry_type
    j["type"] = info.op_type;

    if(info.op_type == symmetry_type::rotation_op ||
       info.op_type == symmetry_type::screw_op ||
       info.op_type == symmetry_type::rotoinversion_op) {
      to_json_array(info.axis.const_cart().normalized(), j["rotation_axis"]["CART"]);
      to_json_array(info.axis.const_frac().normalized(), j["rotation_axis"]["FRAC"]);
      j["rotation_angle"] = info.angle;
    }
    else if(info.op_type == symmetry_type::mirror_op ||
            info.op_type == symmetry_type::glide_op) {
      to_json_array(info.axis.const_cart().normalized(), j["mirror_normal"]["CART"]);
      to_json_array(info.axis.const_frac().normalized(), j["mirror_normal"]["FRAC"]);
    }

    if(info.op_type == symmetry_type::screw_op ||
       info.op_type == symmetry_type::glide_op) {
      to_json_array(info.screw_glide_shift.const_cart(), j["shift"]["CART"]);
      to_json_array(info.screw_glide_shift.const_frac(), j["shift"]["FRAC"]);
    }

    if(info.op_type != symmetry_type::identity_op) {
      to_json_array(info.location.const_cart(), j["invariant_point"]["CART"]);
      to_json_array(info.location.const_frac(), j["invariant_point"]["FRAC"]);
    }
  }

}
