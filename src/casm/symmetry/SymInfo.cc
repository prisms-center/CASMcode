#include "casm/symmetry/SymInfo.hh"

#include "boost/lexical_cast.hpp"
#include "casm/global/enum.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/symmetry/SymGroup.hh"

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

  ENUM_IO_DEF(symmetry_type)
  ENUM_JSON_IO_DEF(symmetry_type)

  jsonParser &to_json(const SymInfoOptions &opt, jsonParser &json) {
    json.put_obj();
    json[traits<COORD_TYPE>::name] = opt.coord_type;
    json["tol"] = opt.tol;
    json["prec"] = opt.prec;
    json["print_matrix_tau"] = opt.print_matrix_tau;
    return json;
  }

  /// \brief Read from JSON
  void from_json(SymInfoOptions &opt, const jsonParser &json) {
    json.get_if(opt.coord_type, traits<COORD_TYPE>::name);
    json.get_if(opt.tol, "tol");
    json.get_if(opt.prec, "prec");
    json.get_if(opt.print_matrix_tau, "print_matrix_tau");
  }

  SymInfoOptions jsonConstructor<SymInfoOptions>::from_json(const jsonParser &json) {
    SymInfoOptions res;
    CASM::from_json(res, json);
    return res;
  }

  SymInfo::SymInfo(const SymOp &op, const Lattice &lat, SymInfoOptions opt) :
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


  /// \brief Print SymInfo
  ///
  /// Of the form:
  /// \code
  /// Mirror Operation with plane Normal = 0.25 0.25 0.0
  /// \endcode
  void print_sym_info(Log &log, const SymInfo &info, SymInfoOptions opt) {

    std::stringstream stream;
    char term(0);
    Eigen::IOFormat format(opt.prec, opt.prec + 1);

    auto print_coord = [&](const Coordinate & coord) {
      coord.print(log, opt.coord_type, term, format);
    };

    auto print_axis = [&](const Coordinate & coord) {
      coord.print_axis(log, opt.coord_type, term, format);
    };

    switch(info.op_type) {
    case symmetry_type::identity_op:
      log << log.indent_str() << "Identity Operation";
      break;

    case symmetry_type::mirror_op:
      log.ostream().setf(std::ios::showpoint);
      log << log.indent_str() << "Mirror Operation with plane normal ";
      print_axis(info.axis);
      break;

    case symmetry_type::glide_op:
      stream.setf(std::ios::showpoint);
      log << log.indent_str() << "Glide Operation with plane normal ";
      print_axis(info.axis);
      log << '\n';
      log << log.indent_str() << "Glide Vector: ";
      print_coord(info.screw_glide_shift);
      break;

    case symmetry_type::rotation_op:
      log.ostream().unsetf(std::ios::showpoint);
      log << log.indent_str() << std::setprecision(opt.prec) << info.angle << "-degree Rotation Operation about axis ";
      log.ostream().setf(std::ios::showpoint);
      print_axis(info.axis);
      break;

    case symmetry_type::screw_op:
      log.ostream().unsetf(std::ios::showpoint);
      log << log.indent_str() << std::setprecision(opt.prec) << info.angle << "-degree Screw Operation along axis ";
      log.ostream().setf(std::ios::showpoint);
      log << log.indent_str() << std::setw(opt.prec + 4);
      print_axis(info.axis);
      log << log.indent_str() << "\n Screw Vector: ";
      print_coord(info.screw_glide_shift);
      break;

    case symmetry_type::inversion_op:
      log << log.indent_str() << "Inversion Operation";
      break;

    case symmetry_type::rotoinversion_op:
      log.ostream().unsetf(std::ios::showpoint);
      log << log.indent_str() << std::setprecision(opt.prec) << info.angle << "-degree Rotoinversion Operation about axis ";
      log.ostream().setf(std::ios::showpoint);
      print_axis(info.axis);
      break;

    case symmetry_type::invalid_op:
      log << log.indent_str() << "Invalid Operation !!!";
      break;

    }
  }

  /// \brief Print SymInfo
  ///
  /// Of the form:
  /// \code
  /// Mirror Operation with plane Normal = 0.25 0.25 0.0
  /// \endcode
  std::string to_string(const SymInfo &info, SymInfoOptions opt) {
    std::stringstream ss;
    Log log(ss);
    print_sym_info(log, info, opt);
    return ss.str();
  }

  struct PolyWriter {

    std::stringstream ss;
    SymInfoOptions opt;
    bool add;

    PolyWriter(SymInfoOptions _opt = SymInfoOptions()) :
      opt(_opt),
      add(false) {}

    void add_constant(double val) {
      if(!almost_zero(val, opt.tol)) {
        if(add) {
          if(val < 0) {
            ss << "-";
            val *= -1.;
          }
          else {
            ss << "+";
          }
        }
        ss << std::setprecision(opt.prec) << val;
        add = true;
      }
    };

    void add_term(double coeff, std::string var) {
      if(!almost_zero(coeff, opt.tol)) {
        if(add && coeff > 0) {
          ss << "+";
        }
        else if(coeff < 0) {
          ss << "-";
          coeff *= -1.;
        }

        if(almost_equal(coeff, 1., opt.tol)) {
          ss << var;
        }
        else {
          ss << std::setprecision(opt.prec) << coeff << "*" << var;
        }
        add = true;
      }
    };

    // normalize vectors and point in direction to make first coeff positive
    static void standardize(Eigen::Vector3d &v, double tol = TOL) {
      v.normalize();
      for(int i = 0; i < 3; ++i) {
        if(almost_zero(v(i), tol)) {
          continue;
        }
        if(v(i) < 0.) {
          v = -1.*v;
        }
        break;
      }
    };

  };

  /// \brief Use axis and invariant point to return line in '0, y, 0'-type notation
  std::string sym_line(const Coordinate &axis, const Coordinate &point, SymInfoOptions opt) {

    std::stringstream ss;
    double tol = opt.tol;

    Eigen::Vector3d a = axis.as_vec(opt.coord_type);
    Eigen::Vector3d p = point.as_vec(opt.coord_type);

    // find first non-zero axis coordinate, to choose 'x', 'y', or 'z' as param
    std::string _param = "xyz";
    int min_i = -1;
    for(int i = 0; i < 3; ++i) {
      if(!almost_zero(a[i], tol)) {
        min_i = i;
        break;
      }
    }
    std::string param = _param.substr(min_i, 1);

    PolyWriter::standardize(a, tol);

    if(opt.coord_type == FRAC) {
      a = scale_to_int(a).cast<double>();
    }

    // avoid lines of form '0, 1.8671, 1.2922+z' in favor of '0, 1.8671, z'
    for(int i = 0; i < 3; ++i) {
      if(almost_zero(a[(i + 1) % 3]) && almost_zero(a[(i + 2) % 3])) {
        p[i] = 0.;
      }
    }


    for(int i = 0; i < 3; ++i) {
      PolyWriter writer(opt);
      writer.add_constant(p[i]);
      writer.add_term(a[i], param);

      if(!writer.ss.str().size()) {
        ss << 0;
      }
      else {
        ss << writer.ss.str();
      }

      if(i != 2) {
        ss << ", ";
      }
    }

    return ss.str();
  }

  /// \brief Use two perpendicular vectors in plane and invariant point to return plane in 'x, y, 0'-type notation
  std::string sym_plane(const Coordinate &v1, const Coordinate &v2, const Coordinate &point, SymInfoOptions opt) {

    std::stringstream ss;
    double tol = opt.tol;

    Eigen::Vector3d a = v1.as_vec(opt.coord_type);
    Eigen::Vector3d b = v2.as_vec(opt.coord_type);
    Eigen::Vector3d p = point.as_vec(opt.coord_type);

    // in most cases, having 'a*x+b*y' in the same term can be avoided,
    // but in some cases it is impossible:
    //   ex: in primitive fcc, the mirror plane: '2*x, 2*y, -x-y'
    //
    // here we mostly avoid such representations
    // and situations 'x' is used for 'a' and neither 'y' or 'z' can be used for 'b'
    //
    for(int i = 0; i < 3; ++i) {
      if(!almost_zero(a[i], tol) && !almost_zero(b[i], tol)) {
        b = b - b[i] / a[i] * a;
        break;
      }
    }

    for(int i = 0; i < 3; ++i) {
      if(!almost_zero(a[i], tol) && !almost_zero(b[i], tol)) {
        a = a - a[i] / b[i] * b;
        break;
      }
    }

    // find first non-zero axis coordinate, to choose 'x', 'y', or 'z' as param
    std::string _param = "xyz";
    int min_i = -1;
    for(int i = 0; i < 3; ++i) {
      if(!almost_zero(a[i], tol)) {
        min_i = i;
        break;
      }
    }
    int min_j = -1;
    for(int j = 1; j < 3; ++j) {
      if(!almost_zero(b[j], tol) && j != min_i) {
        min_j = j;
        break;
      }
    }

    // associate lesser of 'x','y','z' with 'a'
    if(min_j < min_i) {
      using std::swap;
      swap(a, b);
      swap(min_j, min_i);
    }

    std::string param1 = _param.substr(min_i, 1);
    std::string param2 = _param.substr(min_j, 1);

    if(opt.coord_type == FRAC) {
      PolyWriter::standardize(a, tol);
      PolyWriter::standardize(b, tol);

      // fractional coords can be scaled to int
      a = scale_to_int(a).cast<double>();
      b = scale_to_int(b).cast<double>();
    }

    // p + a*x + b*y, where p is a constant point,
    //   a, b are vectors of coeffs, and param1 and param2 are variable names
    for(int i = 0; i < 3; ++i) {
      PolyWriter writer(opt);
      writer.add_constant(p[i]);
      writer.add_term(a[i], param1);
      writer.add_term(b[i], param2);

      if(!writer.ss.str().size()) {
        ss << 0;
      }
      else {
        ss << writer.ss.str();
      }

      if(i != 2) {
        ss << ", ";
      }
    }

    return ss.str();
  }

  /// \brief Use axis and invariant point to return plane in 'x, y, 0'-type notation
  std::string sym_plane(const Coordinate &axis, const Coordinate &point, SymInfoOptions opt) {

    double tol = opt.tol;

    Eigen::Vector3d n = axis.const_cart();
    Eigen::Vector3d _v1;
    Eigen::Vector3d _v2;
    n.normalize();

    // get first vector orthogonal to axis
    if(!almost_equal(std::abs(n(0)), 1., tol)) {
      Eigen::Vector3d x;
      x << 1., 0., 0.;
      _v1 = x - x.dot(n) * n;
    }
    else {
      Eigen::Vector3d y;
      y << 0., 1., 0.;
      _v1 = y - y.dot(n) * n;
    }
    PolyWriter::standardize(_v1, tol);

    // get second vector in symmetry plane
    _v2 = _v1.cross(n);
    PolyWriter::standardize(_v2, tol);

    // use plane defining vectors and invariant point to write the plane equation
    Coordinate v1(_v1, axis.home(), CART);
    Coordinate v2(_v2, axis.home(), CART);
    return sym_plane(v1, v2, point, opt);
  }



  /// \brief Print symmetry symbol to string
  ///
  /// Of the form:
  /// \code
  /// m 0.25 0.25 0.0
  /// \endcode
  std::string to_brief_unicode(const SymInfo &info, SymInfoOptions opt) {

    std::stringstream stream;
    char term(0);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    Eigen::IOFormat format(opt.prec, opt.prec + 1);

    // print rotation or screw integer, and include trailing space
    auto print_r = [&](bool overline) {
      int r = std::lround(360. / std::min(info.angle, 360. - info.angle));
      if(r == 2) {
        stream << boost::lexical_cast<std::string>(r);
        if(overline) stream << "\u0305";
        stream << " ";
      }
      else if(info.angle < 180.) {
        stream << boost::lexical_cast<std::string>(r);
        if(overline) stream << "\u0305";
        // add superscript '+'
        stream << "\u207A ";
      }
      else {
        stream << boost::lexical_cast<std::string>(r);
        if(overline) stream << "\u0305";
        // add superscript '-'
        stream << "\u207B ";
      }

    };

    auto print_coord = [&](const Coordinate & coord) {
      coord.print(stream, opt.coord_type, term, format);
    };

    switch(info.op_type) {
    case symmetry_type::identity_op:
      stream << "1";
      break;

    case symmetry_type::mirror_op:
      stream.setf(std::ios::showpoint);
      stream << "m ";
      stream << sym_plane(info.axis, info.location, opt);
      break;

    case symmetry_type::glide_op:
      stream.setf(std::ios::showpoint);
      stream << "g (";
      print_coord(info.screw_glide_shift); // glide shift
      stream << ") ";
      stream << sym_plane(info.axis, info.location, opt);
      break;

    case symmetry_type::rotation_op:
      print_r(false);
      stream << sym_line(info.axis, info.location, opt);
      break;

    case symmetry_type::screw_op:
      print_r(false);
      stream << "(";
      print_coord(info.screw_glide_shift); // screw shift
      stream << ") " << sym_line(info.axis, info.location, opt);
      break;

    case symmetry_type::inversion_op:
      stream << "1\u0305 ";
      print_coord(info.location); // inversion point
      break;

    case symmetry_type::rotoinversion_op:
      print_r(true);
      stream << sym_line(info.axis, info.location, opt);
      stream << "; ";
      print_coord(info.location); // invariant point
      break;

    case symmetry_type::invalid_op:
      stream << "<error invalid>";
      break;

    }

    return stream.str();
  }

  /// \brief Print SymInfo to string
  std::string description(const SymOp &op, const Lattice &lat, SymInfoOptions opt) {
    return to_string(SymInfo(op, lat, opt), opt);
  }

  /// \brief Print SymGroup info
  void description(Log &log, const SymGroup &g, const Lattice &lat, SymInfoOptions opt) {
    Index i = 1;
    for(const auto &op : g) {
      log << log.indent_str() << i  << ": (" << op.index() << ") " << description(op, lat, opt) << std::endl;
      if(opt.print_matrix_tau) {
        print_matrix_tau_col(log, op, opt.prec);
      }
      ++i;
    }
  }

  /// \brief Print SymInfo to brief string
  std::string brief_description(const SymOp &op, const Lattice &lat, SymInfoOptions opt) {
    return to_brief_unicode(SymInfo(op, lat, opt), opt);
  }

  /// \brief Print SymGroup with brief strings
  void brief_description(Log &log, const SymGroup &g, const Lattice &lat, SymInfoOptions opt) {
    Index i = 1;
    for(const auto &op : g) {
      log << log.indent_str() << i  << ": (" << op.index() << ") " << brief_description(op, lat, opt) << std::endl;
      if(opt.print_matrix_tau) {
        print_matrix_tau_col(log, op, opt.prec);
      }
      ++i;
    }
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
