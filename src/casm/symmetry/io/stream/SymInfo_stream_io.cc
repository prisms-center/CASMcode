#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"

#include "boost/lexical_cast.hpp"
#include "casm/casm_io/Log.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymInfo.hh"

namespace CASM {

const std::string traits<symmetry_type>::name = "symmetry_type";

const std::multimap<symmetry_type, std::vector<std::string> >
    traits<symmetry_type>::strval = {
        {symmetry_type::identity_op, {"identity"}},
        {symmetry_type::mirror_op, {"mirror"}},
        {symmetry_type::glide_op, {"glide"}},
        {symmetry_type::rotation_op, {"rotation"}},
        {symmetry_type::screw_op, {"screw"}},
        {symmetry_type::inversion_op, {"inversion"}},
        {symmetry_type::rotoinversion_op, {"rotoinversion"}},
        {symmetry_type::invalid_op, {"invalid"}}};

ENUM_IO_DEF(symmetry_type)

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

  auto print_coord = [&](const xtal::Coordinate &coord) {
    coord.print(log, opt.coord_type, term, format);
  };

  auto print_axis = [&](const xtal::Coordinate &coord) {
    coord.print_axis(log, opt.coord_type, term, format);
  };

  switch (info.op_type) {
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
      log << log.indent_str() << std::setprecision(opt.prec) << info.angle
          << "-degree Rotation Operation about axis ";
      log.ostream().setf(std::ios::showpoint);
      print_axis(info.axis);
      break;

    case symmetry_type::screw_op:
      log.ostream().unsetf(std::ios::showpoint);
      log << log.indent_str() << std::setprecision(opt.prec) << info.angle
          << "-degree Screw Operation along axis ";
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
      log << log.indent_str() << std::setprecision(opt.prec) << info.angle
          << "-degree Rotoinversion Operation about axis ";
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

  PolyWriter(SymInfoOptions _opt = SymInfoOptions()) : opt(_opt), add(false) {}

  void add_constant(double val) {
    if (!almost_zero(val, opt.tol)) {
      if (add) {
        if (val < 0) {
          ss << "-";
          val *= -1.;
        } else {
          ss << "+";
        }
      }
      ss << std::setprecision(opt.prec) << val;
      add = true;
    }
  };

  void add_term(double coeff, std::string var) {
    if (!almost_zero(coeff, opt.tol)) {
      if (add && coeff > 0) {
        ss << "+";
      } else if (coeff < 0) {
        ss << "-";
        coeff *= -1.;
      }

      if (almost_equal(coeff, 1., opt.tol)) {
        ss << var;
      } else {
        ss << std::setprecision(opt.prec) << coeff << "*" << var;
      }
      add = true;
    }
  };

  // normalize vectors and point in direction to make first coeff positive
  static void standardize(Eigen::Vector3d &v, double tol = TOL) {
    v.normalize();
    for (int i = 0; i < 3; ++i) {
      if (almost_zero(v(i), tol)) {
        continue;
      }
      if (v(i) < 0.) {
        v = -1. * v;
      }
      break;
    }
  };
};

/// \brief Use axis and invariant point to return line in '0, y, 0'-type
/// notation
std::string sym_line(const xtal::Coordinate &axis,
                     const xtal::Coordinate &point, SymInfoOptions opt) {
  std::stringstream ss;
  double tol = opt.tol;

  Eigen::Vector3d a = axis.as_vec(opt.coord_type);
  Eigen::Vector3d p = point.as_vec(opt.coord_type);

  // find first non-zero axis coordinate, to choose 'x', 'y', or 'z' as param
  std::string _param = "xyz";
  int min_i = -1;
  for (int i = 0; i < 3; ++i) {
    if (!almost_zero(a[i], tol)) {
      min_i = i;
      break;
    }
  }
  std::string param = _param.substr(min_i, 1);

  PolyWriter::standardize(a, tol);

  if (opt.coord_type == FRAC) {
    a = scale_to_int(a).cast<double>();
  }

  // avoid lines of form '0, 1.8671, 1.2922+z' in favor of '0, 1.8671, z'
  for (int i = 0; i < 3; ++i) {
    if (almost_zero(a[(i + 1) % 3]) && almost_zero(a[(i + 2) % 3])) {
      p[i] = 0.;
    }
  }

  for (int i = 0; i < 3; ++i) {
    PolyWriter writer(opt);
    writer.add_constant(p[i]);
    writer.add_term(a[i], param);

    if (!writer.ss.str().size()) {
      ss << 0;
    } else {
      ss << writer.ss.str();
    }

    if (i != 2) {
      ss << ", ";
    }
  }

  return ss.str();
}

/// \brief Use two perpendicular vectors in plane and invariant point to return
/// plane in 'x, y, 0'-type notation
std::string sym_plane(const xtal::Coordinate &v1, const xtal::Coordinate &v2,
                      const xtal::Coordinate &point, SymInfoOptions opt) {
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
  // and situations 'x' is used for 'a' and neither 'y' or 'z' can be used for
  // 'b'
  //
  for (int i = 0; i < 3; ++i) {
    if (!almost_zero(a[i], tol) && !almost_zero(b[i], tol)) {
      b = b - b[i] / a[i] * a;
      break;
    }
  }

  for (int i = 0; i < 3; ++i) {
    if (!almost_zero(a[i], tol) && !almost_zero(b[i], tol)) {
      a = a - a[i] / b[i] * b;
      break;
    }
  }

  // find first non-zero axis coordinate, to choose 'x', 'y', or 'z' as param
  std::string _param = "xyz";
  int min_i = -1;
  for (int i = 0; i < 3; ++i) {
    if (!almost_zero(a[i], tol)) {
      min_i = i;
      break;
    }
  }
  int min_j = -1;
  for (int j = 1; j < 3; ++j) {
    if (!almost_zero(b[j], tol) && j != min_i) {
      min_j = j;
      break;
    }
  }

  // associate lesser of 'x','y','z' with 'a'
  if (min_j < min_i) {
    using std::swap;
    swap(a, b);
    swap(min_j, min_i);
  }

  std::string param1 = _param.substr(min_i, 1);
  std::string param2 = _param.substr(min_j, 1);

  if (opt.coord_type == FRAC) {
    PolyWriter::standardize(a, tol);
    PolyWriter::standardize(b, tol);

    // fractional coords can be scaled to int
    a = scale_to_int(a).cast<double>();
    b = scale_to_int(b).cast<double>();
  }

  // p + a*x + b*y, where p is a constant point,
  //   a, b are vectors of coeffs, and param1 and param2 are variable names
  for (int i = 0; i < 3; ++i) {
    PolyWriter writer(opt);
    writer.add_constant(p[i]);
    writer.add_term(a[i], param1);
    writer.add_term(b[i], param2);

    if (!writer.ss.str().size()) {
      ss << 0;
    } else {
      ss << writer.ss.str();
    }

    if (i != 2) {
      ss << ", ";
    }
  }

  return ss.str();
}

/// \brief Use axis and invariant point to return plane in 'x, y, 0'-type
/// notation
std::string sym_plane(const xtal::Coordinate &axis,
                      const xtal::Coordinate &point, SymInfoOptions opt) {
  double tol = opt.tol;

  Eigen::Vector3d n = axis.const_cart();
  Eigen::Vector3d _v1;
  Eigen::Vector3d _v2;
  n.normalize();

  // get first vector orthogonal to axis
  if (!almost_equal(std::abs(n(0)), 1., tol)) {
    Eigen::Vector3d x;
    x << 1., 0., 0.;
    _v1 = x - x.dot(n) * n;
  } else {
    Eigen::Vector3d y;
    y << 0., 1., 0.;
    _v1 = y - y.dot(n) * n;
  }
  PolyWriter::standardize(_v1, tol);

  // get second vector in symmetry plane
  _v2 = _v1.cross(n);
  PolyWriter::standardize(_v2, tol);

  // use plane defining vectors and invariant point to write the plane equation
  xtal::Coordinate v1(_v1, axis.home(), CART);
  xtal::Coordinate v2(_v2, axis.home(), CART);
  return sym_plane(v1, v2, point, opt);
}

/// \brief Print symmetry symbol to string
///
/// Of the form:
/// \code
/// m 0.25 0.25 0.0
/// \endcode
///
/// Uses prime symbol for time reversal symmetry.
std::string to_brief_unicode(const SymInfo &info, SymInfoOptions opt) {
  std::stringstream stream;
  char term(0);
  stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
  Eigen::IOFormat format(opt.prec, opt.prec + 1);

  // print rotation or screw integer, and include trailing space
  auto print_r = [&](bool inversion) {
    int r = std::lround(360. / std::min(info.angle, 360. - info.angle));
    if (r == 2) {
      if (inversion) stream << "-";
      stream << boost::lexical_cast<std::string>(r);
    } else if (info.angle < 180.) {
      if (inversion) stream << "-";
      stream << boost::lexical_cast<std::string>(r);
      // add superscript '+'
      stream << "\u207A";
    } else {
      if (inversion) stream << "-";
      stream << boost::lexical_cast<std::string>(r);
      // add superscript '-'
      stream << "\u207B";
    }
    if (info.time_reversal) stream << "\u2032";
    stream << " ";
  };

  auto print_coord = [&](const xtal::Coordinate &coord) {
    coord.print(stream, opt.coord_type, term, format);
  };

  switch (info.op_type) {
    case symmetry_type::identity_op:
      stream << "1";
      if (info.time_reversal) stream << "\u2032";
      break;

    case symmetry_type::mirror_op:
      stream.setf(std::ios::showpoint);
      stream << "m";
      if (info.time_reversal) stream << "\u2032";
      stream << " ";
      stream << sym_plane(info.axis, info.location, opt);
      break;

    case symmetry_type::glide_op:
      stream.setf(std::ios::showpoint);
      stream << "g";
      if (info.time_reversal) stream << "\u2032";
      stream << " (";
      print_coord(info.screw_glide_shift);  // glide shift
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
      print_coord(info.screw_glide_shift);  // screw shift
      stream << ") " << sym_line(info.axis, info.location, opt);
      break;

    case symmetry_type::inversion_op:
      stream << "-1";
      if (info.time_reversal) stream << "\u2032";
      stream << " ";
      print_coord(info.location);  // inversion point
      break;

    case symmetry_type::rotoinversion_op:
      print_r(true);
      stream << sym_line(info.axis, info.location, opt);
      stream << "; ";
      print_coord(info.location);  // invariant point
      break;

    case symmetry_type::invalid_op:
      stream << "<error invalid>";
      break;
  }

  return stream.str();
}

/// \brief Print SymInfo to string
std::string description(const SymOp &op, const xtal::Lattice &lat,
                        SymInfoOptions opt) {
  return to_string(SymInfo(op, lat), opt);
}

/// \brief Print SymGroup info
///
/// Uses indexing starting at 1
void description(Log &log, const SymGroup &g, const xtal::Lattice &lat,
                 SymInfoOptions opt) {
  bool print_master_index = false;
  Index i = 0;
  for (const auto &op : g) {
    if (op.has_valid_master() && i != op.index()) {
      print_master_index = true;
      continue;
    }
    ++i;
  }

  i = 1;
  for (const auto &op : g) {
    log << log.indent_str() << i;
    if (print_master_index) {
      log << " (" << op.index() + 1 << ")";
    }
    log << ": " << description(op, lat, opt) << std::endl;
    if (opt.print_matrix_tau) {
      print_matrix_tau_col(log, op, opt.prec);
    }
    ++i;
  }
}

/// \brief Print SymInfo to brief string
std::string brief_description(const SymOp &op, const xtal::Lattice &lat,
                              SymInfoOptions opt) {
  return to_brief_unicode(SymInfo(op, lat), opt);
}

/// \brief Print SymGroup with brief strings
///
/// Uses indexing starting at 1
/// Format:
///     <group index> (master group index): <symmetry operation>
/// The master group index is only printed when some differ from the group
/// index.
void brief_description(Log &log, const SymGroup &g, const xtal::Lattice &lat,
                       SymInfoOptions opt) {
  bool print_master_index = false;
  Index i = 0;
  for (const auto &op : g) {
    if (op.has_valid_master() && i != op.index()) {
      print_master_index = true;
      continue;
    }
    ++i;
  }

  i = 1;
  for (const auto &op : g) {
    log << log.indent_str() << i;
    if (print_master_index) {
      log << " (" << op.index() + 1 << ")";
    }
    log << ": " << brief_description(op, lat, opt) << std::endl;
    if (opt.print_matrix_tau) {
      print_matrix_tau_col(log, op, opt.prec);
    }
    ++i;
  }
}

}  // namespace CASM
