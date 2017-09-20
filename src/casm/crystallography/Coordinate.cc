#include "casm/crystallography/Coordinate.hh"

#include "casm/misc/CASM_math.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  Coordinate::Coordinate(const Eigen::Vector3d &init_vec,
                         const Lattice &init_home,
                         COORD_TYPE mode)
    : m_home(&init_home),
      m_basis_ind(-1) {
    if(mode == FRAC)
      _set_frac(init_vec);
    if(mode == CART)
      _set_cart(init_vec);
  }

  //********************************************************************

  Coordinate::Coordinate(double _x, double _y, double _z, const Lattice &init_home, COORD_TYPE mode)
    : m_home(&init_home),
      m_basis_ind(-1) {
    if(mode == FRAC) {
      m_frac_coord << _x, _y, _z;
      _update_cart();
    }
    if(mode == CART) {
      m_cart_coord << _x, _y, _z;
      _update_frac();
    }
  }

  //********************************************************************

  Coordinate &Coordinate::operator +=(const Coordinate &RHS) {
    cart() += RHS.cart();
    return *this;
  }

  //********************************************************************

  Coordinate &Coordinate::operator-=(const Coordinate &RHS) {
    cart() -= RHS.cart();
    return *this;

  }

  //********************************************************************

  Coordinate Coordinate::operator-() const {
    return Coordinate(-frac(), home(), FRAC);
  }

  //********************************************************************

  bool Coordinate::operator==(const Coordinate &RHS) const {
    return CASM::almost_equal(m_cart_coord, RHS.m_cart_coord);
  }


  //********************************************************************

  bool Coordinate::almost_equal(const Coordinate &RHS, double tol) const {
    return dist(RHS) < tol;
  }


  //********************************************************************

  bool Coordinate::compare(const Coordinate &RHS, double compare_tol) const {
    return (compare_type(RHS)) && (min_dist(RHS) < compare_tol);

  }

  //********************************************************************

  bool Coordinate::compare(const Coordinate &RHS, Coordinate &translation, double compare_tol) const {
    return (compare_type(RHS)) && (min_dist(RHS + translation) < compare_tol);
  }

  //********************************************************************

  bool Coordinate::compare_type(const Coordinate &RHS)const {
    return true;
  }

  //********************************************************************

  void Coordinate::read(std::istream &stream, COORD_TYPE mode) {
    if(mode == FRAC) {
      stream >> m_frac_coord;
      _update_cart();
    }
    else if(mode == CART) {
      stream >> m_cart_coord;
      _update_frac();
    }

    return;
  }

  //********************************************************************

  void Coordinate::print(std::ostream &stream, char term, int prec, int pad) const {
    print(stream, COORD_MODE::CHECK(), term, prec, pad);
  }

  //********************************************************************
  void Coordinate::print(std::ostream &stream, COORD_TYPE mode, char term, int prec, int pad) const {

    stream.precision(prec);
    stream.width(prec + pad);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

    if(mode == CART)
      stream << const_cart().transpose();
    else if(mode == FRAC)
      stream << const_frac().transpose();
    if(term) stream << term;
    return;
  }

  //********************************************************************

  /// \brief Print normalized vector
  void Coordinate::print_axis(std::ostream &stream, COORD_TYPE mode, char term, int prec, int pad) const {

    stream.precision(prec);
    stream.width(prec + pad);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

    if(mode == CART)
      stream << const_cart().normalized().transpose();
    else if(mode == FRAC)
      stream << const_frac().normalized().transpose();
    if(term) stream << term;
    return;
  }

  //********************************************************************
  // Finds distance between two coordinates

  double Coordinate::dist(const Coordinate &neighbor) const {
    return (cart() - neighbor.cart()).norm();
  }

  //********************************************************************
  // Finds minimum distance from any periodic image of a coordinate to any
  // periodic image of a neighboring coordinate
  //********************************************************************

  double Coordinate::min_dist(const Coordinate &neighbor) const {
    vector_type tfrac(frac() - neighbor.frac());
    tfrac -= lround(tfrac).cast<double>();

    return (home().lat_column_mat() * tfrac).norm();
  }

  //********************************************************************
  // Finds minimum distance from any periodic image of a coordinate to any
  // periodic image of a neighboring coordinate.  Also passes the
  // calculated translation in the argument.
  //********************************************************************


  double Coordinate::min_dist(const Coordinate &neighbor, Coordinate &translation) const {
    translation = (*this) - neighbor;
    translation.frac() -= lround(translation.const_frac()).cast<double>();
    return translation.const_cart().norm();
  }
  //********************************************************************
  // Finds minimum distance from any periodic image of a coordinate to any
  // periodic image of a neighboring coordinate
  //********************************************************************

  double Coordinate::robust_min_dist(const Coordinate &neighbor) const {
    Coordinate translation(home());
    return robust_min_dist(neighbor, translation);
  }

  //********************************************************************
  // Finds minimum distance from any periodic image of a coordinate to any
  // periodic image of a neighboring coordinate.  Also passes the
  // calculated translation in the argument.
  //********************************************************************

  double Coordinate::robust_min_dist(const Coordinate &neighbor, Coordinate &translation) const {
    double fast_result = min_dist(neighbor, translation);

    if(fast_result < (home().inner_voronoi_radius() + TOL))
      return fast_result;

    translation.voronoi_within();
    return translation.const_cart().norm();
  }

  //********************************************************************
  // Finds minimum distance from any periodic image of a coordinate to any
  // periodic image of a neighboring coordinate
  //********************************************************************

  double Coordinate::min_dist2(const Coordinate &neighbor, const Eigen::Ref<const Eigen::Matrix3d> &metric) const {
    vector_type tfrac(frac() - neighbor.frac());
    tfrac -= lround(tfrac).cast<double>();

    return (home().lat_column_mat() * tfrac).dot(metric * (home().lat_column_mat() * tfrac));
  }

  //********************************************************************
  // Applies symmetry to a coordinate.
  //********************************************************************

  Coordinate &Coordinate::apply_sym(const SymOp &op) {
    cart() = op.matrix() * const_cart() + op.tau();

    return *this;
  }

  //********************************************************************
  // Applies symmetry to a coordinate, but doesn't apply translation
  //********************************************************************

  Coordinate &Coordinate::apply_sym_no_trans(const SymOp &op) {
    cart() = op.matrix() * const_cart();
    return *this;
  }

  //********************************************************************
  // Change the home lattice of the coordinate, selecting one representation (either CART or FRAC)
  // that remains invariant
  //
  // invariant_mode == CART: Cartesian coordinates stay the same, and fractional coordinates are updated
  //    Ex: (my_coord.set_lattice(superlattice, CART); // this is how superlattices get filled.
  //
  // invariant_mode == FRAC: Fractional coordinates stay the same, and Cartesian coordinates are updated
  //    Ex:  you can apply a strain by changing the lattice and keepin FRAC invariant
  //            (my_coord.set_lattice(strained_lattice, FRAC);
  //    Ex:  you can apply a rotation by changing the lattice and keeping FRAC invariant
  //            (my_coord.set_lattice(rotated_lattice, FRAC);
  //********************************************************************

  void Coordinate::set_lattice(const Lattice &new_lat, COORD_TYPE invariant_mode) {
    m_home = &new_lat;

    if(invariant_mode == CART)
      _update_frac();
    else if(invariant_mode == FRAC)
      _update_cart();

    return;
  }

  //********************************************************************
  // Within: if a point is outiside the cell defined by Lattice home
  // then map that point back into the cell via translation by lattice vectors
  //********************************************************************

  bool Coordinate::within() {
    if(PERIODICITY_MODE::IS_LOCAL()) return true;

    bool was_within = true;
    double tshift;
    for(int i = 0; i < 3; i++) {
      tshift = floor(m_frac_coord[i] + 1E-6);
      if(!almost_zero(tshift, TOL)) {
        was_within = false;
        m_frac_coord[i] -= tshift;
      }
    }
    if(!was_within)
      _update_cart();
    return was_within;
  }

  //********************************************************************

  bool Coordinate::within(Coordinate &translation) {
    translation.m_home = m_home;
    if(PERIODICITY_MODE::IS_LOCAL()) return true;

    bool was_within = true;

    for(int i = 0; i < 3; i++) {
      translation.m_frac_coord[i] = -floor(m_frac_coord[i] + 1E-6);
      if(!almost_zero(translation.m_frac_coord[i], TOL)) {
        was_within = false;
      }
    }
    translation._update_cart();

    if(!was_within)
      (*this) += translation;
    return was_within;
  }

  //********************************************************************

  bool Coordinate::is_within() const {

    double tshift;
    for(int i = 0; i < 3; i++) {
      tshift = floor(m_frac_coord[i] + 1E-6);
      if(std::abs(tshift) > TOL) {
        return false;
      }
    }
    return true;
  }

  //********************************************************************
  // Checks to see if coordinate is within voronoi cell of its home lattice.
  //********************************************************************

  int Coordinate::voronoi_number() const {
    return voronoi_number(home());
  }

  //********************************************************************
  // Checks to see if coordinate is within voronoi cell of specified lattice
  //********************************************************************

  int Coordinate::voronoi_number(const Lattice &cell) const {
    return cell.voronoi_number(const_cart());
  }

  //********************************************************************
  // Applies Lattice translations until coordinate is within voronoi cell of its home lattice.
  //********************************************************************

  bool Coordinate::voronoi_within() {
    bool was_within(true);
    Eigen::Vector3d lattice_trans;
    while(home().max_voronoi_measure(cart(), lattice_trans) > (1. + TOL)) {
      was_within = false;
      cart() -= lattice_trans;
    }

    return was_within;
  }

  //********************************************************************
  // Applies Lattice translations until coordinate is within voronoi cell of its home lattice.
  //********************************************************************

  bool Coordinate::voronoi_within(Coordinate &translation) {
    translation.m_home = m_home;
    translation.cart() = vector_type::Zero();

    Eigen::Vector3d lattice_trans;
    bool was_within(true);
    while(home().max_voronoi_measure(cart(), lattice_trans) > (1. + TOL)) {
      was_within = false;
      cart() -= lattice_trans;
      translation.cart() -= lattice_trans;
    }

    return was_within;
  }

  //********************************************************************
  // Checks to see if coordinate describes shift by a general lattice vector l*V1+m*V2+n*V3, where l, m, n are integer
  //********************************************************************

  bool Coordinate::is_lattice_shift(double tol) const {

    //If mode is local, return true only if coordinate describes origin
    if(PERIODICITY_MODE::IS_LOCAL())
      return std::abs(m_frac_coord[0]) < tol && std::abs(m_frac_coord[1]) < tol && std::abs(m_frac_coord[2]) < tol;

    return (std::abs(m_frac_coord[0] - round(m_frac_coord[0])) < tol
            && std::abs(m_frac_coord[1] - round(m_frac_coord[1])) < tol
            && std::abs(m_frac_coord[2] - round(m_frac_coord[2])) < tol);
  }


  //********************************************************************
  // Read/write coordinate to json
  //********************************************************************

  jsonParser &Coordinate::to_json(jsonParser &json) const {
    json.put_obj();

    // mutable Vector3< double > coord[2];
    to_json_array(m_frac_coord, json["FRAC"]);
    to_json_array(m_frac_coord, json["CART"]);

    // mutable int basis_ind;
    json["basis_ind"] = basis_ind();

    return json;
  }

  //********************************************************************

  void Coordinate::from_json(const jsonParser &json) {
    // mutable Vector3< double > coord[2];
    m_frac_coord[0] = json["FRAC"][0].get<double>();
    m_frac_coord[1] = json["FRAC"][1].get<double>();
    m_frac_coord[2] = json["FRAC"][2].get<double>();

    m_cart_coord[0] = json["CART"][0].get<double>();
    m_cart_coord[1] = json["CART"][1].get<double>();
    m_cart_coord[2] = json["CART"][2].get<double>();


    // mutable int basis_ind;
    CASM::from_json(m_basis_ind, json["basis_ind"]);
  }

  //********************************************************************

  jsonParser &to_json(const Coordinate &coord, jsonParser &json) {
    return coord.to_json(json);
  }

  //********************************************************************

  void from_json(Coordinate &coord, const jsonParser &json) {
    coord.from_json(json);
  }

  //********************************************************************
  // Applies symmetry to a coordinate.
  //
  // Overloads the * operator to apply symmetry to a coordinate
  // using its Cartesian representation
  //********************************************************************

  Coordinate operator*(const SymOp &LHS, const Coordinate &RHS) {
    Coordinate tcoord(RHS);
    return tcoord.apply_sym(LHS);
  }

  //********************************************************************
  // Allows direct printing of Coordinate with the << operator.
  //********************************************************************

  std::ostream &operator << (std::ostream &stream, const Coordinate &coord) {
    coord.print(stream);
    return stream;
  }

}

