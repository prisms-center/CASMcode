#include "casm/crystallography/Coordinate.hh"

#include "casm/misc/CASM_math.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {
  namespace Coordinate_impl {
    bool verbose::val = false;
  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  /**
   *
   */
  //********************************************************************

  Coordinate &Coordinate::operator +=(const Coordinate &RHS) {
    cart() += RHS.cart();
    return *this;
  }

  Coordinate &Coordinate::operator -=(const Coordinate &RHS) {
    cart() -= RHS.cart();
    return *this;

  }

  //********************************************************************
  /**
   *
   */
  //********************************************************************

  Coordinate Coordinate::operator-() const {
    return Coordinate(-frac(), home(), FRAC);
  }

  //********************************************************************
  /**
   *
   */
  //********************************************************************

  bool Coordinate::operator ==(const Coordinate &RHS) const {
    return almost_equal(m_cart_coord, RHS.m_cart_coord);
  }


  //********************************************************************

  ///These compares exist to make interface consistent with site
  bool Coordinate::compare(const Coordinate &RHS, double compare_tol) const {
    return (compare_type(RHS)) && (min_dist(RHS) < compare_tol);

  }

  //********************************************************************

  bool Coordinate::compare(const Coordinate &RHS, Coordinate &shift, double compare_tol) const {
    return (compare_type(RHS)) && (min_dist(RHS + shift) < compare_tol);
  }

  //********************************************************************

  bool Coordinate::compare_type(const Coordinate &RHS)const {
    return true;
  }

  //********************************************************************
  /**
   * Applies symmetry to a coordinate.
   *
   * Overloads the * operator to apply symmetry to a coordinate in
   * the current mode unless otherwise specified. Returns the
   * transformed coordinate.
   */
  //********************************************************************

  Coordinate operator*(const SymOp &LHS, const Coordinate &RHS) {
    Coordinate tcoord(RHS);
    return tcoord.apply_sym(LHS);
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
  /**
   * Finds distance between two coordinates
   */
  //********************************************************************

  double Coordinate::dist(const Coordinate &neighbor) const {

    return (cart() - neighbor.cart()).norm();
  }


  //********************************************************************
  /**
   * Finds minimum distance from any periodic image of a coordinate to any
   * periodic image of a neighboring coordinate
   */
  //********************************************************************

  double Coordinate::min_dist(const Coordinate &neighbor) const {
    vector_type tfrac(frac() - neighbor.frac());

    for(int i = 0; i < 3; i++)
      tfrac(i) -= round(tfrac(i));

    return (home().lat_column_mat() * tfrac).norm();
  }

  //********************************************************************
  /**
   * Finds minimum distance from any periodic image of a coordinate to any
   * periodic image of a neighboring coordinate.  Also passes the
   * calculated shift in the argument.
   */
  //********************************************************************


  double Coordinate::min_dist(const Coordinate &neighbor, Coordinate &shift) const {
    shift.m_home = m_home;

    shift.m_frac_coord = (frac() - neighbor.frac());

    for(int i = 0; i < 3; i++)
      shift.m_frac_coord(i) -= round(shift.m_frac_coord(i));

    shift._update_cart();

    return length(shift.const_cart());
  };

  //********************************************************************
  /**
   * Finds minimum distance from any periodic image of a coordinate to any
   * periodic image of a neighboring coordinate
   */
  //********************************************************************

  double Coordinate::min_dist2(const Coordinate &neighbor, const Eigen::Ref<const Eigen::Matrix3d> &metric) const {
    vector_type tfrac(frac() - neighbor.frac());

    for(int i = 0; i < 3; i++)
      tfrac(i) -= round(tfrac(i));

    return (home().lat_column_mat() * tfrac).dot(metric * (home().lat_column_mat() * tfrac));
  }

  //********************************************************************
  /**
   * Applies symmetry to a coordinate.
   */
  //********************************************************************

  Coordinate &Coordinate::apply_sym(const SymOp &op) {
    cart() = op.matrix() * const_cart() + op.tau();

    return *this;
  }

  //********************************************************************
  /**
   * Applies symmetry to a coordinate, but doesn't apply translation
   */
  //********************************************************************

  Coordinate &Coordinate::apply_sym_no_trans(const SymOp &op) {
    cart() = op.matrix() * const_cart();
    return *this;
  }

  //********************************************************************
  //John G. You decide which coordinates (FRAC or CART) to keep the same when you set a new lattice.
  //CART: Adds empty space around atoms (This is what previous set_lattice does)
  //FRAC: Shear atoms with vectors
  //inline
  void Coordinate::set_lattice(const Lattice &new_lat, COORD_TYPE invariant_mode) {
    m_home = &new_lat;

    if(invariant_mode == CART)
      _update_frac();
    else if(invariant_mode == FRAC)
      _update_cart();

    return;
  }


  //********************************************************************
  /**
   * Within: if a point is outiside the cell defined by Lattice home
   * then map that point back into the cell via translation by lattice vectors
   *
   */
  //********************************************************************
  bool Coordinate::within() {
    if(PERIODICITY_MODE::IS_LOCAL()) return true;

    bool is_within = true;
    double tshift;
    for(int i = 0; i < 3; i++) {
      tshift = floor(m_frac_coord[i] + 1E-6);
      if(std::abs(tshift) > TOL) {
        is_within = false;
        m_frac_coord[i] -= tshift;
      }
    }
    if(!is_within)
      _update_cart();
    return is_within;
  };
  //***********************************

  bool Coordinate::is_within() const {

    double tshift;
    for(int i = 0; i < 3; i++) {
      tshift = floor(m_frac_coord[i] + 1E-6);
      if(std::abs(tshift) > TOL) {
        return false;
      }
    }
    return true;
  };

  //********************************************************************

  bool Coordinate::within(Coordinate &translation) {
    translation.m_home = m_home;
    if(PERIODICITY_MODE::IS_LOCAL()) return true;

    bool is_within = true;

    for(int i = 0; i < 3; i++) {
      translation.m_frac_coord[i] = -floor(m_frac_coord[i] + 1E-6);
      if(std::abs(translation.m_frac_coord[i]) > TOL) {
        is_within = false;
        m_frac_coord[i] += translation.m_frac_coord[i];
      }
    }

    translation._update_cart();
    return is_within;
  };

  //********************************************************************
  /**
   * Checks to see if coordinate is within voronoi cell of its home lattice.
   */
  //********************************************************************

  int Coordinate::voronoi_number() const {
    return voronoi_number(home());
  }

  //********************************************************************
  /**
   * Checks to see if coordinate is within voronois cell of specified lattice
   */
  //********************************************************************

  int Coordinate::voronoi_number(const Lattice &cell) const {
    return cell.voronoi_number(m_cart_coord);
  }

  //********************************************************************
  /**
   * Applies Lattice translations until coordinate is within voronois cell of its home lattice.
   * **** HAS NOT BEEN PROPERLY TESTED
   */
  //********************************************************************

  bool Coordinate::voronoi_within() {
    bool was_within(true);

    while(voronoi_number() < 0) {
      was_within = false;
      Coordinate tcoord(home().max_voronoi_vector(m_cart_coord), home(), CART);
      frac() = (round(const_cart().dot(tcoord.const_cart()) / 2) * scale_to_int(tcoord.const_frac())).cast<double>();
    }

    return was_within;
  }

  //********************************************************************
  /**
   * Checks to see if coordinate describes shift by a general lattice vector l*V1+m*V2+n*V3, where l, m, n are integer
   */
  //********************************************************************

  bool Coordinate::is_lattice_shift() {

    //If mode is local, return true only if coordinate describes origin
    if(PERIODICITY_MODE::IS_LOCAL())
      return std::abs(m_frac_coord[0]) < TOL && std::abs(m_frac_coord[1]) < TOL && std::abs(m_frac_coord[2]) < TOL;

    return (std::abs(m_frac_coord[0] - round(m_frac_coord[0])) < TOL
            && std::abs(m_frac_coord[1] - round(m_frac_coord[1])) < TOL
            && std::abs(m_frac_coord[2] - round(m_frac_coord[2])) < TOL);
  }


  //********************************************************************
  /**
   * Read/write coordinate to json
   */
  //********************************************************************

  jsonParser &Coordinate::to_json(jsonParser &json) const {
    json.put_obj();

    // mutable Vector3< double > coord[2];
    json["FRAC"].put_array();
    json["FRAC"].push_back(m_frac_coord[0]);
    json["FRAC"].push_back(m_frac_coord[1]);
    json["FRAC"].push_back(m_frac_coord[2]);

    json["CART"].put_array();
    json["CART"].push_back(m_cart_coord[0]);
    json["CART"].push_back(m_cart_coord[1]);
    json["CART"].push_back(m_cart_coord[2]);

    // mutable int basis_ind;
    json["basis_ind"] = basis_ind();

    return json;
  }

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
  };

  jsonParser &to_json(const Coordinate &coord, jsonParser &json) {
    return coord.to_json(json);
  }

  void from_json(Coordinate &coord, const jsonParser &json) {
    coord.from_json(json);
  };


  //********************************************************************
  /**
   * Allows direct printing of Coordinate with the << operator.
   */
  //********************************************************************

  std::ostream &operator << (std::ostream &stream, const Coordinate &coord) {
    coord.print(stream);
    return stream;
  }

};

