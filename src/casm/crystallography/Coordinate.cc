#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/SymTools.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {
  namespace xtal {

    Coordinate::Coordinate(const Eigen::Ref<const Coordinate::vector_type> &init_vec,
                           const Lattice &init_home,
                           COORD_TYPE mode)
      : m_home(&init_home) {
      if(mode == FRAC)
        _set_frac(init_vec);
      if(mode == CART)
        _set_cart(init_vec);
    }

    //********************************************************************

    Coordinate::Coordinate(double _x, double _y, double _z, const Lattice &init_home, COORD_TYPE mode)
      : m_home(&init_home) {
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

    bool Coordinate::almost_equal(const Coordinate &RHS) const {
      return dist(RHS) < lattice().tol();
    }

    //********************************************************************

    bool Coordinate::compare(const Coordinate &RHS) const {
      return (compare_type(RHS)) && (min_dist(RHS) < lattice().tol());

    }

    //********************************************************************

    bool Coordinate::compare(const Coordinate &RHS, Coordinate &translation) const {
      return (compare_type(RHS)) && (min_dist(RHS + translation) < lattice().tol());
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

    void Coordinate::print(std::ostream &stream, char term, Eigen::IOFormat format) const {
      print(stream, COORD_MODE::CHECK(), term, format);
    }

    //********************************************************************

    void _formatted_print(std::ostream &stream, Eigen::Vector3d vec, COORD_TYPE mode, char term, Eigen::IOFormat format) {
      //    stream.precision(prec);
      //    stream.width(prec + pad);
      //    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

      stream << vec.transpose().format(format);
      if(term) stream << term;
      return;
    }

    //********************************************************************
    void Coordinate::print(std::ostream &stream, COORD_TYPE mode, char term, Eigen::IOFormat format) const {
      Eigen::Vector3d vec;

      if(mode == CART)
        vec = const_cart();
      else if(mode == FRAC)
        vec = const_frac();
      _formatted_print(stream, vec, mode, term, format);
      return;
    }

    //********************************************************************

    void Coordinate::print_axis(std::ostream &stream, COORD_TYPE mode, char term, Eigen::IOFormat format) const {

      Eigen::Vector3d vec;
      if(mode == CART)
        vec = const_cart().normalized();
      else if(mode == FRAC)
        vec = const_frac().normalized();
      _formatted_print(stream, vec, mode, term, format);
      return;
    }

    //********************************************************************

    double Coordinate::dist(const Coordinate &neighbor) const {
      return (cart() - neighbor.cart()).norm();
    }

    double Coordinate::min_dist(const Coordinate &neighbor) const {
      return this->min_translation(neighbor).const_cart().norm();
    }

    Coordinate Coordinate::min_translation(const Coordinate &neighbor) const {
      Coordinate translation = (*this) - neighbor;
      translation.frac() -= lround(translation.const_frac()).cast<double>();
      return translation;
    }

    double Coordinate::robust_min_dist(const Coordinate &neighbor) const {
      Coordinate translation = this->min_translation(neighbor);
      double fast_result = translation.const_cart().norm();

      if(fast_result < (home().inner_voronoi_radius() + TOL))
        return fast_result;

      translation.voronoi_within();
      return translation.const_cart().norm();
    }

    double Coordinate::min_dist2(const Coordinate &neighbor, const Eigen::Ref<const Eigen::Matrix3d> &metric) const {
      vector_type tfrac(frac() - neighbor.frac());
      tfrac -= lround(tfrac).cast<double>();

      return (home().lat_column_mat() * tfrac).dot(metric * (home().lat_column_mat() * tfrac));
    }

    //********************************************************************

    /* Coordinate &Coordinate::apply_sym(const SymOp &op) { */
    /*   cart() = get_matrix(op) * this->const_cart() + get_translation(op); */

    /*   return *this; */
    /* } */

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
    // Applies symmetry to a coordinate.
    //
    // Overloads the * operator to apply symmetry to a coordinate
    // using its Cartesian representation
    //********************************************************************

    Coordinate operator*(const SymOp &LHS, const Coordinate &RHS) {
      Coordinate tcoord(RHS);
      sym::apply(LHS, tcoord);
      return tcoord;
    }

    //********************************************************************
    // Allows direct printing of Coordinate with the << operator.
    //********************************************************************

    std::ostream &operator << (std::ostream &stream, const Coordinate &coord) {
      coord.print(stream);
      return stream;
    }

  }

  namespace sym {
    xtal::Coordinate &apply(const xtal::SymOp &op, xtal::Coordinate &mutating_coord) {
      mutating_coord.cart() = get_matrix(op) * mutating_coord.const_cart() + get_translation(op);
      return mutating_coord;
    }

    xtal::Coordinate copy_apply(const xtal::SymOp &op, xtal::Coordinate coord) {
      apply(op, coord);
      return coord;
    }
  } // namespace sym
}

