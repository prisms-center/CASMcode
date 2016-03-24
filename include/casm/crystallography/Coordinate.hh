#ifndef COORDINATE_HH
#define COORDINATE_HH

#include <iostream>
#include <cmath>
#include <cassert>

#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/container/LinearAlgebra.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class Lattice;
  class SymOp;

  namespace Coordinate_impl {

    class FracCoordinate;
    class FracCoordinateComponent;
    class CartCoordinate;
    class CartCoordinateComponent;

  }


  class Coordinate {
  public:
    typedef Eigen::Vector3d vector_type;
    typedef vector_type::Index size_type;

    /**NOTE: Coordinate does not have a default constructor
       e.g: this is not allowed-> Coordinate() : home(nullptr) { is_current[FRAC]=false; is_current[CART]=false;}; **/

    ///Minimal constructor only takes a lattice
    explicit Coordinate(const Lattice &init_home) :
      m_home(&init_home),
      m_frac_coord(vector_type::Zero()),
      m_cart_coord(vector_type::Zero()),
      m_basis_ind(-1) {

    }

    Coordinate(const vector_type &init_vec, const Lattice &init_home, COORD_TYPE mode);

    Coordinate(double _x, double _y, double _z, const Lattice &init_home, COORD_TYPE mode);

    /// \brief Set the fractional coordinate vector
    Coordinate_impl::FracCoordinate frac();

    /// \brief const Access the fractional coordinate vector
    inline
    const vector_type &frac() const {
      return m_frac_coord;
    }

    /// \brief user override to force const Access the fractional coordinate vector
    inline
    const vector_type &const_frac() const {
      return m_frac_coord;
    }

    /// \brief Set a component of the fractional coordinate vector
    Coordinate_impl::FracCoordinateComponent frac(size_type index);

    /// \brief const Access a component of the fractional coordinate vector
    const double &frac(size_type index) const;

    /// \brief Set Cartesian coordinate vector and update fractional coordiante vector
    Coordinate_impl::CartCoordinate cart();

    /// \brief const Access the Cartesian coordinate vector
    inline
    const vector_type &cart() const {
      return m_cart_coord;
    }

    /// \brief user override to force const Access the Cartesian coordinate vector
    inline
    const vector_type &const_cart() const {
      return m_cart_coord;
    }

    /// \brief Set a component of the Cartesian coordinate vector
    Coordinate_impl::CartCoordinateComponent cart(size_type index);

    /// \brief const Access a component of the Cartesian coordinate vector
    const double &cart(size_type index) const;

    Coordinate &operator+=(const Coordinate &RHS);
    Coordinate &operator-=(const Coordinate &RHS);

    Coordinate operator-() const;

    bool operator==(const Coordinate &RHS) const; //Ivy
    bool operator!=(const Coordinate &RHS) const {
      return !(*this == RHS);
    }


    ///These compares exist to make interface consistent with site
    bool compare(const Coordinate &RHS, double tol = TOL) const;
    bool compare(const Coordinate &RHS, Coordinate &shift, double tol = TOL) const;
    bool compare_type(const Coordinate &RHS)const;

    ///Map coordinate into the unit cell using a lattice translation
    bool within();

    ///Same as within(), but lattice translation is stored in Coordinate translation
    bool within(Coordinate &translation);

    ///Checks to see if coordinate is in the unit cell, but does not translate it
    bool is_within() const;

    int voronoi_number() const;
    int voronoi_number(const Lattice &cell) const;

    ///Map coordinate into the voronoi cell using a lattice translation
    bool voronoi_within();

    ///Same as voronoi_within(), but lattice translation is stored in Coordinate translation
    bool voronoi_within(Coordinate &translation);

    ///Checks to see if coordinate is at a lattice translation with respect to the origin
    bool is_lattice_shift();

    Coordinate get_normal_vector(Coordinate coord2, Coordinate coord3);

    ///update the home lattice of a coordinate, keeping original cartesian representation
    void set_lattice(const Lattice &new_lat);
    ///update the home lattice of a coordinate, keeping representation specified mode
    void set_lattice(const Lattice &new_lat, COORD_TYPE mode); //John G, use to specify whether to keep CART or FRAC the same, when assigning a new lattice

    void set_basis_ind(Index _basis_ind) {
      m_basis_ind = _basis_ind;
    }

    Index basis_ind() const {
      return m_basis_ind;
    }

    ///Check the home lattice of the coordinate
    const Lattice &home() const {
      assert(m_home && "Coordinate doesn't have valid home lattice");
      return *m_home;
    }

    //term is terminal character, prec is precision, pad is field width - precision  (should be greater than 3)
    void read(std::istream &stream);
    void read(std::istream &stream, COORD_TYPE mode); //not yet implemented
    void print(std::ostream &stream, COORD_TYPE mode, char term = 0, int prec = 7, int pad = 5) const;
    void print(std::ostream &stream, char term = 0, int prec = 7, int pad = 5) const;

    double dist(const Coordinate &neighbor) const;

    ///Returns distance from periodic image of neighbor that is closest
    ///Is unsafe if min_dist is comparable to half a lattice vector in length
    double min_dist(const Coordinate &neighbor) const;
    double min_dist(const Coordinate &neighbor, Coordinate &shift)const; //Added by Ivy 11/05/12

    ///Finds same shift as min_dist but returns shift(CART).transpose()*metric*shift(CART)
    double min_dist2(const Coordinate &neighbor, const Eigen::Ref<const Eigen::Matrix3d> &metric) const;

    ///Transform coordinate by symmetry operation (including translation)
    Coordinate &apply_sym(const SymOp &op); //AAB
    ///Transform coordinate by symmetry operation (without translation)
    Coordinate &apply_sym_no_trans(const SymOp &op); //AAB
    //static bool SelfTest();

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);

  private:
    inline
    void _update_cart() {
      m_cart_coord = home().lat_column_mat() * m_frac_coord;
    }

    inline
    void _update_frac() {
      m_frac_coord = home().inv_lat_column_mat() * m_cart_coord;
    }

    inline
    void _set_frac(const Eigen::Ref<const vector_type> &f) {
      m_frac_coord = f;
      _update_cart();
    }


    inline
    void _set_frac(size_type ind, double val) {
      m_frac_coord[ind] = val;
      _update_cart();
    }

    inline
    void _set_cart(const Eigen::Ref<const vector_type> &c) {
      m_cart_coord = c;
      _update_frac();
    }

    inline
    void _set_cart(size_type ind, double val) {
      m_cart_coord[ind] = val;
      _update_frac();
    }

    friend Coordinate_impl::FracCoordinate;
    friend Coordinate_impl::CartCoordinate;
    friend Coordinate_impl::FracCoordinateComponent;
    friend Coordinate_impl::CartCoordinateComponent;

    Lattice const *m_home;

    vector_type m_frac_coord, m_cart_coord;

    Index m_basis_ind;
  };

  jsonParser &to_json(const Coordinate &value, jsonParser &json);

  // Lattice must be set already
  void from_json(Coordinate &value, const jsonParser &json);

  Coordinate operator*(const SymOp &LHS, const Coordinate &RHS); //AAB
  Coordinate operator+(const Coordinate &LHS, const Coordinate &RHS) {
    return Coordinate(LHS) += RHS;
  }

  Coordinate operator-(const Coordinate &LHS, const Coordinate &RHS) {
    return Coordinate(LHS) -= RHS;
  }

  namespace Coordinate_impl {

    /// \brief A class to enable vector assignment to the fractional vector of a Coordinate
    ///
    /// Typically only used indirectly as a temporary when performing
    /// \code
    /// Coordinate coord;
    /// coord.frac() = Coordinate::vector_type(a,b,c);
    /// \endcode
    ///
    class FracCoordinate {

    public:

      explicit FracCoordinate(Coordinate &coord) :
        m_coord(&coord) {}


      FracCoordinate &operator=(const Eigen::Ref<const Coordinate::vector_type> &vec) {
        m_coord->_set_frac(vec);
        return *this;
      }

      FracCoordinate &operator+=(const Eigen::Ref<const Coordinate::vector_type> &vec) {
        (m_coord->m_frac_coord) += vec;
        m_coord->_update_cart();
        return *this;
      }

      FracCoordinate &operator-=(const Eigen::Ref<const Coordinate::vector_type> &vec) {
        (m_coord->m_frac_coord) -= vec;
        m_coord->_update_cart();
        return *this;
      }

      FracCoordinate &operator*=(double val) {
        (m_coord->m_frac_coord) *= val;
        (m_coord->m_cart_coord) *= val;
        return *this;
      }

      FracCoordinate &operator/=(double val) {
        (m_coord->m_frac_coord) /= val;
        (m_coord->m_cart_coord) /= val;
        return *this;
      }


      operator const Eigen::MatrixBase<Eigen::Vector3d> &() const {
        return m_coord->m_frac_coord;
      }

      operator const Eigen::Vector3d &() const {
        return m_coord->m_frac_coord;
      }

      operator Eigen::Ref<const Eigen::Vector3d> () const {
        return m_coord->m_frac_coord;
      }

    private:

      Coordinate *m_coord;
    };

    /// \brief A class to enable assignment to a component of the fractional vector of a Coordinate
    ///
    /// Typically only used indirectly as a temporary when performing
    /// \code
    /// Coordinate coord;
    /// double a, b, c;
    /// coord.frac(0) = a;
    /// coord.frac(1) = b;
    /// coord.frac(2) = c;
    /// \endcode
    ///
    class FracCoordinateComponent {

    public:

      explicit FracCoordinateComponent(Coordinate &coord, Coordinate::size_type index) :
        m_coord(&coord), m_index(index) {}

      FracCoordinateComponent &operator=(double val) {
        m_coord->_set_frac(m_index, val);
        return *this;
      }

      FracCoordinateComponent &operator+=(double val) {
        m_coord->_set_frac(m_index, m_coord->m_frac_coord(m_index) + val);
        return *this;
      }

      FracCoordinateComponent &operator-=(double val) {
        m_coord->_set_frac(m_index, m_coord->m_frac_coord(m_index) - val);
        return *this;
      }

      FracCoordinateComponent &operator*=(double val) {
        m_coord->_set_frac(m_index, m_coord->m_frac_coord(m_index)*val);
        return *this;
      }

      FracCoordinateComponent &operator/=(double val) {
        m_coord->_set_frac(m_index, m_coord->m_frac_coord(m_index) / val);
        return *this;
      }

      operator const double &() const {
        return m_coord->m_frac_coord(m_index);
      }

    private:

      Coordinate *m_coord;
      Coordinate::size_type m_index;
    };

    /// \brief A class to enable vector assignment to the Cartesian vector of a Coordinate
    ///
    /// Typically only used indirectly as a temporary when performing
    /// \code
    /// Coordinate coord;
    /// coord.cart() = Coordinate::vector_type(a,b,c);
    /// \endcode
    ///
    class CartCoordinate {

    public:

      explicit CartCoordinate(Coordinate &coord) :
        m_coord(&coord) {}

      CartCoordinate &operator=(const Eigen::Ref<const Coordinate::vector_type> &vec) {
        m_coord->_set_cart(vec);
        return *this;
      }

      CartCoordinate &operator+=(const Eigen::Ref<const Coordinate::vector_type> &vec) {
        (m_coord->m_cart_coord) += vec;
        m_coord->_update_frac();
        return *this;
      }

      CartCoordinate &operator-=(const Eigen::Ref<const Coordinate::vector_type> &vec) {
        (m_coord->m_cart_coord) -= vec;
        m_coord->_update_frac();
        return *this;
      }

      CartCoordinate &operator*=(double val) {
        (m_coord->m_frac_coord) *= val;
        (m_coord->m_cart_coord) *= val;
        return *this;
      }

      CartCoordinate &operator/=(double val) {
        (m_coord->m_frac_coord) /= val;
        (m_coord->m_cart_coord) /= val;
        return *this;
      }

      operator const Eigen::MatrixBase<Eigen::Vector3d> &() const {
        return m_coord->m_cart_coord;
      }

      operator const Eigen::Vector3d &() const {
        return m_coord->m_cart_coord;
      }

      operator Eigen::Ref<const Eigen::Vector3d> () const {
        return m_coord->m_cart_coord;
      }


    private:

      Coordinate *m_coord;
    };

    /// \brief A class to enable assignment to a component of the Cartesian vector of a Coordinate
    ///
    /// Typically only used indirectly as a temporary when performing
    /// \code
    /// Coordinate coord;
    /// double a, b, c;
    /// coord.cart(0) = a;
    /// coord.cart(1) = b;
    /// coord.cart(2) = c;
    /// \endcode
    ///
    class CartCoordinateComponent {

    public:

      explicit CartCoordinateComponent(Coordinate &coord, Coordinate::size_type index) :
        m_coord(&coord), m_index(index) {}

      CartCoordinateComponent &operator=(double val) {
        m_coord->_set_cart(m_index, val);
        return *this;
      }

      CartCoordinateComponent &operator+=(double val) {
        m_coord->_set_cart(m_index, m_coord->m_cart_coord(m_index) + val);
        return *this;
      }

      CartCoordinateComponent &operator-=(double val) {
        m_coord->_set_cart(m_index, m_coord->m_cart_coord(m_index) - val);
        return *this;
      }

      CartCoordinateComponent &operator*=(double val) {
        m_coord->_set_cart(m_index, m_coord->m_cart_coord(m_index)*val);
        return *this;
      }

      CartCoordinateComponent &operator/=(double val) {
        m_coord->_set_cart(m_index, m_coord->m_cart_coord(m_index) / val);
        return *this;
      }

      operator const double &() const {
        return m_coord->m_cart_coord(m_index);
      }

    private:

      Coordinate *m_coord;
      Coordinate::size_type m_index;
    };
  }

}

namespace std {
  template<>
  struct is_floating_point<CASM::Coordinate_impl::FracCoordinateComponent> {
    static const bool value = true;
  };

  template<>
  struct is_floating_point<CASM::Coordinate_impl::CartCoordinateComponent> {
    static const bool value = true;
  };
}
#endif
