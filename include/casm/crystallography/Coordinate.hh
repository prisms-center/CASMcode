#ifndef COORDINATE_HH
#define COORDINATE_HH

#include <iostream>
#include <cmath>
#include <cassert>

#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/container/LinearAlgebra.hh"

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
  protected:

    Lattice const *home;

  private:
    Eigen::Vector3d m_frac_coord, m_cart_coord;

    Index m_basis_ind;
  public:

    /**NOTE: Coordinate does not have a default constructor
       e.g: this is not allowed-> Coordinate() : home(nullptr) { is_current[FRAC]=false; is_current[CART]=false;}; **/

    ///Minimal constructor only takes a lattice
    explicit Coordinate(const Lattice &init_home) : home(&init_home), m_basis_ind(-1) {

    };

    Coordinate(const Eigen::Vector3d &init_vec, const Lattice &init_home, COORD_TYPE mode) : home(&init_home), m_basis_ind(-1) {
      if(mode == FRAC)
        frac() = init_vec;
      if(mode == CART)
        cart() = init_vec;
    };

    /// \brief Set the fractional coordinate vector
    Coordinate_impl::FracCoordinate frac();

    /// \brief const Access the fractional coordinate vector
    const Eigen::Vector3d &frac() const;

    /// \brief Set a component of the fractional coordinate vector
    Coordinate_impl::FracCoordinateComponent frac(size_type index);

    /// \brief const Access a component of the fractional coordinate vector
    const double &frac(size_type index) const;



    /// \brief Set Cartesian coordinate vector and update fractional coordiante vector
    Coordinate_impl::CartCoordinate cart();

    /// \brief const Access the Cartesian coordinate vector
    const Eigen::Vector3d &cart() const;

    /// \brief Set a component of the Cartesian coordinate vector
    Coordinate_impl::CartCoordinateComponent cart(size_type index);

    /// \brief const Access a component of the Cartesian coordinate vector
    const double &cart(size_type index) const;

    Coordinate &operator +=(const Coordinate &RHS);
    Coordinate &operator -=(const Coordinate &RHS);

    Coordinate operator -() const;

    bool operator ==(const Coordinate &RHS) const; //Ivy
    bool operator !=(const Coordinate &RHS) const {
      return !(*this == RHS);
    };


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

    double &operator[](int ind);  //Access element of current coordinate mode

    void set_basis_ind(Index _basis_ind) {
      m_basis_ind = _basis_ind;
    };

    Index basis_ind() const {
      return m_basis_ind;
    };

    // pseudoconstant access of element in current coordinate mode
    // can't be const, since it calls calc_coord()
    double get(int ind) const;
    double get(int ind, COORD_TYPE mode) const;
    double &at(int ind);
    double &at(int ind, COORD_TYPE mode);

    ///Check the home lattice of the coordinate
    Lattice const *get_home() const {
      return home;
    };

    ///Cast as vector in the current coordinate mode
    operator Eigen::Vector3d();


    /**Retur vector in the current coordinate mode
       example: Coordinate my_coord(prim);
                Eigen::Vector3d my_cart_coord = my_coord(CART); **/
    Eigen::Vector3d &operator()();  //Cast as vector in the current coordinate mode
    Eigen::Vector3d &operator()(COORD_TYPE mode);

    const Eigen::Vector3d &operator()() const;  //Cast as vector in the current coordinate mode
    const Eigen::Vector3d &operator()(COORD_TYPE mode) const;

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
    double min_dist2(const Coordinate &neighbor, const Matrix3<double> &metric) const;

    ///Transform coordinate by symmetry operation (including translation)
    Coordinate &apply_sym(const SymOp &op); //AAB
    ///Transform coordinate by symmetry operation (without translation)
    Coordinate &apply_sym_no_trans(const SymOp &op); //AAB
    //static bool SelfTest();

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);
  };

  jsonParser &to_json(const Coordinate &value, jsonParser &json);

  // Lattice must be set already
  void from_json(Coordinate &value, const jsonParser &json);

  Coordinate operator*(const SymOp &LHS, const Coordinate &RHS); //AAB

  namespace Coordinate_impl {

    /// \brief A class to enable vector assignment to the fractional vector of a Coordinate
    ///
    /// Typically only used indirectly as a temporary when performing
    /// \code
    /// Coordinate coord;
    /// coord.frac() = Eigen::Vector3d(a,b,c);
    /// \endcode
    ///
    class FracCoordinate {

    public:

      explicit FracCoordinate(Coordinate &coord) :
        m_coord(&coord) {}

      template<typename Derived>
      FracCoordinate &operator=(const Eigen::MatrixBase<Derived> &vec) {
        m_coord->_set_frac(vec);
        return *this;
      }

      template<typename Derived>
      FracCoordinate &operator+=(const Eigen::MatrixBase<Derived> &vec) {
        Eigen::Vector3d tmp = m_coord->m_vec[FRAC] + vec;
        m_coord->_set_frac(tmp);
        return *this;
      }

      template<typename Derived>
      FracCoordinate &operator-=(const Eigen::MatrixBase<Derived> &vec) {
        m_coord->_set_frac(m_coord->m_vec[FRAC] - vec);
        return *this;
      }

      FracCoordinate &operator*=(double val) {
        m_coord->_set_frac(m_coord->m_vec[FRAC]*val);
        return *this;
      }

      FracCoordinate &operator/=(double val) {
        m_coord->_set_frac(m_coord->m_vec[FRAC] / val);
        return *this;
      }

      operator const Eigen::Vector3d &() const {
        return m_coord->m_vec[FRAC];
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
        m_coord->_set_frac(m_index, m_coord->m_vec[FRAC](m_index) + val);
        return *this;
      }

      FracCoordinateComponent &operator-=(double val) {
        m_coord->_set_frac(m_index, m_coord->m_vec[FRAC](m_index) - val);
        return *this;
      }

      FracCoordinateComponent &operator*=(double val) {
        m_coord->_set_frac(m_index, m_coord->m_vec[FRAC](m_index)*val);
        return *this;
      }

      FracCoordinateComponent &operator/=(double val) {
        m_coord->_set_frac(m_index, m_coord->m_vec[FRAC](m_index) / val);
        return *this;
      }

      operator const double &() const {
        return m_coord->m_vec[FRAC](m_index);
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
    /// coord.cart() = Eigen::Vector3d(a,b,c);
    /// \endcode
    ///
    class CartCoordinate {

    public:

      explicit CartCoordinate(Coordinate &coord) :
        m_coord(&coord) {}

      template<typename Derived>
      CartCoordinate &operator=(const Eigen::MatrixBase<Derived> &vec) {
        m_coord->_set_cart(vec);
        return *this;
      }

      template<typename Derived>
      CartCoordinate &operator+=(const Eigen::MatrixBase<Derived> &vec) {
        m_coord->_set_cart(m_coord->m_vec[CART] + vec);
        return *this;
      }

      template<typename Derived>
      CartCoordinate &operator-=(const Eigen::MatrixBase<Derived> &vec) {
        m_coord->_set_cart(m_coord->m_vec[CART] - vec);
        return *this;
      }

      CartCoordinate &operator*=(double val) {
        m_coord->_set_cart(m_coord->m_vec[CART]*val);
        return *this;
      }

      CartCoordinate &operator/=(double val) {
        m_coord->_set_cart(m_coord->m_vec[CART] / val);
        return *this;
      }

      operator const Eigen::Vector3d &() const {
        return m_coord->m_vec[CART];
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
        m_coord->_set_cart(m_index, m_coord->m_vec[CART](m_index) + val);
        return *this;
      }

      CartCoordinateComponent &operator-=(double val) {
        m_coord->_set_cart(m_index, m_coord->m_vec[CART](m_index) - val);
        return *this;
      }

      CartCoordinateComponent &operator*=(double val) {
        m_coord->_set_cart(m_index, m_coord->m_vec[CART](m_index)*val);
        return *this;
      }

      CartCoordinateComponent &operator/=(double val) {
        m_coord->_set_cart(m_index, m_coord->m_vec[CART](m_index) / val);
        return *this;
      }

      operator const double &() const {
        return m_coord->m_vec[CART](m_index);
      }

    private:

      Coordinate *m_coord;
      Coordinate::size_type m_index;
    };
  }

};
#endif
