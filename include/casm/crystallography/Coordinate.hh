#ifndef COORDINATE_HH
#define COORDINATE_HH

#include <iostream>
#include <cmath>
#include <cassert>

#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/container/LinearAlgebra.hh"
#include "casm/crystallography/Lattice.hh"

namespace CASM {

  class Lattice;
  class SymOp;

  namespace Coordinate_impl {
    class FracCoordinate;
    class FracCoordinateComponent;
    class CartCoordinate;
    class CartCoordinateComponent;

  }

  /** \defgroup Coordinate
   *  \ingroup Crystallography
   *  \brief Relates to coordinates
   *
   *  @{
   */

  /// \brief Represents cartesian and fractional coordinates
  ///
  class Coordinate {
  public:
    typedef Eigen::Vector3d vector_type;
    typedef vector_type::Index size_type;

    /// \brief construct a coordinate describing origin of _home lattice
    static Coordinate origin(const Lattice &_home);

    // NOTE: Coordinate does not have a default constructor
    // e.g: this is not allowed-> Coordinate() : home(nullptr) { is_current[FRAC]=false; is_current[CART]=false;};

    ///Minimal constructor only takes a lattice
    explicit Coordinate(const Lattice &_home) :
      m_home(&_home),
      m_frac_coord(vector_type::Zero()),
      m_cart_coord(vector_type::Zero()),
      m_basis_ind(-1) {
    }

    Coordinate(const vector_type &_vec, const Lattice &_home, COORD_TYPE _mode);

    Coordinate(double _x, double _y, double _z, const Lattice &_home, COORD_TYPE _mode);

    /// \brief Set the fractional coordinate vector
    Coordinate_impl::FracCoordinate frac();

    /// \brief const Access the fractional coordinate vector
    const vector_type &frac() const {
      return m_frac_coord;
    }

    /// \brief user override to force const Access the fractional coordinate vector
    const vector_type &const_frac() const {
      return m_frac_coord;
    }

    /// \brief Set a component of the fractional coordinate vector
    Coordinate_impl::FracCoordinateComponent frac(size_type index);

    /// \brief const Access a component of the fractional coordinate vector
    const double &frac(size_type index) const {
      return m_frac_coord(index);
    }

    /// \brief user override to force const Access the fractional coordinate vector
    const double &const_frac(size_type index) const {
      return m_frac_coord(index);
    }

    /// \brief Set Cartesian coordinate vector and update fractional coordinate vector
    Coordinate_impl::CartCoordinate cart();

    /// \brief const Access the Cartesian coordinate vector
    const vector_type &cart() const {
      return m_cart_coord;
    }

    /// \brief user override to force const Access the Cartesian coordinate vector
    const vector_type &const_cart() const {
      return m_cart_coord;
    }

    /// \brief Set a component of the Cartesian coordinate vector
    Coordinate_impl::CartCoordinateComponent cart(size_type index);

    /// \brief const Access a component of the Cartesian coordinate vector
    const double &cart(size_type index) const {
      return m_cart_coord(index);
    }

    /// \brief const Access a component of the Cartesian coordinate vector
    const double &const_cart(size_type index) const {
      return m_cart_coord(index);
    }

    /// \brief Positive translation of this coordinate by RHS.cart()
    Coordinate &operator+=(const Coordinate &RHS);

    /// \brief Negative translation of this coordinate by RHS.cart()
    Coordinate &operator-=(const Coordinate &RHS);

    /// \brief unary minus of this coordinate
    Coordinate operator-() const;

    bool operator==(const Coordinate &RHS) const;

    bool operator!=(const Coordinate &RHS) const {
      return !(*this == RHS);
    }

    bool almost_equal(const Coordinate &RHS, double tol) const;

    /// Returns true if this->min_dist(RHS)<tol
    bool compare(const Coordinate &RHS, double tol = TOL) const;

    /// Returns true if this->min_dist(RHS)<tol
    /// if true, calculates @param translation such that
    /// *this = (RHS+translation)
    bool compare(const Coordinate &RHS, Coordinate &translation, double tol = TOL) const;

    /// Return true -- Exists to allow duck-typing with Site
    bool compare_type(const Coordinate &RHS)const;

    /// Map coordinate into the unit cell using a lattice translation
    /// returns true if *this was already within the unit cell
    bool within();

    ///Same as within(), but lattice translation is stored in Coordinate translation, such that
    /// coord_after = coord_before + translation
    /// returns true if *this was already within the unit cell
    bool within(Coordinate &translation);

    /// Checks to see if coordinate is in the unit cell, but does not translate it
    bool is_within() const;

    /// Number of periodic images of this coordinate that are on voronoi cell boundary
    /// of its home lattice
    int voronoi_number() const;

    /// Number of periodic images of this coordinate that are on voronoi cell boundary
    /// of 'cell'
    int voronoi_number(const Lattice &cell) const;

    ///Map coordinate into the voronoi cell using a lattice translation
    bool voronoi_within();

    ///Same as voronoi_within(), but lattice translation is stored in Coordinate translation such that
    /// coord_after = coord_before + translation
    /// returns true if *this was already within the voronoi cell
    bool voronoi_within(Coordinate &translation);

    ///Checks to see if coordinate is at a lattice translation with respect to the origin
    bool is_lattice_shift(double tol = TOL) const;

    /// \brief  Change the home lattice of the coordinate, selecting one representation
    ///         (either CART or FRAC) that remains invariant
    ///
    /// \param invariant_mode
    /// \parblock
    /// invariant_mode == CART: Cartesian coordinates stay the same, and fractional coordinates are updated
    ///    Ex: (my_coord.set_lattice(superlattice, CART); // this is how superlattices get filled.
    ///
    /// invariant_mode == FRAC: Fractional coordinates stay the same, and Cartesian coordinates are updated
    ///    Ex:  you can apply a strain by changing the lattice and keepin FRAC invariant
    ///            (my_coord.set_lattice(strained_lattice, FRAC);
    ///    Ex:  you can apply a rotation by changing the lattice and keeping FRAC invariant
    ///            (my_coord.set_lattice(rotated_lattice, FRAC);
    /// \endparblock
    void set_lattice(const Lattice &new_lat, COORD_TYPE mode); //John G, use to specify whether to keep CART or FRAC the same, when assigning a new lattice

    /// \brief Set basis Index
    void set_basis_ind(Index _basis_ind) {
      m_basis_ind = _basis_ind;
    }

    /// \brief Access basis Index
    Index basis_ind() const {
      return m_basis_ind;
    }

    /// \brief Access the home lattice of the coordinate
    const Lattice &home() const {
      assert(m_home && "Coordinate doesn't have valid home lattice");
      return *m_home;
    }

    //term is terminal character, prec is precision, pad is field width - precision  (should be greater than 3)
    void read(std::istream &stream, COORD_TYPE mode);
    void print(std::ostream &stream, COORD_TYPE mode, char term = 0, int prec = 7, int pad = 5) const;
    void print(std::ostream &stream, char term = 0, int prec = 7, int pad = 5) const;

    /// \brief Print normalized vector
    void print_axis(std::ostream &stream, COORD_TYPE mode, char term = 0, int prec = 7, int pad = 5) const;

    /// \brief distance (in Angstr) of neighbor from *this
    double dist(const Coordinate &neighbor) const;

    /// \brief Returns distance (in Angstr) to nearest periodic image of neighbor
    ///
    /// Is unsafe if min_dist is comparable to half a lattice vector in length
    double min_dist(const Coordinate &neighbor) const;

    /// \brief Returns distance (in Angstr) to nearest periodic image of neighbor
    ///
    /// This version calculates the translation such that
    /// (neighbor+translation) is the periodic nearest period image of neighbor
    double min_dist(const Coordinate &neighbor, Coordinate &translation)const;

    /// \brief Returns distance (in Angstr) to nearest periodic image of neighbor
    ///
    /// It is safe in all cases, because it uses the lattice Wigner-Seitz cell to
    /// determine the nearest image, but this makes it slower than min_dist
    double robust_min_dist(const Coordinate &neighbor) const;

    /// \brief Returns distance (in Angstr) to nearest periodic image of neighbor
    ///
    /// This version calculates the translation such that
    /// (neighbor+translation) is the periodic nearest period image of neighbor
    /// It is safe in all cases, because it uses the lattice Wigner-Seitz cell to
    /// determine the nearest image, but this makes it slower than min_dist
    double robust_min_dist(const Coordinate &neighbor, Coordinate &translation)const;

    ///Finds same shift as min_dist but returns shift(CART).transpose()*metric*shift(CART)
    double min_dist2(const Coordinate &neighbor, const Eigen::Ref<const Eigen::Matrix3d> &metric) const;

    ///Transform coordinate by symmetry operation (including translation)
    Coordinate &apply_sym(const SymOp &op); //AAB

    ///Transform coordinate by symmetry operation (without translation)
    Coordinate &apply_sym_no_trans(const SymOp &op); //AAB

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);

  private:

    void _update_cart() {
      m_cart_coord = home().lat_column_mat() * m_frac_coord;
    }


    void _update_frac() {
      m_frac_coord = home().inv_lat_column_mat() * m_cart_coord;
    }


    void _set_frac(const Eigen::Ref<const vector_type> &f) {
      m_frac_coord = f;
      _update_cart();
    }


    void _set_frac(size_type ind, double val) {
      m_frac_coord[ind] = val;
      _update_cart();
    }


    void _set_cart(const Eigen::Ref<const vector_type> &c) {
      m_cart_coord = c;
      _update_frac();
    }


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

  inline
  Coordinate Coordinate::origin(const Lattice &_home) {
    return Coordinate(_home);
  }

  jsonParser &to_json(const Coordinate &value, jsonParser &json);

  // Lattice must be set already
  void from_json(Coordinate &value, const jsonParser &json);

  Coordinate operator*(const SymOp &LHS, const Coordinate &RHS); //AAB

  inline
  Coordinate operator+(const Coordinate &LHS, const Coordinate &RHS) {
    return Coordinate(LHS) += RHS;
  }

  inline
  Coordinate operator-(const Coordinate &LHS, const Coordinate &RHS) {
    return Coordinate(LHS) -= RHS;
  }

  /** @} */

  namespace Coordinate_impl {

    /// \brief A class to enable vector assignment to the fractional vector of a Coordinate
    ///
    /// Typically only used indirectly as a temporary when performing
    /// \code
    /// Coordinate coord;
    /// coord.frac() = Coordinate::vector_type(a,b,c);
    /// \endcode
    ///
    /// \relates Coordinate
    ///
    class FracCoordinate {
    public:

      explicit FracCoordinate(Coordinate &coord) :
        m_coord(&coord) {}


      FracCoordinate &operator=(const Eigen::Ref<const Coordinate::vector_type> &vec) {
        m_coord->_set_frac(vec);
        return *this;
      }

      FracCoordinate &operator=(const FracCoordinate &RHS) {
        m_coord->_set_frac(RHS.m_coord->const_frac());
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
    /// \relates Coordinate
    ///
    class FracCoordinateComponent {
    public:

      explicit FracCoordinateComponent(Coordinate &coord, Coordinate::size_type index) :
        m_coord(&coord), m_index(index) {}

      FracCoordinateComponent &operator=(double val) {
        m_coord->_set_frac(m_index, val);
        return *this;
      }

      FracCoordinateComponent &operator=(const FracCoordinateComponent &RHS) {
        m_coord->_set_frac(m_index, RHS.m_coord->const_frac(m_index));
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
    /// \relates Coordinate
    ///
    class CartCoordinate {
    public:

      explicit CartCoordinate(Coordinate &coord) :
        m_coord(&coord) {}

      CartCoordinate &operator=(const Eigen::Ref<const Coordinate::vector_type> &vec) {
        m_coord->_set_cart(vec);
        return *this;
      }

      CartCoordinate &operator=(const CartCoordinate &RHS) {
        m_coord->_set_cart(RHS.m_coord->const_cart());
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
    /// \relates Coordinate
    ///
    class CartCoordinateComponent {
    public:

      explicit CartCoordinateComponent(Coordinate &coord, Coordinate::size_type index) :
        m_coord(&coord), m_index(index) {}

      CartCoordinateComponent &operator=(double val) {
        m_coord->_set_cart(m_index, val);
        return *this;
      }

      CartCoordinateComponent &operator=(const CartCoordinateComponent &RHS) {
        m_coord->_set_cart(m_index, RHS.m_coord->const_cart(m_index));
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


  /// \brief Set the fractional coordinate vector
  inline
  Coordinate_impl::FracCoordinate Coordinate::frac() {
    return Coordinate_impl::FracCoordinate(*this);
  }

  /// \brief Set a component of the fractional coordinate vector
  inline
  Coordinate_impl::FracCoordinateComponent Coordinate::frac(Coordinate::size_type index) {
    return Coordinate_impl::FracCoordinateComponent(*this, index);
  }

  /// \brief Set Cartesian coordinate vector and update fractional coordinate vector
  inline
  Coordinate_impl::CartCoordinate Coordinate::cart() {
    return Coordinate_impl::CartCoordinate(*this);
  }

  /// \brief Set a component of the Cartesian coordinate vector
  inline
  Coordinate_impl::CartCoordinateComponent Coordinate::cart(Coordinate::size_type index) {
    return Coordinate_impl::CartCoordinateComponent(*this, index);
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
