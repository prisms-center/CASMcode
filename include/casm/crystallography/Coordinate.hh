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


  class Coordinate {
  protected:

    Lattice const *home;

  private:
    //coords[FRAC] is fractional coordinate; coord[CART] is cartesian coordinate
    mutable Vector3< double > coord[2];
    mutable bool is_current[2];


    //Following routines added by Ivy to be used in SelfTest
    bool switch_test();
    bool calc_F_test();
    bool calc_C_test();
    Index m_basis_ind;
  public:

    /**NOTE: Coordinate does not have a default constructor
       e.g: this is not allowed-> Coordinate() : home(NULL) { is_current[FRAC]=false; is_current[CART]=false;}; **/

    ///Minimal constructor only takes a lattice
    explicit Coordinate(const Lattice &init_home) : home(&init_home), m_basis_ind(-1) {
      is_current[FRAC] = false, is_current[CART] = false;
    };

    Coordinate(const Coordinate &init_coord, const Lattice &init_home) : home(&init_home), m_basis_ind(-1) {
      coord[CART] = init_coord(CART);
      is_current[CART] = true;
      is_current[FRAC] = false;
    };

    Coordinate(const Vector3<double> &init_vec, const Lattice &init_home) : home(&init_home), m_basis_ind(-1) {
      coord[mode_ind()] = init_vec;
      is_current[mode_ind()] = true;
      is_current[!mode_ind()] = false;
    };

    Coordinate(const Vector3<double> &init_vec, const Lattice &init_home, COORD_TYPE mode) : home(&init_home), m_basis_ind(-1) {
      coord[mode] = init_vec;
      is_current[mode] = true;
      is_current[!mode] = false;
    };

    Coordinate &operator +=(const Coordinate &RHS); //Ivy
    Coordinate operator +(const Coordinate &RHS) const;
    Coordinate &operator -=(const Coordinate &RHS);
    Coordinate operator -(const Coordinate &RHS) const; //Ivy
    Coordinate operator -() const; //Ivy

    bool operator ==(const Coordinate &RHS) const; //Ivy
    bool operator !=(const Coordinate &RHS) const {
      return !(*this == RHS);
    };

    /** check if two Coordinates are the same without updating CART/FRAC
    use following only if you are certain of initialization status of coordinates **/
    bool unsafe_compare(const Coordinate &RHS, COORD_TYPE mode) const;

    ///These compares exist to make interface consistent with site
    bool compare(const Coordinate &RHS, double tol = TOL) const;
    bool compare(const Coordinate &RHS, Coordinate &shift, double tol = TOL) const;
    bool compare_type(const Coordinate &RHS)const;

    ///Methods to check on the global COORD_MODE state
    bool is_cart() const {
      return COORD_MODE::IS_CART();
    };
    bool is_frac() const {
      return COORD_MODE::IS_FRAC();
    };
    int  mode_ind() const {
      return int(COORD_MODE::CHECK());
    }; //this will index our arrays

    ///calc() makes sure that that coordinate has up-to-date value in current mode, as determined by is_current[mode]
    ///returns true as long as a valid value exists or can be calculated
    bool calc() const;
    bool calc(COORD_TYPE mode) const; //same as calc(), but for specified mode

    //update() calls calc(FRAC) and calc(CART) and returns true if successful
    //update(mode) refreshes the value of the specified mode, regardless of whether it is current
    //returns true unless Coordinate is improperly initialized
    bool update() const;
    bool update(COORD_TYPE mode) const;

    //invalidate(mode) does calc(!mode) and then marks value in specified mode as not current
    bool invalidate(COORD_TYPE mode);

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
    operator Vector3< double >();


    /**Retur vector in the current coordinate mode
       example: Coordinate my_coord(prim);
                Vector3<double> my_cart_coord = my_coord(CART); **/
    Vector3< double > &operator()();  //Cast as vector in the current coordinate mode
    Vector3< double > &operator()(COORD_TYPE mode);

    const Vector3< double > &operator()() const;  //Cast as vector in the current coordinate mode
    const Vector3< double > &operator()(COORD_TYPE mode) const;

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

};
#endif
