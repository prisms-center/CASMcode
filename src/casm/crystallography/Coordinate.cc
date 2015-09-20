#include "casm/crystallography/Coordinate.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  //********************************************************************
  /**
   *
   */
  //********************************************************************

  Coordinate &Coordinate::operator +=(const Coordinate &RHS) {
    if(home == RHS.home && calc(FRAC)) {
      coord[FRAC] += RHS(FRAC);
      is_current[CART] = false;
    }
    else if(calc(CART)) {
      coord[CART] += RHS(CART);
      is_current[FRAC] = false;
    }

    else {
      std::cerr << "WARNING: Attempting to add to a Coordinate that has been initialized improperly!\n";
      assert(false);
    }
    return *this;
  }

  Coordinate &Coordinate::operator -=(const Coordinate &RHS) {
    if(home == RHS.home && calc(FRAC)) {
      coord[FRAC] -= RHS(FRAC);
      is_current[CART] = false;
    }
    else if(calc(CART)) {
      coord[CART] -= RHS(CART);
      is_current[FRAC] = false;
    }
    else {
      std::cerr << "WARNING: Attempting to subtract from a Coordinate that has been initialized improperly!\n";
      assert(false);
    }
    return *this;

  }

  //********************************************************************
  /**
   *
   */
  //********************************************************************


  Coordinate Coordinate::operator +(const Coordinate &RHS) const {
    Coordinate tcoord(*this);
    return tcoord += RHS;
  }

  //********************************************************************
  /**
   *
   */
  //********************************************************************

  Coordinate Coordinate::operator -(const Coordinate &RHS) const {
    Coordinate tcoord(*this);
    return tcoord -= RHS;
  }

  //********************************************************************
  /**
   *
   */
  //********************************************************************

  Coordinate Coordinate::operator-() const {
    if(is_current[CART]) return Coordinate(-coord[CART], *home, CART);
    return Coordinate(-coord[FRAC], *home, FRAC);
  }

  //********************************************************************
  /**
   *
   */
  //********************************************************************

  bool Coordinate::operator ==(const Coordinate &RHS) const {
    return (*this - RHS)().is_zero();
  }

  //********************************************************************

  bool Coordinate::unsafe_compare(const Coordinate &RHS, COORD_TYPE mode) const {
    return (std::abs(coord[mode][0] - RHS.coord[mode][0]) < TOL &&
            std::abs(coord[mode][1] - RHS.coord[mode][1]) < TOL &&
            std::abs(coord[mode][2] - RHS.coord[mode][2]) < TOL);
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
  /*
  Coordinate operator*(const Coordinate &RHS, int m) {
    Coordinate tcoord(RHS);
    for(int i = 0; i < 3 ; i++) {
      tcoord[i] = RHS[i] * m;
    }
    return tcoord;
  }
  */

  //********************************************************************

  void Coordinate::read(std::istream &stream) {
    stream >> coord[mode_ind()];
    is_current[mode_ind()] = true;
    is_current[!mode_ind()] = false;
    return;
  }

  //********************************************************************
  void Coordinate::read(std::istream &stream, COORD_TYPE mode) {
    stream >> coord[mode];
    is_current[mode] = true;
    is_current[!mode] = false;
    return;
  }

  //********************************************************************
  void Coordinate::print(std::ostream &stream, COORD_TYPE mode, char term, int prec, int pad) const {

    stream.precision(prec);
    stream.width(prec + pad);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

    if(mode != COORD_DEFAULT)
      stream << (*this)(mode);
    else
      stream << (*this)();
    if(term) stream << term;
    return;
  }

  //********************************************************************

  void Coordinate::print(std::ostream &stream, char term, int prec, int pad) const {

    stream.precision(prec);
    stream.width(prec + pad);
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    stream << (*this)();
    if(term) stream << term;
    return;
  }


  //********************************************************************
  /**
   * Finds distance between two coordinates
   */
  //********************************************************************

  double Coordinate::dist(const Coordinate &neighbor) const {

    return ((*this)(CART) - neighbor(CART)).length();

  }


  //********************************************************************
  /**
   * Finds minimum distance from any periodic image of a coordinate to any
   * periodic image of a neighboring coordinate
   */
  //********************************************************************

  double Coordinate::min_dist(const Coordinate &neighbor) const {
    Coordinate tcoord(*this);
    tcoord -= neighbor;

    if(!tcoord.calc(FRAC)) {
      std::cerr << "Attempting to find minimum distance between two points that have been initialized inproperly.\n";
      return NAN;
    }

    for(int i = 0; i < 3; i++)
      tcoord.at(i, FRAC) -= round(tcoord.get(i, FRAC));

    return tcoord(CART).length();

  }

  //********************************************************************
  /**
   * Finds minimum distance from any periodic image of a coordinate to any
   * periodic image of a neighboring coordinate.  Also passes the
   * calculated shift in the argument.
   */
  //********************************************************************


  double Coordinate::min_dist(const Coordinate &neighbor, Coordinate &shift) const {
    Coordinate tcoord(*this);
    tcoord -= neighbor;

    if(!tcoord.calc(FRAC)) {
      std::cerr << "Attempting to find minimum distance between two points that have been initialized inproperly.\n";
      return NAN;
    }

    for(int i = 0; i < 3; i++)
      tcoord.at(i, FRAC) -= round(tcoord.get(i, FRAC));

    shift = tcoord;
    return tcoord(CART).length();


  };

  //********************************************************************
  /**
   * Applies symmetry to a coordinate.
   */
  //********************************************************************

  Coordinate &Coordinate::apply_sym(const SymOp &op) {
    if(calc()) {
      coord[mode_ind()] = op.get_matrix() * coord[mode_ind()] + op.tau()();
      is_current[!mode_ind()] = false;
    }
    else {
      std::cerr << "WARNING: Failed to apply symmetry to Coordinate" << std::endl;
      assert(false);
    }
    return *this;
  }

  //********************************************************************
  /**
   * Applies symmetry to a coordinate, but doesn't apply translation
   */
  //********************************************************************

  Coordinate &Coordinate::apply_sym_no_trans(const SymOp &op) {
    if(calc()) {
      coord[mode_ind()] = op.get_matrix() * coord[mode_ind()];
      is_current[!mode_ind()] = false;
    }
    else
      std::cerr << "WARNING: Failed to apply symmetry to Coordinate" << std::endl;

    return *this;
  }

  //********************************************************************
  /**
   * Calculates the coordinate of the current mode.
   *
   * If the coordinate of the current mode has not been calculated, and
   * the cooordinate of the other mode and the coordinate transformation
   * matrix both exist, then the coordinates are calculated and the
   * corresponding is_current is set to true.  On the other hand, if the
   * indicated mode has been calculated or if the coordinate of the other
   * mode does not exist, or the coordinate transformation matrix does not
   * exist, then the function simlpy returns the is_current boolean value.
   */
  //********************************************************************

  //inline
  bool Coordinate::calc() const {
    if(!is_current[mode_ind()] && is_current[!mode_ind()] && home) {
      coord[mode_ind()] = (home->coord_trans(!mode_ind())) * coord[!mode_ind()];
      is_current[mode_ind()] = true;
    }

    return is_current[mode_ind()];
  }


  //********************************************************************
  /**
   * Calculates the coordinates of the specified mode
   *
   * Returns true or false depending on whether the coordinates of the
   * specified mode have been calculated.
   * @param COORD_TYPE mode
   */
  //********************************************************************

  //inline
  bool Coordinate::calc(COORD_TYPE mode) const {
    if(!is_current[mode] && is_current[!mode] && home) {
      coord[mode] = (home->coord_trans(!mode)) * coord[!mode];
      //Have to acces home's members using '->' since it's a pointer
      is_current[mode] = true;
    }

    return is_current[mode];
  }

  //********************************************************************
  Coordinate Coordinate::get_normal_vector(Coordinate coord_2, Coordinate coord_3) {
    if(home != coord_2.get_home() || home != coord_3.get_home()) {
      std::cout << "WARNING: Your coordinates don't live in the same home!" << std::endl;
      std::cout << "See Coordinate::get_normal_vector" << std::endl << std::endl;
    }
    Vector3<double> normal;
    Vector3< double > tvec_1;
    Vector3< double > tvec_2;

    tvec_1 = coord_2(CART) - (*this)(CART);
    tvec_2 = coord_3(CART) - (*this)(CART);
    normal = tvec_1.cross(tvec_2);
    normal.normalize();

    return Coordinate(normal, *home, CART);
  }



  //********************************************************************
  /**
   * Assigns lattice associated with Coordinate.
   *
   * @param Lattice new_lat
   */
  //********************************************************************

  //inline
  void Coordinate::set_lattice(const Lattice &new_lat) {
    if(home == &new_lat)
      return;

    calc(CART);
    home = &new_lat;
    is_current[FRAC] = false;
    return;
  }

  //********************************************************************
  //John G. You decide which coordinates (FRAC or CART) to keep the same when you set a new lattice.
  //CART: Adds empty space around atoms (This is what previous set_lattice does)
  //FRAC: Shear atoms with vectors
  //inline
  void Coordinate::set_lattice(const Lattice &new_lat, COORD_TYPE mode) {
    if(home == &new_lat)
      return;

    calc(mode);
    home = &new_lat;
    is_current[(mode + 1) % 2] = false;
    return;
  }

  //********************************************************************
  /**
   * Returns the value of a Coordinate specified by index ind.
   *
   * @param index ind
   * @return coordinate values corresponding to index
   */
  //********************************************************************

  //inline
  double &Coordinate::operator[](int ind) {

    calc();
    is_current[mode_ind()] = true;
    is_current[!mode_ind()] = false;
    return coord[mode_ind()][ind];

  }

  //********************************************************************
  /**
   * Retrieves one of the coordinate values specified by index of the current mode.
   * @param integer index
   * @return coordinate value of index
   * @return NAN
   */
  //********************************************************************
  //inline
  double Coordinate::get(int ind) const {
    if(calc())
      return coord[mode_ind()][ind];
    else
      return NAN;
  }

  //********************************************************************
  /**
   * Retrieves one of the coordinate values specified by index and mode.
   * @param integer index
   * @return coordinate value of index
   * @return NAN
   */
  //********************************************************************

  //inline
  double Coordinate::get(int ind, COORD_TYPE mode) const {
    if(calc(mode))
      return coord[mode][ind];
    else
      return NAN;
  }

  //********************************************************************
  /**
   * Retrieves one of the coordinate values specified by index of the current mode.
   * @param integer index
   * @return coordinate value of index
   * @return NAN
   */
  //********************************************************************

  //inline
  double &Coordinate::at(int ind) {
    calc();
    is_current[mode_ind()] = true;
    is_current[!mode_ind()] = false;
    return coord[mode_ind()][ind];
  }

  //********************************************************************
  /**
   * Retrieves one of the coordinate values specified by index and mode.
   * @param integer index
   * @return coordinate value of index
   * @return NAN
   */
  //********************************************************************

  //inline
  double &Coordinate::at(int ind, COORD_TYPE mode) {
    calc(mode);
    is_current[mode] = true;
    is_current[!mode] = false;
    return coord[mode][ind];
  }

  //********************************************************************
  /**
   * Overloads () operator to allow casting as Vector3.
   *
   * Using () operator, casts the Coordinate as Vector3 in its current
   * coordinate mode.
   * @return coordinate in vector form
   */
  //********************************************************************
  //inline
  Vector3< double > &Coordinate::operator()() {
    if(!calc()) {
      std::cerr << "WARNING: Accessing uninitialized coordinate!" << std::endl;
      assert(0);
    }
    is_current[mode_ind()] = true;
    is_current[!mode_ind()] = false;
    return coord[mode_ind()];
  }

  //inline
  const Vector3< double > &Coordinate::operator()() const {
    if(!calc()) {
      std::cerr << "WARNING: Accessing uninitialized coordinate!" << std::endl;
      assert(0);
    }
    return coord[mode_ind()];
  }

  //********************************************************************
  /**
   * Overloads () operator to allow casting as Vector3 for specified mode.
   *
   * Using () operator, casts the Coordinate as Vector3 in its current
   * coordinate mode.
   * @return coordinate in vector form
   */
  //********************************************************************

  //inline
  Vector3< double > &Coordinate::operator()(COORD_TYPE mode) {

    if(!calc(mode) && !home) {
      std::cerr << "WARNING:  Performing lookup of undefined coordinate." << std::endl;
      assert(0);
    }
    is_current[mode] = true;
    is_current[!mode] = false;
    return coord[mode];
  }

  //inline
  const Vector3< double > &Coordinate::operator()(COORD_TYPE mode) const {
    if(!calc(mode)) {
      std::cerr << "WARNING:  Performing lookup of undefined coordinate." << std::endl;
      assert(0);
    }
    return coord[mode];
  }
  //********************************************************************
  /**
   * Overloads () operator to allow casting as Vector3.
   *
   * Casts as Vector3 in the current coordinate mode.
   * @return coordinate in vector form
   */
  //********************************************************************

  //inline
  Coordinate::operator Vector3< double >() {
    if(!calc()) {
      std::cerr << "WARNING: Accessing uninitialized Coordinate!" << std::endl;
      assert(0);
    }
    is_current[mode_ind()] = true;
    is_current[!mode_ind()] = false;
    return coord[mode_ind()];
  }

  //********************************************************************
  /**
   * Updates both cartesian and fractional coordinates
   *
   * Calls the calc(COORD_TYPE mode) function for both the fractional
   * and cartesian modes.
   */
  //********************************************************************

  //inline
  bool Coordinate::update() const {
    return calc(FRAC) && calc(CART);
  }

  //********************************************************************
  /**
   * Updates specified mode
   *
   * Refreshes value, regardless of whether it is already current
   */
  //********************************************************************

  //inline
  bool Coordinate::update(COORD_TYPE mode) const {
    if(is_current[!mode] && home) {
      coord[mode] = (home->coord_trans(!mode)) * coord[!mode];
      //Have to acces home's members using '->' since it's a pointer
      return is_current[mode] = true;
    }
    return false;

  }


  //********************************************************************
  /**
   * erases value of specifies mode
   *
   * First calculate value in !mode, and then set is_current[mode] to false
   * returns true if at conclusion, !mode has valid value
   */
  //********************************************************************

  //inline
  bool Coordinate::invalidate(COORD_TYPE mode) {
    if(is_current[!mode]) {
      is_current[mode] = false;
      return true;
    }
    update();
    is_current[mode] = false;
    return is_current[!mode];
  }

  //********************************************************************
  /**
   * Within: if a point is outiside the cell defined by Lattice home
   * then map that point back into the cell via translation by lattice vectors
   *
   */
  //********************************************************************
  bool Coordinate::within() {
    if(!calc(FRAC)) {
      std::cerr << " Coordinate::within() called on a coordinate that has no valid fractional coordinate." << std::endl;
      std::cerr << "Frac: " << std::setw(10) <<  coord[FRAC] << " " << is_current[FRAC] << "\n Cart: " << std::setw(10) << coord[CART] << " " << is_current[CART] << '\n' << "Lattice: " << home;

      exit(1);
    }
    if(PERIODICITY_MODE::IS_LOCAL()) return true;

    bool is_within = true;
    double tshift;
    for(int i = 0; i < 3; i++) {
      tshift = floor(coord[FRAC][i] + 1E-6);
      if(std::abs(tshift) > TOL) {
        is_within = false;
        coord[FRAC][i] -= tshift;
        is_current[CART] = false;
      }
    }
    return is_within;
  };
  //***********************************

  bool Coordinate::is_within() const {
    if(!calc(FRAC)) {
      std::cerr << " Coordinate::within() called on a coordinate that has no valid fractional coordinate." << std::endl;
      std::cerr << "Frac: " << std::setw(10) <<  coord[FRAC] << " " << is_current[FRAC] << "\n Cart: " << std::setw(10) << coord[CART] << " " << is_current[CART] << '\n' << "Lattice: " << home;

      exit(1);
    }

    double tshift;
    for(int i = 0; i < 3; i++) {
      tshift = floor(coord[FRAC][i] + 1E-6);
      if(std::abs(tshift) > TOL) {
        return false;
      }
    }
    return true;
  };



  //********************************************************************


  //********************************************************************

  bool Coordinate::within(Coordinate &translation) {

    translation.set_lattice(*home);
    if(!calc(FRAC)) {
      std::cerr << " Coordinate::within called on a coordinate that has no valid fractional coordinate." << std::endl;
      exit(1);
    }
    if(PERIODICITY_MODE::IS_LOCAL()) return true;

    bool is_within = true;
    double tshift;
    for(int i = 0; i < 3; i++) {
      tshift = floor(coord[FRAC][i] + 1E-6);
      if(std::abs(tshift) > TOL) {
        is_within = false;
        coord[FRAC][i] -= tshift;
        is_current[CART] = false;
        translation.at(i, FRAC) = -tshift;
      }
      else
        translation.at(i, FRAC) = 0.0;
    }
    return is_within;

  };


  //********************************************************************
  /**
   * Checks to see if coordinate is within voronoi cell of its home lattice.
   */
  //********************************************************************

  int Coordinate::voronoi_number() const {
    return voronoi_number(*home);
  }

  //********************************************************************
  /**
   * Checks to see if coordinate is within voronois cell of specified lattice
   */
  //********************************************************************

  int Coordinate::voronoi_number(const Lattice &cell) const {
    if(!calc(CART)) {
      std::cerr << "WARNING: Attempting to find Voronoi number of improperly initialized coordinate! Exiting... \n";
      return 0;
    }
    return cell.voronoi_number(coord[CART]);

  }


  //********************************************************************
  /**
   * Applies Lattice translations until coordinate is within voronois cell of its home lattice.
   * **** HAS NOT BEEN PROPERLY TESTED
   */
  //********************************************************************

  bool Coordinate::voronoi_within() {
    bool was_within(true);

    if(!calc(CART)) {
      std::cerr << "WARNING: Attempting to find Voronoi number of improperly initialized coordinate! Exiting... \n";
      return 0;
    }

    while(voronoi_number() < 0) {
      was_within = false;
      Coordinate tcoord(home->max_voronoi_vector(coord[CART]), *home, CART);

      Vector3<int> int_vec(round(coord[CART].dot(tcoord(CART)) / 2)*tcoord(FRAC).scale_to_int());

      at(0, FRAC) -= int_vec[0];
      at(1, FRAC) -= int_vec[1];
      at(2, FRAC) -= int_vec[2];
    }

    return was_within;
  }

  //********************************************************************
  /**
   * Checks to see if coordinate describes shift by a general lattice vector l*V1+m*V2+n*V3, where l, m, n are integer
   */
  //********************************************************************

  bool Coordinate::is_lattice_shift() {
    if(!calc(FRAC)) {
      //Return true if coordinate describes origin
      if(calc(CART))
        return std::abs(coord[CART][0]) < TOL && std::abs(coord[CART][1]) < TOL && std::abs(coord[CART][2]) < TOL;
    }

    //If mode is local, return true only if coordinate describes origin
    if(PERIODICITY_MODE::IS_LOCAL())
      return std::abs(coord[FRAC][0]) < TOL && std::abs(coord[FRAC][1]) < TOL && std::abs(coord[FRAC][2]) < TOL;

    return (std::abs(coord[FRAC][0] - round(coord[FRAC][0])) < TOL
            && std::abs(coord[FRAC][1] - round(coord[FRAC][1])) < TOL
            && std::abs(coord[FRAC][2] - round(coord[FRAC][2])) < TOL);
  }

  //********************************************************************
  /**
   * Used in Coordinate SelfTest() to test switching of modes.
   */
  //********************************************************************

  bool Coordinate::switch_test() {
    // Default COORD_MODE is fractional so there should be fractional coordinates
    // and no cartesian coordinates
    if(is_frac() && !is_cart() && is_current[FRAC] && is_current[CART] == false) {
      COORD_MODE new_mode(CART);
      calc(CART);

      // Now that COORD_MODE is set to cart and cart has been calculated,
      // test the mode and see if the cart component of is_currents is now true
      if(is_cart() && !is_frac() && is_current[CART] && is_current[FRAC]) {
        COORD_MODE new_mode(FRAC);
        // Changed one of the frac coordinate values
        (*this)[1] = 0.4;

        // Because the value of frac was edited, the original cartesian coordinates are no longer valid
        // so is_current[CART] would be false
        if(is_frac() && !is_cart() && is_current[FRAC] && !is_current[CART]) {
          (*this)[1] = 0.2; // Changed back to original value so subsequent tests would make sense
          return true;
        }
      }
    }

    return false;
  }


  //********************************************************************
  /**
   * Used in Coordinate SelfTest() to test that F is calculated properly
   */
  //********************************************************************

  bool Coordinate::calc_F_test() {
    if((((*this).get(0, FRAC) - 0.5) < TOL)
       && (((*this).get(1, FRAC) - 0.2) < TOL)
       && (((*this).get(2, FRAC) - 0.1) < TOL)) {
      return true;
    }
    else return false;
  }

  //********************************************************************
  /**
   * Used in Coordinate SelfTest() to test that C is calculated properly
   */
  //********************************************************************

  bool Coordinate::calc_C_test() {
    if((((*this).get(0, CART) - 1.0) < TOL)
       && (((*this).get(1, CART) - 0.2) < TOL)
       && (((*this).get(2, CART) - 0.3) < TOL)) {
      return true;
    }
    else return false;

  }

  //********************************************************************
  /**
   * Read/write coordinate to json
   */
  //********************************************************************

  jsonParser &Coordinate::to_json(jsonParser &json) const {
    json.put_obj();

    // mutable Vector3< double > coord[2];
    update();
    json["FRAC"] = coord[FRAC];
    json["CART"] = coord[CART];

    // mutable bool is_current[2];
    json["is_current"] = jsonParser::object();
    json["is_current"]["FRAC"] = is_current[FRAC];
    json["is_current"]["CART"] = is_current[CART];

    // mutable int basis_ind;
    json["basis_ind"] = basis_ind();

    return json;
  }

  void Coordinate::from_json(const jsonParser &json) {
    try {
      // mutable Vector3< double > coord[2];
      CASM::from_json(coord[FRAC], json["FRAC"]);
      CASM::from_json(coord[CART], json["CART"]);

      // mutable bool is_current[2];
      CASM::from_json(is_current[FRAC], json["is_current"]["FRAC"]);
      CASM::from_json(is_current[CART], json["is_current"]["CART"]);

      // mutable int basis_ind;
      CASM::from_json(m_basis_ind, json["basis_ind"]);

    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  };

  jsonParser &to_json(const Coordinate &coord, jsonParser &json) {
    return coord.to_json(json);
  }

  void from_json(Coordinate &coord, const jsonParser &json) {
    try {
      coord.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
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


  //********************************************************************

  /*
  bool Coordinate::SelfTest() {
    //set lattice parameters
    Vector3< double > vec1(2, 0, 0), vec2(0, 1, 0), vec3(0, 0, 3);

    Lattice my_lat(vec1, vec2, vec3);

    Coordinate my_coord(my_lat);
    my_coord.set_lattice(my_lat);
    my_coord[0] = 0.5;
    my_coord[1] = 0.2;
    my_coord[2] = 0.1;


    using namespace SelfTestable;
    TestSet("Coordinate");

    // Test that switches from F to C to F produces corresponding switch
    // in is_currents
    if(!Test(my_coord.switch_test() , "Switch Mode Test")) return false;

    // Test get_F works
    if(!Test(my_coord.calc_F_test() , "Calc Frac Test")) return false;

    // Test get_C works
    if(!Test(my_coord.calc_C_test() , "Calc Cart Test")) return false;

    return true;
  }; // end of SelfTest for Coordinate
  */


};

