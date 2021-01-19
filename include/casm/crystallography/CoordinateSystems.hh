#ifndef COORDINATESYSTEMS_HH
#define COORDINATESYSTEMS_HH

#include <iostream>

#include "casm/global/enum.hh"

namespace CASM {
namespace xtal {

/** \ingroup Coordinate
 *  @{
 */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// COORD_MODE specifies the current coordinate mode (Fractional or Cartesian)

/// COORD_MODE is a container for the global coordinate system mode
/// i.e., fractional or cartesian.  The global coordinate mode is
/// contained in ACTIVE_MODE.  In order for ACTIVE_MODE to be initialized
/// properly, the following line must be added before 'int main()':
///  COORD_TYPE COORD_MODE::ACTIVE_MODE = FRAC;

///@code
/// COORD_MODE new_mode(FRAC);     //creates new_mode and sets current mode to
/// FRAC Coordinate my_coord=something; my_coord();                    //
/// accesses fractional vector of my_coord do_something_with_frac(); { // Curly
/// brackets limit the scope of anything created inside
///
///    COORD_MODE new_new_mode(COORD_MODE::CART)
///    new_new_mode.check();  //Returns CART
///    COORD_MODE::CHECK();   //Returns CART
///    do_something_with_cart();
/// } //new_new_mode is destroyed here
/// COORD_MODE::CHECK();  //Returns FRAC
///@endcode

class COORD_MODE {
 private:
  /// ACTIVE_MODE is a hidden global variable that specifies what coordinate
  /// system (CART or FRAC)
  /// is currently in use
  static COORD_TYPE ACTIVE_MODE;

  /// old_mode specifies the value of ACTIVE_MODE when this COORD_MODE object
  /// was instantiated
  COORD_TYPE old_mode;

 public:
  /// Static method to check if mode is CART (call using COORD_MODE::IS_CART() )
  static bool IS_CART() { return ACTIVE_MODE == CART; };
  /// Static method to check if mode is FRAC (call using COORD_MODE::IS_FRAC() )
  static bool IS_FRAC() { return ACTIVE_MODE == FRAC; };

  /// get the current mode (call using COORD_MODE::CHECK())
  static COORD_TYPE CHECK() { return ACTIVE_MODE; };

  /// get a string with the name of the active mode
  static std::string NAME() { return NAME(ACTIVE_MODE); };
  static std::string NAME(COORD_TYPE mode) {
    if (mode == FRAC) return "Direct";
    if (mode == CART) return "Cartesian";
    return NAME(ACTIVE_MODE);
  };

  COORD_MODE(COORD_TYPE new_mode) {
    old_mode = ACTIVE_MODE;
    ACTIVE_MODE = new_mode;
  };
  ~COORD_MODE() { ACTIVE_MODE = old_mode; };

  void set(const COORD_TYPE new_mode) { ACTIVE_MODE = new_mode; };
  COORD_TYPE check() { return ACTIVE_MODE; };

  // static bool SelfTest();
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class PERIODICITY_MODE {
 private:
  static PERIODICITY_TYPE ACTIVE_MODE;
  PERIODICITY_TYPE old_mode;

 public:
  //    static bool IS_PERIODIC() { return (ACTIVE_MODE == PERIODIC); }
  // TODO: Revert back
  //    static bool IS_LOCAL() { return (ACTIVE_MODE == LOCAL); }
  static bool IS_LOCAL() { return (ACTIVE_MODE == LOCAL); }
  static bool IS_PERIODIC() { return (ACTIVE_MODE == PERIODIC); }
  static PERIODICITY_TYPE CHECK() { return ACTIVE_MODE; }

  PERIODICITY_MODE(PERIODICITY_TYPE new_mode) {
    old_mode = ACTIVE_MODE;
    ACTIVE_MODE = new_mode;
  };
  ~PERIODICITY_MODE() { ACTIVE_MODE = old_mode; };

  void set(const PERIODICITY_TYPE new_mode) { ACTIVE_MODE = new_mode; };
  PERIODICITY_TYPE check() { return ACTIVE_MODE; };
};

/** @} */

}  // namespace xtal
};  // namespace CASM
#endif
