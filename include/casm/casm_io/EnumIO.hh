#ifndef CASM_EnumIO
#define CASM_EnumIO

#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  /// \defgroup casmIO IO
  ///
  /// \brief Input/output classes and functions

  /// \defgroup EnumIO
  ///
  /// \brief Helpers for enum class IO
  ///
  /// \ingroup casmIO

  template <typename T>
  struct traits {};

  /// \brief Print help message describing recognized strings for allowed enum values
  ///
  /// Of form:
  /// \code
  /// Options are:
  ///   CART or cart
  ///   FRAC or frac or DIRECT or direct
  ///   INTEGRAL or integral
  /// \endcode
  ///
  /// \ingroup EnumIO
  ///
  template<typename ENUM>
  std::string help() {
    std::stringstream ss;
    ss << "Options are:\n";
    for(auto it = traits<ENUM>::strval.begin(); it != traits<ENUM>::strval.end(); ++it) {
      ss << "  ";
      for(auto sit = it->second.begin(); sit != it->second.end(); sit++) {
        if(sit != it->second.begin()) {
          ss << " or " << *sit;
        }
        else {
          ss << *sit;
        }
      }
      ss << "\n";
    }
    return ss.str();
  }

  /// \brief Throw invalid_argument error for unrecognized strings
  ///
  /// Prints to serr:
  /// \code
  /// Invalid CoordType: mistakevalue
  /// Options are:
  ///   CART or cart
  ///   FRAC or frac or DIRECT or direct
  ///   INTEGRAL or integral
  /// \endcode
  ///
  /// \ingroup EnumIO
  ///
  template<typename ENUM>
  void invalid_enum_string(std::string val, std::ostream &serr) {
    std::stringstream s;
    s << "Invalid " << traits<ENUM>::name << ": " << val;
    serr << s.str() << "\n";
    serr << help<ENUM>();
    throw std::invalid_argument(std::string("ERROR: ") + s.str());
  }

  /// \brief Return string representation of enum class
  ///
  /// \ingroup EnumIO
  ///
  template<typename ENUM>
  std::string to_string(ENUM val) {
    return traits<ENUM>::strval.find(val)->second[0];
  }

  /// \brief Return enum class object from string representation
  ///
  /// \ingroup EnumIO
  ///
  template<typename ENUM>
  ENUM from_string(const std::string &val) {
    for(auto it = traits<ENUM>::strval.begin(); it != traits<ENUM>::strval.end(); ++it) {
      for(auto sit = it->second.begin(); sit != it->second.end(); sit++) {
        if(*sit == val) {
          return it->first;
        }
      }
    }

    invalid_enum_string<ENUM>(val, std::cerr);
    return traits<ENUM>::strval.begin()->first;
  }

#define ENUM_TRAITS(ENUM) \
  template<> \
  struct traits<ENUM> { \
  \
    static const std::string name; \
  \
    static const std::multimap<ENUM, std::vector<std::string> > strval; \
  \
  }; \
 
#define ENUM_IO(ENUM) \
  inline std::ostream &operator<<(std::ostream &sout, const ENUM& val) { \
    sout << to_string<ENUM>(val); \
    return sout; \
  } \
  \
  inline std::istream &operator>>(std::istream &sin, ENUM& val) { \
    std::string s; \
    sin >> s; \
    val = from_string<ENUM>(s); \
    return sin; \
  } \
  \
  inline jsonParser &to_json(const ENUM &val, jsonParser &json) { \
    return to_json(to_string<ENUM>(val), json); \
  } \
  \
  inline void from_json(ENUM& val, const jsonParser& json) { \
    val = from_string<ENUM>(json.get<std::string>()); \
  } \
 

}

#endif
