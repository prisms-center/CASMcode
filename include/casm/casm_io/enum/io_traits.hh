#ifndef CASM_support_enum_io_traits
#define CASM_support_enum_io_traits

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "casm/casm_io/Help.hh"

namespace CASM {

class jsonParser;
template <typename T>
struct traits;

/// \defgroup casmIO IO
///
/// \brief Input/output classes and functions

/// \defgroup EnumIO
///
/// \brief Helpers for enum class IO
///
/// \ingroup casmIO

template <typename ENUM>
std::string to_string(ENUM val);

/// \brief Print help message describing recognized strings for allowed enum
/// values
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
template <typename ENUM>
std::string multiline_enum_help() {
  std::stringstream ss;
  ss << "Options are:\n";
  for (auto it = traits<ENUM>::strval.begin(); it != traits<ENUM>::strval.end();
       ++it) {
    ss << "  ";
    for (auto sit = it->second.begin(); sit != it->second.end(); sit++) {
      if (sit != it->second.begin()) {
        ss << " or " << *sit;
      } else {
        ss << *sit;
      }
    }
    ss << "\n";
  }
  return ss.str();
}

/// \brief Print short help message describing recognized strings for allowed
/// enum values
///
/// \param options List of options
/// \param _default If not empty, the option to indicate as default
///
/// Of form:
/// \code
/// Options are: {<other>, 'FRAC' (default), 'CART', 'INTEGRAL'}
/// \endcode
/// or (if other.empty()):
/// \code
/// Options are: {'FRAC' (default), 'CART', 'INTEGRAL'}
/// \endcode
///
/// \ingroup EnumIO
///
template <typename StringContainer>
std::string standard_singleline_help(StringContainer options,
                                     std::string _default = "") {
  std::stringstream ss;
  ss << "Options are: {";
  for (auto it = options.begin(); it != options.end(); ++it) {
    if (it != options.begin()) {
      ss << ", ";
    }
    ss << "'" << *it << "'";
    if (*it == _default) {
      ss << " (default)";
    }
  }
  ss << "}";
  return ss.str();
}

/// \brief Print short help message describing recognized strings for allowed
/// enum values
///
/// Of form:
/// \code
/// Options are: {<other>, 'FRAC' (default), 'CART', 'INTEGRAL'}
/// \endcode
/// or (if other.empty()):
/// \code
/// Options are: {'FRAC' (default), 'CART', 'INTEGRAL'}
/// \endcode
///
/// \ingroup EnumIO
///
template <typename ENUM>
std::string standard_singleline_enum_help(std::string _default = "",
                                          std::string other = "") {
  std::vector<std::string> options;
  if (!other.empty()) {
    options.push_back(std::string("<") + other + ">");
  }
  for (auto it = traits<ENUM>::strval.begin(); it != traits<ENUM>::strval.end();
       ++it) {
    options.push_back(to_string(it->first));
  }
  return standard_singleline_help(options, _default);
}

/// \brief Print short help message describing recognized strings for allowed
/// enum values
///
/// Of form:
/// \code
/// Options are: {'FRAC' (default), 'CART', 'INTEGRAL'}
/// \endcode
///
/// \ingroup EnumIO
///
template <typename ENUM>
std::string singleline_enum_help() {
  return standard_singleline_enum_help<ENUM>(
      traits<ENUM>::strval.begin()->second[0]);
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
template <typename ENUM>
void invalid_enum_string(std::string val) {
  std::stringstream s;
  s << "Invalid " << traits<ENUM>::name << ": " << val << ". "
    << singleline_help<ENUM>();
  // serr << s.str() << "\n";
  // serr << help<ENUM>();
  throw std::invalid_argument(std::string("ERROR: ") + s.str());
}

/// \brief Return string representation of enum class
///
/// \ingroup EnumIO
///
template <typename ENUM>
std::string to_string(ENUM val) {
  return traits<ENUM>::strval.find(val)->second[0];
}

/// \brief Return all matching enum class members from string representation
///
/// Matches using substr of length val.size(), so "FRAC" matches for "F", "FR",
/// "FRA", "FRAC", "FRAC.*".
///
/// \ingroup EnumIO
///
template <typename ENUM>
std::set<ENUM> matches(const std::string &val) {
  std::set<ENUM> res;
  for (auto it = traits<ENUM>::strval.begin(); it != traits<ENUM>::strval.end();
       ++it) {
    for (auto sit = it->second.begin(); sit != it->second.end(); sit++) {
      if (sit->substr(0, val.size()) == val) {
        res.insert(it->first);
      }
    }
  }
  return res;
}

/// \brief Return enum class object from string representation
///
/// \ingroup EnumIO
///
template <typename ENUM>
ENUM from_string(const std::string &val) {
  std::set<ENUM> _matches = matches<ENUM>(val);

  if (_matches.size() == 1) {
    return *_matches.begin();
  }

  invalid_enum_string<ENUM>(val);              // throws
  return traits<ENUM>::strval.begin()->first;  // never reached
}

#define ENUM_TRAITS(ENUM)                                               \
  template <>                                                           \
  struct traits<ENUM> {                                                 \
    static const std::string name;                                      \
                                                                        \
    static const std::multimap<ENUM, std::vector<std::string> > strval; \
  };                                                                    \
  template <>                                                           \
  inline std::string singleline_help<ENUM>() {                          \
    return singleline_enum_help<ENUM>();                                \
  }                                                                     \
  template <>                                                           \
  inline std::string multiline_help<ENUM>() {                           \
    return multiline_enum_help<ENUM>();                                 \
  }

}  // namespace CASM

#endif
