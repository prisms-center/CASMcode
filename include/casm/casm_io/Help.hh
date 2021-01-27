#ifndef CASM_Help
#define CASM_Help

#include "casm/misc/TypeInfo.hh"

namespace CASM {

/// Help message for reading a type T from input
template <typename T>
std::string multiline_help() {
  return std::string("No help found for type '") + type_name<T>() + "'";
}

/// Help message for reading a type T from input
template <typename T>
std::string singleline_help() {
  return std::string("No help found for type '") + type_name<T>() + "'";
}

/// Uses 'multiline_help<T>()' by default
template <typename T>
std::string help() {
  return multiline_help<T>();
}
}  // namespace CASM

#endif
