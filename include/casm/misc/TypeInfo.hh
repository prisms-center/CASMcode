#ifndef CASM_TypeInfo
#define CASM_TypeInfo

#include <cxxabi.h>

#include <cstdlib>
#include <memory>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>

namespace CASM {

/// Get type name as string, via type_name<T>()
///
/// - Thanks to: https://stackoverflow.com/a/20170989
template <typename T>
std::string type_name() {
  typedef typename std::remove_reference<T>::type TR;
  std::unique_ptr<char, void (*)(void *)> name(
      abi::__cxa_demangle(typeid(TR).name(), nullptr, nullptr, nullptr),
      std::free);
  std::string res = name.get();
  if (std::is_const<TR>::value) {
    res += " const";
  }
  if (std::is_volatile<TR>::value) {
    res += " volatile";
  }
  if (std::is_lvalue_reference<T>::value) {
    res += "&";
  } else if (std::is_rvalue_reference<T>::value) {
    res += "&&";
  }
  return res;
}
}  // namespace CASM

#endif
